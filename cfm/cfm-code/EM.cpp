/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# em_ft_train.cpp
#
# Description: 	Train the bayesian network parameters using EM applied using
#				the inference algorithms within the LibDai package.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "mpi.h"
#include "EM.h"
#include "IPFP.h"
#include "Comms.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#include "time.h"

static const int MAX_GRAD_ASCENT_ITERATIONS = 1000000;

EM::EM(config_t *a_cfg, FeatureCalculator *an_fc, std::string &a_status_filename){
	cfg = a_cfg; 
	fc = an_fc;
	status_filename = a_status_filename;
	int num_energies_to_include = cfg->map_d_to_energy.back() + 1;
	param = new Param( fc->getFeatureNames(), num_energies_to_include );

	//Initialise the communicator
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );
	if( mpi_rank == MASTER ) comm = new MasterComms();
	else comm = new WorkerComms();

};

EM::~EM(){
  delete param;
  delete comm;
};

void EM::writeStatus( const char *msg ){

	std::ofstream out;
	out.open(status_filename.c_str(), std::fstream::out | std::fstream::app);
	out << msg << std::endl;
	out.close();
}

void EM::writeParamsToFile( std::string &filename ){
	param->saveToFile( filename );
}

double EM::run( std::vector<MolData> &data, int group, bool fullZeroInit ){

	unused_zeroed = 0;
	int iter = 0;

	// Initialise Parameters (randomly) - and share between threads
	if( fullZeroInit ) param->fullZeroInit();
	else param->randomInit();

	comm->broadcastInitialParams( param );

	// EM
	iter = 0;
	double Q, prevQ = -10000000.0;
	while( iter < MAX_EM_ITERATIONS ){

		std::string msg = "EM Iteration " + boost::lexical_cast<std::string>(iter);
		if( comm->isMaster() ) writeStatus( msg.c_str() );
		comm->printToMasterOnly(msg.c_str());

		//Reset sufficient counts
		suft_counts_t suft;
		initSuft(suft, data);

		int num_converged = 0, num_nonconverged =0;
		int tot_numc = 0, total_numnonc = 0;

		time_t before, after;
		before = time( NULL );

		//Do the inference part (E-step)
		std::vector<MolData>::iterator itdata = data.begin();
		for( int molidx = 0; itdata != data.end(); ++itdata, molidx++ ){
		
			if( itdata->getGroup() == group ) continue;
			MolData *moldata = &(*itdata);
			
			//Compute the transition probabilities
			moldata->computeTransitionThetas( *param );
		    moldata->computeTransitionProbabilities();

			//Apply the peak evidence and compute the beliefs
			IPFP ipfp( moldata, cfg); 
			beliefs_t *beliefs = ipfp.calculateBeliefs();
			int status = ipfp.status;
			if( status == NON_CONVERGE || status == OSC_CONVERGE ) num_nonconverged++;
			else if( status == COMPLETE_CONVERGE || status == CONVERGE_AFTER_MOD ) num_converged++;

			//Record the sufficient statistics
			recordSufficientStatistics( suft, molidx, moldata, beliefs);
		}
		after = time( NULL );
		std::string time_msg = "Completed IPFP processing: Time Elapsed = " + boost::lexical_cast<std::string>(after - before) + " seconds";
		comm->printWithWorkerId(time_msg.c_str());
		writeStatus(time_msg.c_str());
		MPI_Barrier(MPI_COMM_WORLD);   	//All threads wait for master

		total_numnonc = comm->collectSumInMaster( num_nonconverged );
		tot_numc = comm->collectSumInMaster( num_converged );
		std::string cvg_msg = "Num Converged: " + boost::lexical_cast<std::string>(tot_numc);
		std::string noncvg_msg = "Num Non-Converged: " + boost::lexical_cast<std::string>(total_numnonc);
		if( comm->isMaster() ) writeStatus( cvg_msg.c_str() );		
		if( comm->isMaster() ) writeStatus( noncvg_msg.c_str() );	
		comm->printToMasterOnly(cvg_msg.c_str());
		comm->printToMasterOnly(noncvg_msg.c_str());

		MPI_Barrier(MPI_COMM_WORLD);   	//All threads wait for master

		//Find a new set of parameters to maximize the expected log likelihood (M-step)
		Q = updateParameters(data, suft);

		//Check for convergence
		double Qdif = fabs( prevQ - Q );
		std::string qdif_str;
		qdif_str += boost::lexical_cast<std::string>(Qdif);
		qdif_str += " ";
		qdif_str += boost::lexical_cast<std::string>(prevQ);
		qdif_str += " ";
		qdif_str += boost::lexical_cast<std::string>(Q);
		writeStatus(qdif_str.c_str());

		if( Q < prevQ - cfg->em_converge_thresh ){ 
			std::string msg = "Warning: Expected likeihood went down!! :(";
			msg += boost::lexical_cast<std::string>(prevQ);
			msg += " ";
			msg += boost::lexical_cast<std::string>(Q);
			writeStatus(msg.c_str());
			std::cout << msg << std::endl;
		}
		prevQ = Q;
		if( Qdif < cfg->em_converge_thresh ){ 
			comm->printToMasterOnly(("EM Converged after " + boost::lexical_cast<std::string>(iter) + " iterations").c_str());
			break;
		}
		iter++;
	}

	if( iter >= MAX_EM_ITERATIONS )
		comm->printToMasterOnly(("Warning: EM did not converge after " + boost::lexical_cast<std::string>(iter) + " iterations.").c_str());

	return Q;		
}

void EM::initSuft(suft_counts_t &suft, std::vector<MolData> &data ){

	//Resize the suft structure for each molecule
	unsigned int num_mols = data.size();
	suft.values.resize(num_mols);
	for( unsigned int i = 0; i < num_mols; i++){
		const FragmentGraph *fg = data[i].getFragmentGraph();
		int len =  fg->getNumTransitions() + fg->getNumFragments();
		int num_spectra = data[i].getNumSpectra();
		suft.values[i].resize( len*num_spectra );
	}
}

void EM::recordSufficientStatistics( suft_counts_t &suft, int molidx, MolData *moldata, beliefs_t *beliefs ){

	const FragmentGraph *fg = moldata->getFragmentGraph();

	unsigned int num_transitions = fg->getNumTransitions();
	unsigned int num_fragments = fg->getNumFragments();

	int len_offset = num_transitions + num_fragments;

	//Accumulate the Sufficient Statistics
	for( unsigned int i = 0; i < num_transitions; i++ ){
		
		const Transition *t = fg->getTransitionAtIdx(i);

		double belief = 0.0;
		int energy = cfg->map_d_to_energy[0];
		if( t->getFromId() == 0 )	//main ion is always id = 0
			belief += exp(beliefs->tn[i][0]);
		for( unsigned int d = 1; d < cfg->model_depth; d++ ){
			energy = cfg->map_d_to_energy[d];
			if( energy != cfg->map_d_to_energy[d-1] ){
				suft.values[molidx][i + cfg->map_d_to_energy[d-1]*len_offset] = belief;
				belief = 0.0;
			}
			belief += exp(beliefs->tn[i][d]);
		}
		suft.values[molidx][i + energy*len_offset] = belief;
	}

	//Accumulate the persistence terms
	int offset = num_transitions;
	for( unsigned int i = 0; i < num_fragments; i++ ){
		
		double belief = 0.0;
		int energy = cfg->map_d_to_energy[0];
		if( i == 0 )	//main ion is always id = 0
			belief += exp(beliefs->ps[i][0]);
		for( unsigned int d = 1; d < cfg->model_depth; d++ ){
			energy = cfg->map_d_to_energy[d];
			if( energy != cfg->map_d_to_energy[d-1] ){
				suft.values[molidx][i + offset + cfg->map_d_to_energy[d-1]*len_offset] = belief;
				belief = 0.0;
			}
			belief += exp(beliefs->ps[i][d]);
		}
		suft.values[molidx][i + offset + energy*len_offset] = belief;
	}
}

double EM::updateParameters( std::vector<MolData> &data, suft_counts_t &suft ){

	unsigned int l;

	//Use Gradient Ascent, starting from the current parameters, to find the new parameters
	double prevQ=100000.0;
	int iter = 0, converged = 0;

	//Initialise gradients
	std::vector<double> grads(param->getNumWeights());
	for( l = 0; l < grads.size(); l++ ) grads[l] = 0.0; 

	//Initial Q and gradient calculation
	double Q = 0.0;
	std::vector<MolData>::iterator itdata = data.begin();
	for( int molidx = 0; itdata != data.end(); ++itdata, molidx++ )
		Q += computeAndAccumulateGradient(grads, molidx, *itdata, param, suft, comm->used_idxs);

	//Collect the used_idxs from all the processors into the MASTER
	comm->setMasterUsedIdxs();

	if( !unused_zeroed && comm->isMaster() ){ 
		zeroUnusedParams();
		unused_zeroed = 1;
	}

	//Add the regularizer
	if( comm->isMaster() ) Q += addRegularizers( grads, param );

	//Accumulate the Q and gradient terms from all the processors in the MASTER
	Q = comm->collectQInMaster( Q );
	comm->printToMasterOnly(("Initial Q=" + boost::lexical_cast<std::string>(Q)).c_str());
	comm->collectGradsInMaster( grads );

	//Main Loop
	double tmp_step = cfg->starting_step_size;
	std::vector<double> tmp_grads(grads.size());
	for( l = 0; l < grads.size(); l++ ) tmp_grads[l] = 0.0;
	while( converged < cfg->converge_count_thresh && iter < MAX_GRAD_ASCENT_ITERATIONS ){

		//Backtracking Line Search
		Param tmp_params(*param);
		double tmpQsum = 0.0, grad_sq = 0.0;
		int line_search_converged = 0;
		if( comm->isMaster() ) grad_sq = computeGradientSquaredTerm( grads );
		
		int search_count = 0;
		double tmpQ, prev_tmpQ = -10000000000000000.0;
		while( !line_search_converged && search_count < cfg->max_search_count ){

			//Reset the Tmp Gradients
			tmpQ = 0.0;
			resetTmpGradients( tmp_grads );
		
			if( comm->isMaster() ) stepAndSetTemporaryParams( tmp_params, grads, tmp_step );
			comm->broadcastParams( &tmp_params );

			//Iterate training data molecules (always batch processing)
			itdata = data.begin();
			for( int molidx = 0; itdata != data.end(); ++itdata, molidx++ ){
				tmpQ += computeAndAccumulateGradient(tmp_grads, molidx, *itdata, &tmp_params, suft, comm->used_idxs);
			}
			if( comm->isMaster() )
				tmpQ += addRegularizers( tmp_grads, &tmp_params );
			
			//Accumulate the Q terms from all the processors in the MASTER
			tmpQ = comm->collectQInMaster( tmpQ );

			if( comm->isMaster() )
				line_search_converged = updateLineSearchConverged(Q, tmpQ, prev_tmpQ, tmp_step, grad_sq );
			line_search_converged = comm->broadcastConverged( line_search_converged );
			
			tmp_step *= cfg->line_search_beta;
			search_count++;
			prev_tmpQ = tmpQ;
		}

		//Update final Q and grads
		Q = tmpQ;
		comm->collectGradsInMaster( tmp_grads );		
		if( comm->isMaster() ){ 
			copyGrads( tmp_grads, grads );
			copyTmpToParamWeights( tmp_params );
		}

		//Check for convergence
		if( comm->isMaster() ){
			if( fabs(prevQ - Q) < cfg->ga_converge_thresh ) converged++;
			else converged = 0;
			std::string msg ="Q = " + boost::lexical_cast<std::string>(Q);
			msg += (" step=" + boost::lexical_cast<std::string>(tmp_step/cfg->line_search_beta));
			msg += (" grad_sq=" + boost::lexical_cast<std::string>(grad_sq));
			std::cout << msg << std::endl;
			writeStatus( msg.c_str());
			prevQ = Q;
		}

		//Move the step back one higher for the next loop (saves spending so long searching
		//for each step if the starting step is too high)
		tmp_step /= (cfg->line_search_beta*cfg->line_search_beta);

		//Master broadcasts converged flag to all
		converged = comm->broadcastConverged( converged );
		iter++;
	}

	//Master computes and broadcasts Qdif to all
	if( comm->isMaster() ){ 
		if( iter < MAX_GRAD_ASCENT_ITERATIONS )
			std::cout << "Gradient Descent converged after " << iter << " iterations." << std::endl;
		else
			std::cout << "Gradient Descent did not converge!" << std::endl;
	}

	//Master broadcasts final Q and param weights to all
	Q = comm->broadcastQ( Q );
	comm->broadcastParams( param );

	return( Q );
}

void EM::stepAndSetTemporaryParams( Param &tmp_params, std::vector<double> &grads_sum, double tmp_step){

	//Make a step - update the (used) weights
	std::set<unsigned int>::iterator it = ((MasterComms *)comm)->master_used_idxs.begin();
	for( ; it != ((MasterComms *)comm)->master_used_idxs.end(); ++it )
		tmp_params.setWeightAtIdx( param->getWeightAtIdx(*it) + grads_sum[*it]*tmp_step, *it );	
}

double EM::computeAndAccumulateGradient(std::vector<double> &grads, int molidx, MolData &moldata, Param *tmp_params, suft_counts_t &suft, std::set<unsigned int> &used_idxs){

	double Q = 0.0;
	const FragmentGraph *fg = moldata.getFragmentGraph();
	unsigned int num_transitions = fg->getNumTransitions();
	unsigned int num_fragments = fg->getNumFragments();
	
	int offset = num_transitions;

	//Compute the latest transition thetas
	moldata.computeTransitionThetas( *tmp_params );
	suft_t *suft_values = &(suft.values[molidx]);
	
	for( unsigned int energy = 0; energy < tmp_params->getNumEnergyLevels(); energy++ ){

		unsigned int grad_offset = energy * tmp_params->getNumWeightsPerEnergyLevel();
		unsigned int suft_offset = energy * (num_transitions + num_fragments);

		//Iterate over from_id (i)
		tmap_t::const_iterator it = fg->getFromIdTMap()->begin();
		for( int from_idx=0; it != fg->getFromIdTMap()->end(); ++it, from_idx++ ){

			//Calculate the denominator of the sum terms
			double denom = 1.0;
			std::vector<int>::const_iterator itt = it->begin();
			for( ; itt != it->end(); ++itt )
				denom += exp( moldata.getThetaForIdx(energy, *itt) );

			//Complete the innermost sum terms	(sum over j')	
			std::map<unsigned int, double> sum_terms;
			for( itt = it->begin(); itt != it->end(); ++itt ){
				const FeatureVector *fv = moldata.getFeatureVectorForIdx(*itt);
				std::vector<feature_t>::const_iterator fvit = fv->getFeatureBegin();				
				for( ; fvit != fv->getFeatureEnd(); ++fvit ){
					double val = exp( moldata.getThetaForIdx(energy, *itt))/denom;
					if( sum_terms.find(*fvit) != sum_terms.end() )
						sum_terms[*fvit] += val;
					else
						sum_terms[*fvit] = val;
				}
			}

			//Accumulate the transition (i \neq j) terms of the gradient (sum over j)
			double nu_sum = 0.0;
			for( itt = it->begin(); itt != it->end(); ++itt ){
				double nu = (*suft_values)[*itt + suft_offset];
				nu_sum += nu;
				const FeatureVector *fv = moldata.getFeatureVectorForIdx(*itt);
				std::vector<feature_t>::const_iterator fvit = fv->getFeatureBegin();				
				for( ; fvit != fv->getFeatureEnd(); ++fvit ){
					grads[*fvit + grad_offset] += nu;
					used_idxs.insert(*fvit + grad_offset);
				}
				Q += nu*(moldata.getThetaForIdx(energy, *itt) - log(denom));
			}

			//Accumulate the last term of each transition and the 
			//persistence (i = j) terms of the gradient and Q
			std::map<unsigned int, double>::iterator sit = sum_terms.begin();
			double nu = (*suft_values)[offset + from_idx + suft_offset];	//persistence (i=j)
			for( ; sit != sum_terms.end(); ++sit ){ 
				grads[sit->first + grad_offset] -= (nu_sum + nu)*sit->second;
				used_idxs.insert(sit->first + grad_offset);
			}
			Q -= nu*log(denom);
		}

	}
	return Q;
}

double EM::addRegularizers( std::vector<double> &grads, Param *tmp_param ){

	double Q = 0.0;
	std::set<unsigned int>::iterator it = ((MasterComms *)comm)->master_used_idxs.begin();
	for( ; it != ((MasterComms *)comm)->master_used_idxs.end(); ++it ){
		
		double weight = tmp_param->getWeightAtIdx(*it);
		Q -= 0.5*cfg->lambda*weight*weight;
		grads[*it] -= cfg->lambda*weight;		
	}

	//Remove the Bias terms (don't regularize the bias terms!)
	unsigned int weights_per_energy = tmp_param->getNumWeightsPerEnergyLevel();
	for( unsigned int energy = 0; energy < tmp_param->getNumEnergyLevels(); energy++ ){
		double bias = tmp_param->getWeightAtIdx(energy * weights_per_energy);
		Q += 0.5*cfg->lambda*bias*bias; 
		grads[energy * weights_per_energy] += cfg->lambda*bias;
	}
	return Q;
}

void EM::zeroUnusedParams(){

	unsigned int i;
	for( i = 0; i < param->getNumWeights(); i++ ){
		if( ((MasterComms *)comm)->master_used_idxs.find(i) == ((MasterComms *)comm)->master_used_idxs.end() )
			param->setWeightAtIdx(0.0, i);
	}
	
}

int EM::updateLineSearchConverged(double Qsum, double tmpQsum, double prev_tmpQsum, double tmp_step, double grad_sq){
	
	//Usual convergence rule
	if( tmpQsum >= Qsum + cfg->line_search_alpha*tmp_step*grad_sq ) return 1;
	
	//Rule that seems necessary due numerical imprecision...?
	if( tmpQsum < prev_tmpQsum ) return 1;
	
	return 0;
}

double EM::computeGradientSquaredTerm( std::vector<double> &grads){

	double grad_sq = 0.0;
	std::set<unsigned int>::iterator it = ((MasterComms *)comm)->master_used_idxs.begin();
	for( ; it != ((MasterComms *)comm)->master_used_idxs.end(); ++it ) 
		grad_sq += (grads[*it]*grads[*it]);
	return grad_sq;
}

void EM::copyTmpToParamWeights( Param &tmp_params ){

	std::set<unsigned int>::iterator it = ((MasterComms *)comm)->master_used_idxs.begin();
	for( ; it != ((MasterComms *)comm)->master_used_idxs.end(); ++it )
		param->setWeightAtIdx( tmp_params.getWeightAtIdx(*it), *it );
}

void EM::copyGrads( std::vector<double> &tmp_grads, std::vector<double> &grads ){
	
	std::set<unsigned int>::iterator it = ((MasterComms *)comm)->master_used_idxs.begin();
	for( ; it != ((MasterComms *)comm)->master_used_idxs.end(); ++it )
		grads[*it] = tmp_grads[*it];
}

void EM::resetTmpGradients( std::vector<double> &tmp_grads ){
	
	if( comm->isMaster() ){
		std::set<unsigned int>::iterator it = ((MasterComms *)comm)->master_used_idxs.begin();
		for( ; it != ((MasterComms *)comm)->master_used_idxs.end(); ++it ) 
			tmp_grads[*it] = 0.0;
	}
	else{
		std::set<unsigned int>::iterator it = comm->used_idxs.begin();
		for( ; it != comm->used_idxs.end(); ++it ) tmp_grads[*it] = 0.0;
	}
}
