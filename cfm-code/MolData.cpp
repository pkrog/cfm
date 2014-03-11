/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# MolData.cpp
#
# Description: 	Class to hold the input data belonging to a molecule:
#					 - An ID and smiles/inchi
#					 - (optional) A computed fragmentation graph
#					 - (optional) A computed set of features corresponding to that graph
#					 - (optional) A computed set of theta values for that graph
#				     - (optional) A computed set of transition probabilities using those thetas.
#					 - (optional) A set of spectra
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "MolData.h"
#include "FragmentGraphGenerator.h"
#include "Inference.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>

#include <string>
#include <iostream>
#include <fstream>

void MolData::computeFragmentGraph( int depth ){
	
	FragmentGraphGenerator fgen;
	fg = fgen.createNewGraph();
	FragmentTreeNode *startnode = fgen.createStartNode( smiles_or_inchi );
	fgen.compute( *startnode, depth, -1 );
	delete startnode;
	graph_computed = 1;

}

void MolData::computeFragmentGraphAndReplaceMolsWithFVs( int depth, FeatureCalculator *fc, bool retain_smiles){
	
	//Compute the fragment graph, replacing the transition molecules with feature vectors
	FragmentGraphGenerator fgen( fc );
	fg = fgen.createNewGraph();
	FragmentTreeNode *startnode = fgen.createStartNode( smiles_or_inchi );
	fgen.compute( *startnode, depth, -1 );
	delete startnode;
	graph_computed = 1;

	//Copy all the feature vector pointers up into the mol data
	fvs.resize( fg->getNumTransitions() );
	for( unsigned int i = 0; i < fg->getNumTransitions(); i++ )
		fvs[i] = fg->getTransitionAtIdx(i)->getTmpFV();

	//Delete all the fragment smiles (we only need these while we're computing the graph)
	if(!retain_smiles) fg->clearAllSmiles();

}

void MolData::computeLikelyFragmentGraphAndSetThetas( Param *param, config_t *cfg, double prob_thresh_for_prune ){
	
	//Compute the fragment graph, replacing the transition molecules with feature vectors
	LikelyFragmentGraphGenerator fgen(param, cfg, prob_thresh_for_prune ); 
	fg = fgen.createNewGraph();
	FragmentTreeNode *startnode = fgen.createStartNode( smiles_or_inchi );
	fgen.compute( *startnode, cfg->fg_depth, -1, 0.0 );
	delete startnode;
	graph_computed = 1;

	//Copy all the theta values up into the mol data and delete the tmp thetas
	const unsigned int num_levels = param->getNumEnergyLevels();
	thetas.resize(num_levels);
	for( unsigned int energy = 0; energy < num_levels; energy++ ){
		thetas[energy].resize( fg->getNumTransitions() );
		for( unsigned int i = 0; i < fg->getNumTransitions(); i++ )
			thetas[energy][i] = (*(fg->getTransitionAtIdx(i)->getTmpThetas()))[energy];
	}
	
	//Delete all the fragment smiles (we only need these while we're computing the graph)
	fg->clearAllSmiles();

}

//Function to compute a much reduced fragment graph containing only those
//fragmentations as actually occur in the spectra, based on a computed set of beliefs
//thresholding inclusion in the graph by the provided belief_thresh value (log domain)
void MolData::computeEvidenceFragmentGraph( beliefs_t *beliefs, double log_belief_thresh, config_t *cfg ){

	unsigned int num_transitions = fg->getNumTransitions();
	unsigned int num_fragments = fg->getNumFragments();
	ev_fg = new EvidenceFragmentGraph();
	
	//Add the root fragment
	Transition null_t;
	std::vector<int> id_lookup(num_fragments,-1);
	std::vector<double> main_ev;
	computeFragmentEvidenceValues(main_ev, 0, beliefs, cfg );
	id_lookup[0] = ev_fg->addToGraphDirectNoCheck(EvidenceFragment(*(fg->getFragmentAtIdx(0)),-1,main_ev),&null_t,-1);



	//Add the fragments and transitions if the beliefs are high enough
	std::vector<int> t_added_flags(num_transitions, 0);
	for( int depth = 0; depth < beliefs->tn[0].size(); depth++ ){

		for( unsigned int i = 0; i < num_transitions; i++ ){
			if(t_added_flags[i]) continue;

			if( beliefs->tn[i][depth] >= log_belief_thresh ){
				const Transition *t = fg->getTransitionAtIdx(i);
				const Fragment *f = fg->getFragmentAtIdx( t->getToId() );
				int from_id = id_lookup[t->getFromId()];
				if( id_lookup[t->getToId()] == -1 ){
					std::vector<double> evidence;
					computeFragmentEvidenceValues(evidence, t->getToId(), beliefs, cfg );
					int to_id = ev_fg->addToGraphDirectNoCheck(EvidenceFragment(*f,-1,evidence),t,from_id);
					id_lookup[t->getToId()] = to_id;
				}
				else if( id_lookup[t->getFromId()] != -1)
					ev_fg->addTransition(from_id, id_lookup[t->getToId()], t->getNLSmiles());
				t_added_flags[i] = 1;
			}
		}
	}
	ev_graph_computed = 1;
}

void MolData::computeFragmentEvidenceValues(std::vector<double> &evidence, int frag_idx, const beliefs_t *beliefs, config_t *cfg ){

	const tmap_t *tomap = fg->getToIdTMap();

	evidence.resize(cfg->dv_spectrum_depths.size());
	for( int energy = 0; energy < cfg->dv_spectrum_depths.size(); energy++ ){
	
		//Find out the depth of the spectrum of interest in the beliefs
		int depth = cfg->dv_spectrum_depths[energy] - 1;
		if(cfg->use_single_energy_cfm )
			for( int i = 0; i < energy;  i++ ) 
				depth += cfg->dv_spectrum_depths[i];
		
		//Compute the accumulated belief for the fragment of interest at the depth of interest
		evidence[energy] = beliefs->ps[frag_idx][depth];
		std::vector<int>::const_iterator it = (*tomap)[frag_idx].begin();
		for( ; it != (*tomap)[frag_idx].end(); ++it )
			evidence[energy] = logAdd( evidence[energy], beliefs->tn[*it][depth]);
	}

}

void MolData::annotatePeaks( double abs_tol, double ppm_tol, bool prune_deadends){

	if(!ev_graph_computed){
		std::cout  << "Warning: called annotate peaks without first computing evidence graph" << std::endl;
		return;
	}

	unsigned int num_fragments = ev_fg->getNumFragments();
	std::vector<int> frag_flags(num_fragments,0);

	//Assign fragment annotations to peaks
	for( int energy = 0; energy < spectra.size(); energy++ ){
		Spectrum::iterator it = spectra[energy].begin();
		for( ; it != spectra[energy].end();  ++it ){
		
			double mass_tol = getMassTol( abs_tol, ppm_tol, it->mass );
			int num_fragments = ev_fg->getNumFragments();
			for( int fidx = 0; fidx < num_fragments; fidx++ ){
			
				const EvidenceFragment *f = ev_fg->getFragmentAtIdx(fidx);
				if( fabs( f->getMass() - it->mass ) <= mass_tol ){
					it->annotations.push_back(annotation_t(f->getId(),f->getEvidence(energy)));	
					frag_flags[fidx] = 1;
				}
			}
		}
	}

	//Sort the annotations by evidence score
	for( int energy = 0; energy < spectra.size(); energy++ ){
		Spectrum::iterator it = spectra[energy].begin();
		for( ; it != spectra[energy].end();  ++it ){
			std::sort( it->annotations.begin(), it->annotations.end(), sort_annotations_by_score );
		}
	}

	if(!prune_deadends) return;

	std::vector<int> delete_flags(num_fragments);
	for(unsigned int i = 0; i < num_fragments; i++ ){	
		//If there is no peak for this fragment, check there is
		//an annotated descendent, else flag for deletion
		std::vector<int> direct_flags(num_fragments);
		ev_fg->setFlagsForDirectPaths( direct_flags, 0, frag_flags );
		if( ev_fg->fragmentIsRedundant(i, frag_flags, direct_flags) )
			delete_flags[i] = 1;
		else delete_flags[i] = 0;
	}

	//Rebuild the graph again without deleted fragments
	Transition null_t;
	EvidenceFragmentGraph *new_ev_fg = new EvidenceFragmentGraph();
	std::vector<int> id_lookup(num_fragments, -1);
	id_lookup[0] = new_ev_fg->addToGraphDirectNoCheck( *ev_fg->getFragmentAtIdx(0), &null_t, -1);
	for(unsigned int i = 0; i < ev_fg->getNumTransitions(); i++ ){	
		const Transition *t = ev_fg->getTransitionAtIdx(i);
		int fromid = t->getFromId();
		int toid = t->getToId();
		if( delete_flags[toid] || delete_flags[fromid] ) continue;
		if( id_lookup[toid] == -1 ) 
			id_lookup[toid] = new_ev_fg->addToGraphDirectNoCheck(*ev_fg->getFragmentAtIdx(toid), t, id_lookup[fromid]);
		else
			new_ev_fg->addTransition(id_lookup[fromid], id_lookup[toid], t->getNLSmiles() );
	}
	delete ev_fg;
	ev_fg = new_ev_fg;

	//Update annotations with new ids
	for( int energy = 0; energy < spectra.size(); energy++ ){
		Spectrum::iterator it = spectra[energy].begin();
		for( ; it != spectra[energy].end();  ++it ){
			for( unsigned int i = 0; i < it->annotations.size(); i++)
				it->annotations[i].first = id_lookup[it->annotations[i].first];
		}
	}

}



void MolData::computeFeatureVectors( FeatureCalculator *fc, bool deleteMols ){

	fvs.resize( fg->getNumTransitions() );
	for( unsigned int i = 0; i < fg->getNumTransitions(); i++ ){
		const Transition *t = fg->getTransitionAtIdx(i);
		fvs[i] = fc->computeFV( t->getIon(), t->getNeutralLoss() );
		if( deleteMols ) fg->deleteMolsForTransitionAtIdx(i);
	}
	if( !deleteMols ) return;
	fg->clearAllSmiles();

}

void MolData::computeTransitionThetas( Param &param ){

	const unsigned int num_levels = param.getNumEnergyLevels();
	thetas.resize(num_levels);
	for( unsigned int energy = 0; energy < num_levels; energy++ ){
		
		//Compute the theta value for each feature vector
		thetas[energy].resize( fg->getNumTransitions() );
		for( unsigned int i = 0; i < fg->getNumTransitions(); i++ )
			thetas[energy][i] = param.computeTheta( *fvs[i], energy );		
	}
}

void MolData::computeTransitionProbabilities(){

	log_probs.resize( thetas.size() );
	for( unsigned int energy = 0; energy < thetas.size(); energy++ ){
		
		//Will store an entry for all transitions, followed by all persistences
		log_probs[energy].resize( fg->getNumTransitions() + fg->getNumFragments() );

		//Compute all the denominators
		std::vector<double> denom_cache( fg->getNumFragments() );
		const tmap_t *from_id_map = fg->getFromIdTMap();
		for( unsigned int i = 0; i < fg->getNumFragments(); i++ ){
		
			double denom = 0.0;
			std::vector<int>::const_iterator it = (*from_id_map)[i].begin();
			for( ; it != (*from_id_map)[i].end(); ++it )
				denom = logAdd( denom, thetas[energy][*it] );
			denom_cache[i] = denom;
		}
	
		//Set the transition log probabilities
		for( unsigned int i = 0; i < fg->getNumTransitions(); i++ ){
			const Transition *t = fg->getTransitionAtIdx(i);
			log_probs[energy][i] = thetas[energy][i] - denom_cache[ t->getFromId() ];
		}

		//Set the persistence log probabilities
		int offset = fg->getNumTransitions();
		for( unsigned int i = 0; i < fg->getNumFragments(); i++ )
			log_probs[energy][offset + i] = -denom_cache[i];

	}
	return;
}

void MolData::readInSpectraFromFile( const std::string &peak_filename, bool readToPredicted ){

	//Set the spectra that we are writing to
	std::vector<Spectrum> *spec_dest = &spectra;
	if( readToPredicted ) spec_dest = &predicted_spectra;

	//Clear any existing entries
	spec_dest->clear();
	spec_dest->resize(0);

	Spectrum *curr_spec;
	std::string line;
	std::ifstream ifs ( peak_filename.c_str() , std::ifstream::in );
	if(!ifs) std::cout << "Warning: Could not open file " << peak_filename << std::endl;
	int first = 1;
	while( ifs.good() ){
		getline( ifs, line );
		if( line == "" ) continue;
		
		//Check for the energy specifier - start a new spectrum if found 
		//or start one anyway if there is no energy specifier
		if( line.substr(0,3) == "low" 
			|| line.substr(0,3) == "med"
			|| line.substr(0,4) == "high" 
			|| line.substr(0,6) == "energy" ){
				spec_dest->push_back( Spectrum() );
				curr_spec = &(spec_dest->back());
				first = 0;
				continue;
		}
		else if( first ){
			spec_dest->push_back( Spectrum() );
			curr_spec = &(spec_dest->back());	
		}
		first = 0;

		//Otherwise allocate peak to the current spectrum
		std::stringstream ss(line);
		double mass, intensity;
		ss >> mass >> intensity;
		curr_spec->push_back( Peak(mass, intensity) );
	}
	std::vector<Spectrum>::iterator it = spec_dest->begin();
	for( ; it != spec_dest->end(); ++it ) normalizeAndSortSpectrum(*it);
	ifs.close();
}

void MolData::removePeaksWithNoFragment( double abs_tol, double ppm_tol ){

	std::vector<Spectrum>::iterator it = spectra.begin();
	for( ; it != spectra.end(); ++it ){ 
		
		//Remove any peaks more than mass_tol away from any fragment
		Spectrum::iterator itp = it->begin();
		for( ; itp != it->end(); ){
			
			double mass_tol = getMassTol( abs_tol, ppm_tol, itp->mass );
			bool found = false;
			for( unsigned int i = 0; i < fg->getNumFragments(); i++ ){
				const Fragment *f = fg->getFragmentAtIdx(i);
				if( fabs( f->getMass() - itp->mass ) < mass_tol ){ 
					found = true;
					break;
				}
			}
			if(!found) itp = it->erase( itp );
			else ++itp;
		}

		//Renormalise
		normalizeAndSortSpectrum(*it);
	}
}

void MolData::createInterpolatedSpectra( config_t &cfg ){
	
	//Actual spectra
	spectra.resize(cfg.model_depth);
	for( int i = cfg.spectrum_depths.size()-1; i >= 0 ; i-- )
		spectra[ cfg.spectrum_depths[i] - 1 ] = spectra[i];
		
	//Create a preliminary spectrum with just the main peak
	Spectrum initial_spec;
	initial_spec.push_back( Peak(fg->getFragmentAtIdx(0)->getMass(), 100.0)  );

	//Interpolated spectra
	int prev_d = 0;
	std::vector<int>::iterator it = cfg.spectrum_depths.begin();
	for( ; it != cfg.spectrum_depths.end(); ++it ){
		int d = *it; 
		for( int i = prev_d; i < d-1; i++ ){
			double ratio = (double)(i+1-prev_d)/(d-prev_d);
			Spectrum *low_spec;
			if( prev_d > 0 ) low_spec = &(spectra[prev_d-1]);
			else low_spec = &initial_spec;
			interpolateSpectra( spectra[i], *low_spec, spectra[d-1], ratio, cfg);
		}
		prev_d = d; 
	}
}

void MolData::interpolateSpectra( Spectrum &output, Spectrum &lower_spec, Spectrum &higher_spec, double ratio, config_t &cfg  ){
	
	//Ensure the spectra are sorted
	std::sort( lower_spec.begin(), lower_spec.end(), sort_peaks_by_mass );
	std::sort( higher_spec.begin(), higher_spec.end(), sort_peaks_by_mass );

	//Now do the interpolation, combining peaks within mass_tol in the ratio given
	output.clear();
	Spectrum::iterator it_low = lower_spec.begin();
	Spectrum::iterator it_high = higher_spec.begin();
	while( it_low != lower_spec.end() && it_high != higher_spec.end() ){
		double mass = 0.5*(it_low->mass + it_high->mass);
		double mass_tol = getMassTol( cfg.abs_mass_tol, cfg.ppm_mass_tol, mass );
		if( fabs( it_low->mass - it_high->mass ) < mass_tol ){
			double intensity = (1-ratio)*it_low->intensity + ratio*it_high->intensity;
			output.push_back( Peak(mass, intensity) );
			++it_low;
			++it_high;
		} 
		else if( it_low->mass > it_high->mass ){
			output.push_back( *it_high );
			output.back().intensity *= ratio;
			++it_high;
		}
		else{
			output.push_back( *it_low );
			output.back().intensity *= (1-ratio);
			++it_low;		
		}
	}
	while( it_low != lower_spec.end() ){
		output.push_back( *it_low );
		output.back().intensity *= (1-ratio);
		++it_low;		
	}
	while( it_high != higher_spec.end() ){
		output.push_back( *it_high );
		output.back().intensity *= ratio;
		++it_high;		
	}

}

	
void MolData::computePredictedSpectra( Param &param, config_t &cfg, bool postprocess, bool use_existing_thetas ){

	//Divert to other function if doing single energy CFM
	if( cfg.use_single_energy_cfm ){ 
		computePredictedSingleEnergySpectra( param, cfg, postprocess, use_existing_thetas );
		return;
	}

	//Compute the transition probabilities using this parameter set
	if(!use_existing_thetas) computeTransitionThetas( param );
	computeTransitionProbabilities();

	//Run forward inference
	std::vector<Message> msgs;
	runInferenceDownwardPass( msgs, this, cfg, cfg.model_depth );

	//Generate and collect the peak results
	std::string outfilename;
	predicted_spectra.resize( cfg.spectrum_depths.size() );
	for( unsigned int energy = 0; energy < cfg.spectrum_depths.size(); energy++ ){
		predicted_spectra[energy].clear();

		//Extract the peaks from the relevant message
		int msg_depth = cfg.spectrum_depths[energy]-1;		
		Message *msg = &(msgs[msg_depth]);
		translatePeaksFromMsgToSpectra( predicted_spectra[energy], msg);

	}

	if( postprocess ) postprocessPredictedSpectra();
}

void MolData::computePredictedSingleEnergySpectra( Param &param, config_t &cfg, bool postprocess, bool use_existing_thetas ){

	//Compute the transition probabilities using this parameter set
	if( !use_existing_thetas ) computeTransitionThetas( param );
	computeTransitionProbabilities();

	//Generate and collect the peak results
	std::string outfilename;
	predicted_spectra.resize( cfg.spectrum_depths.size() );
	for( unsigned int energy = 0; energy < cfg.spectrum_depths.size(); energy++ ){
		predicted_spectra[energy].clear();

		config_t se_cfg;
		initSingleEnergyConfig( se_cfg, cfg, energy );

		//Run forward inference
		std::vector<Message> msgs;
		runInferenceDownwardPass( msgs, this, se_cfg, se_cfg.model_depth );

		//Extract the peaks from the relevant message
		int msg_depth = se_cfg.spectrum_depths[0]-1;		
		Message *msg = &(msgs[msg_depth]);
		translatePeaksFromMsgToSpectra( predicted_spectra[energy], msg );
	}

	if( postprocess ) postprocessPredictedSpectra();
}

void MolData::translatePeaksFromMsgToSpectra( Spectrum &out_spec, Message *msg ){

	//Create the peaks
	std::map<double, double> peak_probs;
	Message::const_iterator itt = msg->begin();
	for( ; itt != msg->end() ; ++itt ){
		double mass = fg->getFragmentAtIdx(itt.index())->getMass();
		if( peak_probs.find(mass) != peak_probs.end() )
			peak_probs[mass] += exp(*itt);
		else
			peak_probs[mass] = exp(*itt);
	}		
		
	//Add the peaks to the spectra
	std::map<double, double>::iterator itm = peak_probs.begin();
	for( ; itm != peak_probs.end(); ++itm )
		out_spec.push_back( Peak(itm->first, itm->second*100.0) );
}

void MolData::writePredictedSpectraToFile( std::string &filename ){

	std::ofstream of;
	of.open(filename.c_str());
	if( !of.is_open() ){
		std::cout << "Warning: Trouble opening predicted spectrum file" << std::endl;
		return;
	}
	std::streambuf *buf = of.rdbuf();	
	std::ostream out(buf);
	outputSpectra( out, "Predicted" );
	of.close();
}

void MolData::writeFullEnumerationSpectrumToFile( std::string &filename ){
	
	std::ofstream of;
	of.open(filename.c_str());
	if( !of.is_open() ){
		std::cout << "Warning: Trouble opening enumerated spectrum file" << std::endl;
		return;
	}
	std::streambuf *buf = of.rdbuf();	
	std::ostream out(buf);

	//Get and sort the fragment masses
	unsigned int numf = fg->getNumFragments();
	std::vector<double> all_masses(numf);
	for( unsigned int i = 0; i < numf; i++ ){
		const Fragment *f = fg->getFragmentAtIdx(i);
		all_masses[i] = f->getMass();
	}
	std::sort( all_masses.begin(), all_masses.end() );

	//Remove repeats
	std::vector<double> all_unique_masses(numf);
	int i_unique = 0;
	double prev_mass = -1.0, tol = 1e-10;
	for( int i = 0; i < numf; i++ ){
		if( fabs(all_masses[i] - prev_mass) > tol )
			all_unique_masses[i_unique++] = all_masses[i];
		prev_mass = all_masses[i];
	}
	int num_unique = i_unique;

	//All peaks of uniform intensity at each possible fragment mass
	for( unsigned int energy = 0; energy < spectra.size(); energy++ ){
		out << "energy" << energy << std::endl;

		out << std::setprecision(10);
		double peak_height = 100.0/(double)num_unique;
		for( int i = 0; i < num_unique; i++ )
			out << all_unique_masses[i] << " " << peak_height << std::endl;
	}
	of.close();
}


void MolData::outputSpectra( std::ostream &out, const char*spec_type ){
	
	std::vector<Spectrum> *spectra_to_output;
	if(std::string(spec_type) == "Predicted") spectra_to_output = &predicted_spectra;
	else if(std::string(spec_type) == "Annotated") spectra_to_output = &spectra;
	else std::cout << "Unknown spectrum type to output: " << spec_type << std::endl;

	std::vector<Spectrum>::iterator it = spectra_to_output->begin();
	for( int energy = 0; it != spectra_to_output->end(); ++it, energy++ ){
		out << "energy" << energy << std::endl;

		Spectrum::iterator itp = it->begin();
		out << std::setprecision(10);
		for( ; itp != it->end(); ++itp ){
			out << itp->mass << " " << itp->intensity;
			std::vector<annotation_t>::iterator ita = itp->annotations.begin();
			for( ; ita != itp->annotations.end(); ++ita )
				out << " " << ita->first; // << " (" << ita->second << ") ";
			out << std::endl;
		}
	}
}

void MolData::postprocessPredictedSpectra(){
	
	std::vector<Spectrum>::iterator it = predicted_spectra.begin();
	for( ; it != predicted_spectra.end(); ++it){
		postprocessSpectrum( *it );
		normalizeAndSortSpectrum( *it );
	}
}

void MolData::postprocessSpectrum( Spectrum &spectrum ){

	std::sort( spectrum.begin(), spectrum.end(), sort_peaks_by_intensity );
	double total = 0.0;
	Spectrum::iterator it = spectrum.begin();
	int count = 0;
	for( ; it != spectrum.end(); ++it ){
		total += it->intensity;
		count++;
		//Take the top 80% of energy (assuming at least 5 peaks), 
		//or the highest 30 peaks (whichever comes first)
		if( (total > 80.0 && count > 5 ) || count > 30 ) break;	
	}
	spectrum.resize( count );
	std::sort( spectrum.begin(), spectrum.end(), sort_peaks_by_mass );
}

void MolData::normalizeAndSortSpectrum( Spectrum &spectrum){
		
	//Compute the normalizer
	double sum = 0.0;
	Spectrum::iterator itp = spectrum.begin();
	for( ; itp != spectrum.end(); ++itp )
		sum += itp->intensity;
	double norm = 100.0/sum;

	//Adjust the values
	for( itp = spectrum.begin(); itp != spectrum.end(); ++itp )
		itp->intensity *= norm;

	//Ensure the peaks are sorted by mass
	std::sort( spectrum.begin(), spectrum.end(), sort_peaks_by_mass );
}

MolData::~MolData(){

	if( graph_computed ) delete fg;
	if( ev_graph_computed ) delete ev_fg;

	//Delete any computed feature vectors
	std::vector<FeatureVector *>::iterator it = fvs.begin();
	for( ; it != fvs.end(); ++it ) delete *it;

}
