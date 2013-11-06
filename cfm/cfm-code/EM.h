/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# EM.h
#
# Description: 	Class to apply Expectation Maximization algorithm to derive 
#				model parameters.
#					E-step: IPFP or equivalent.
#					M-step: Gradient Ascent
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __EM_TRAIN_H__
#define __EM_TRAIN_H__

static const int MAX_EM_ITERATIONS = 100;

#include "Config.h"
#include "MolData.h"
#include "Param.h"
#include "Comms.h"
#include "IPFP.h"

//Sufficient stats, access by transition index for a given molecule
typedef std::vector<double> suft_t;

struct suft_counts_t{
	//Access each by molecule index
	std::vector<suft_t> values;
};

class EM{
public:
	//Constructor
	//Note: To include the group in the status filename, include _GRP_ in the name
	EM(config_t *a_cfg, FeatureCalculator *an_fc, std::string &a_status_filename);
	~EM();

	//Run the EM algorithm on the supplied data (except the specified group), 
	//return the final likelihood value.
	double run( std::vector<MolData> &data, int group, bool fullZeroInit = false );
	
	//After running EM, the final params can be written out to file
	void writeParamsToFile( std::string &filename );

	//This is public so the test can access it....there must be a better way?
	double computeAndAccumulateGradient(std::vector<double> &grads, int molidx, MolData &moldata, Param *tmp_params, suft_counts_t &suft, std::set<unsigned int> &used_idxs);

private:
	//The feature calculator to use - preconfigured with feature spec
	FeatureCalculator *fc;

	//The current parameters
	Param *param;

	//Configuration data
	config_t *cfg;

	//Communicator (for exchanging data between threads)
	Comms *comm;
	int unused_zeroed;	//Use to note when unused parameters have been zeroed (so we don't
						//do it more times than we need to)

	//For writing status messages to a log file
	std::string status_filename;
	void writeStatus( const char *msg );

	//Initialise sufficient statistics
	void initSuft( suft_counts_t &suft, std::vector<MolData> &data );

	//Update sufficient statistics based on beliefs
	void recordSufficientStatistics( suft_counts_t &suft, int molidx, MolData *moldata, beliefs_t *beliefs );

	//M-step: Run gradient 
	double updateParameters( std::vector<MolData> &data, suft_counts_t &suft );

	//Helper functions - used within UpdateParameters
	void stepAndSetTemporaryParams( Param &tmp_params, std::vector<double> &grads_sum, double tmp_step);
	double addRegularizers( std::vector<double> &grads, Param *tmp_param );
	void zeroUnusedParams();
	int updateLineSearchConverged(double Qsum, double tmpQsum, double prev_tmpQsum, double tmp_step, double grad_sq);
	double computeGradientSquaredTerm( std::vector<double> &grads);
	void copyTmpToParamWeights( Param &tmp_params );
	void copyGrads( std::vector<double> &tmp_grads, std::vector<double> &grads );
	void resetTmpGradients( std::vector<double> &tmp_grads );
};

#endif // __EM_TRAIN_H__
