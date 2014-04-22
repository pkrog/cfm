/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Config.cpp
#
# Description: 	Structs, functions, defaults for setting general 
#			    configuration data.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Config.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

void initDefaultConfig( config_t &cfg ){

	cfg.lambda = DEFAULT_LAMBDA;
	cfg.converge_count_thresh = DEFAULT_CONVERGE_COUNT_THRESH;
	cfg.em_converge_thresh = DEFAULT_EM_CONVERGE_THRESH;
	cfg.ga_converge_thresh = DEFAULT_GA_CONVERGE_THRESH;
	cfg.model_depth = DEFAULT_MODEL_DEPTH;
	cfg.abs_mass_tol = DEFAULT_ABS_MASS_TOL;
	cfg.ppm_mass_tol = DEFAULT_PPM_MASS_TOL;
	cfg.num_em_restarts = DEFAULT_NUM_EM_RESTARTS;
	cfg.lambda_hold = DEFAULT_LAMBDA_HOLD;
	cfg.line_search_alpha = DEFAULT_LINE_SEARCH_ALPHA;
	cfg.line_search_beta = DEFAULT_LINE_SEARCH_BETA;
	cfg.starting_step_size = 1.0;
	cfg.max_search_count = DEFAULT_MAX_SEARCH_COUNT;
	cfg.spectrum_depths.clear();
	cfg.spectrum_weights.clear();
	cfg.dv_spectrum_indexes.clear();
	cfg.interpolate_spectra = 0;
	cfg.intermediate_weights = 0.0;
	cfg.fg_depth = DEFAULT_FRAGGRAPH_DEPTH;
	cfg.ipfp_algorithm = DEFAULT_IPFP_ALGORITHM;
	cfg.ipfp_converge_thresh = DEFAULT_IPFP_CONVERGE_THRESH;
	cfg.osc_ipfp_converge_thresh = DEFAULT_IPFP_OSC_CONVERGE_THRESH;
	cfg.use_single_energy_cfm = 0;
	cfg.ionization_mode = DEFAULT_IONIZATION_MODE;
}


void initConfig( config_t &cfg, std::string &filename ){

	std::string line, name;
	double value;
	std::ifstream ifs ( filename.c_str(), std::ifstream::in  );

	initDefaultConfig( cfg );

	//Read the config file into the paramater update config structure
	if(!ifs) std::cout << "Could not open file " << filename << std::endl;
	while( ifs.good() ){

		getline( ifs, line );
		if( line.size() < 3 ) continue;	//in case of empty line

		std::stringstream ss1(line);
		ss1 >> name >> value;

		if( name == "lambda" ) cfg.lambda = value;
		else if( name == "ionization_mode" ) cfg.ionization_mode = (int)value;
		else if( name == "converge_count_thresh" ) cfg.converge_count_thresh = (int)value;
		else if( name == "em_converge_thresh" ) cfg.em_converge_thresh = value;
		else if( name == "ga_converge_thresh" ) cfg.ga_converge_thresh = value;
		else if( name == "model_depth" ) cfg.model_depth = (unsigned int)value;
		else if( name == "spectrum_depth" ) cfg.spectrum_depths.push_back( (unsigned int)value );
		else if( name == "spectrum_weight" ) cfg.spectrum_weights.push_back( (double)value );
		else if( name == "interpolate_spectra" ) cfg.interpolate_spectra = (int)value;
		else if( name == "intermediate_weights" ) cfg.intermediate_weights = (double)value;
		else if( name == "abs_mass_tol" ) cfg.abs_mass_tol = (double)value;
		else if( name == "ppm_mass_tol" ) cfg.ppm_mass_tol = (double)value;
		else if( name == "num_em_restarts" ) cfg.num_em_restarts = (int)value;
		else if( name == "lambda_hold" ) cfg.lambda_hold = (double)value;
		else if( name == "line_search_alpha" ) cfg.line_search_alpha = (double)value;
		else if( name == "line_search_beta" ) cfg.line_search_beta = (double)value;
		else if( name == "starting_step_size" ) cfg.starting_step_size = (double)value;
		else if( name == "max_search_count" ) cfg.max_search_count = (int)value;
		else if( name == "fg_depth" ) cfg.fg_depth = (int)value;
		else if( name == "ipfp_algorithm" ) cfg.ipfp_algorithm = (int)value;
		else if( name == "ipfp_converge_thresh" ) cfg.ipfp_converge_thresh = (double)value;
		else if( name == "osc_ipfp_converge_thresh" ) cfg.osc_ipfp_converge_thresh = (double)value;
		else if( name == "use_single_energy_cfm" ) cfg.use_single_energy_cfm = (int)value;
		else std::cout << "Warning: Unknown paramater configuration identifier " << name << std::endl;
	}
	ifs.close();

	if( cfg.spectrum_depths.size() != cfg.spectrum_weights.size() )
		std::cout << "Warning: Mismatch between size of spectrum depths and weights" << std::endl;
	
	initDerivedConfig(cfg);

	//Report config parameters
	if( cfg.use_single_energy_cfm ) std::cout << "Using Single Energy CFM" << std::endl;
	else std::cout << "Using Combined Energy CFM" << std::endl;
	if( cfg.ionization_mode == POSITIVE_IONIZATION_MODE ) std::cout << "Positive Ionization Mode" << std::endl;
	else if( cfg.ionization_mode == NEGATIVE_IONIZATION_MODE ) std::cout << "Negative Ionization Mode" << std::endl;
	else{ 
		std::cout << "Warning: Unknown Ionization Mode, reverting to default mode (positive)" << std::endl;
		cfg.ionization_mode = DEFAULT_IONIZATION_MODE;
	}
	std::cout << "Using Lambda " << cfg.lambda << std::endl;
	std::cout << "Using Lambda Hold " << cfg.lambda_hold << std::endl;
	std::cout << "Using Converge Count Threshold " << cfg.converge_count_thresh << std::endl;
	std::cout << "Using EM Convergence Threshold " << cfg.em_converge_thresh << std::endl;
	std::cout << "Using GA Convergence Threshold " << cfg.ga_converge_thresh << std::endl;
	std::cout << "Using Fragmentation Graph Depth " << cfg.fg_depth << std::endl;
	std::cout << "Using Model Depth " << cfg.model_depth << std::endl;
	std::cout << "Using Spectrum Depths and Weights: ";
	for( unsigned int i =0; i < cfg.spectrum_depths.size(); i++ )
		std::cout << "(" << cfg.spectrum_depths[i] << "," << cfg.spectrum_weights[i] << ") ";
	std::cout << std::endl;
	if(cfg.interpolate_spectra) std::cout << "Using interpolated spectra with intermediate weights=" << cfg.intermediate_weights << std::endl;
	else std::cout << "Not interpolated spectra" << std::endl;
	std::cout << "Using Absolute mass tolerance " << cfg.abs_mass_tol << std::endl;
	std::cout << "Using PPM mass tolerance" << cfg.ppm_mass_tol << std::endl;
	std::cout << "Using Line Search Alpha Beta " << cfg.line_search_alpha << " " << cfg.line_search_beta << std::endl;
	std::cout << "Using Starting Step Size " << cfg.starting_step_size << std::endl;
	std::cout << "Using Max Line Search Count " << cfg.max_search_count << std::endl;
	if( cfg.ipfp_algorithm == 0 ) std::cout << "Using standard IPFP" << std::endl;
	else if( cfg.ipfp_algorithm == 1 ) std::cout << "Using GEMA" << std::endl;
	else if( cfg.ipfp_algorithm == 2 ) std::cout << "Using IPFP with Oscillatory Adjustment" << std::endl;
	else std::cout << "Warning: Unknown IPFP algorithm id" << std::endl;
	std::cout << "Using IPFP Converge Thresh " << cfg.ipfp_converge_thresh << std::endl;
	std::cout << "Using IPFP Oscillatory Converge Thresh " << cfg.osc_ipfp_converge_thresh << std::endl;

}

void initDerivedConfig( config_t &cfg, int se_energy ){

	//Derived Parameters
	cfg.map_d_to_energy.resize( cfg.model_depth );
	int energy = 0;
	if( se_energy > 0 ) energy = se_energy;
	for( unsigned int d = 0; d < cfg.model_depth; d++ ){
		std::vector<int>::iterator it = cfg.spectrum_depths.begin();
		for( ; it != cfg.spectrum_depths.end(); ++it ){
			if( d == *it ) energy++;
		}
		cfg.map_d_to_energy[d] = energy;
	}

	if(cfg.interpolate_spectra){
		cfg.dv_spectrum_depths.resize( cfg.model_depth );
		cfg.dv_spectrum_weights.resize( cfg.model_depth );

		for( unsigned int d = 0; d < cfg.model_depth; d++ ){
			cfg.dv_spectrum_depths[d] = d+1;
			int energy = cfg.map_d_to_energy[d];
			if( cfg.spectrum_depths[energy] == d+1 )
				cfg.dv_spectrum_weights[d] = cfg.spectrum_weights[energy];
			else
				cfg.dv_spectrum_weights[d] = cfg.intermediate_weights;	
		}
	}
	else{
		cfg.dv_spectrum_depths = cfg.spectrum_depths;
		cfg.dv_spectrum_weights = cfg.spectrum_weights;
	}

	//Re-normalise weights
	double sum = 0.0;
	std::vector<double>::iterator it = cfg.dv_spectrum_weights.begin();
	for( ; it != cfg.dv_spectrum_weights.end(); ++it ) sum += *it;
	it = cfg.dv_spectrum_weights.begin();
	for( ; it != cfg.dv_spectrum_weights.end(); ++it ) *it = *it/sum;

	//Set spectrum indexes
	cfg.dv_spectrum_indexes.clear();
	if( se_energy >= 0 ) cfg.dv_spectrum_indexes.push_back(se_energy);
	else{
		for( int i = 0; i < cfg.dv_spectrum_depths.size(); i++ )
			cfg.dv_spectrum_indexes.push_back(i);
	}

}

void initSingleEnergyConfig( config_t &se_cfg, config_t &cfg, int energy ){

	se_cfg = cfg;

	//Adjust the depth parameters to include only one spectrum
	se_cfg.model_depth = cfg.spectrum_depths[energy];	
	se_cfg.spectrum_depths.resize(1);
	se_cfg.spectrum_depths[0] = cfg.spectrum_depths[energy];
	se_cfg.spectrum_weights.resize(1);
	se_cfg.spectrum_weights[0] = 1.0;

	//Re-derive the derived parameters
	initDerivedConfig(se_cfg, energy);

}