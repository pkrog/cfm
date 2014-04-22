/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# main.cpp
#
# Description: 	Learn the parameters of a model for the mass spec 
#				fragmentation process using EM, then use it to predict
#				the mass spectra.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "mpi.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#include <fstream>
#include <sstream>
#include <iostream>

#include "EM.h"
#include "Features.h"
#include "Config.h"
#include "MolData.h"
#include "FragmentGraphGenerator.h"

void parseInputFile(std::vector<MolData> &data, std::string &input_filename, int mpi_rank, int mpi_nump );
void trainCombinedEnergyCFM( std::string &param_filename, config_t &cfg, FeatureCalculator &fc, std::string &status_filename, int group, std::vector<MolData> &data);
void trainSingleEnergyCFM( std::string &param_filename, config_t &cfg, FeatureCalculator &fc, std::string &status_filename, int group, std::vector<MolData> &data, int start_energy, int no_train);

int main(int argc, char *argv[])
{
	int mpi_rank, mpi_nump;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if (argc < 4 || argc > 9)
	{
		std::cout << std::endl << std::endl;
		std::cout << std::endl << "Usage: cfm-train.exe <input_filename> <feature_filename> <config_filename> <peakfile_dir> <group> <status_filename> <no_train> <start_energy>" << std::endl << std::endl << std::endl;
		std::cout << std::endl << "input_filename:" << std::endl << "Text file with number of mols on first line, then " << std::endl << "id smiles_or_inchi cross_validation_group" << std::endl << "on each line after that." << std::endl;
		std::cout << std::endl << "feature_filename:" << std::endl << "Text file with list of feature names to include, line separated:" << std::endl << "BreakAtomPair" << std::endl << "IonRootPairs...etc" << std::endl;
		std::cout << std::endl<< "config_filename:" << std::endl << "Text file listing configuration parameters. Line separated 'name value'." << std::endl;
		std::cout << std::endl << "peakfile_dir:" << std::endl << "Directory containing files with spectra. Each file should be called <id>.txt, where <id> is the id specified in the input file, and contains a list of peaks 'mass intensity' on each line, with either 'low','med' and 'high' lines beginning spectra of different energy levels, or 'energy0', 'energy1', etc. ";
		std::cout << "e.g." << std::endl << "energy0" << std::endl << "65.02 40.0" << std::endl << "86.11 60.0" << std::endl << "energy1" << std::endl << "65.02 100.0 ... etc" << std::endl;
		std::cout << std::endl << "group (opt):" << std::endl << "Cross validation group to run. Otherwise will assume 10 groups and run all of them." << std::endl;
		std::cout << std::endl << "status_filename (opt):" << std::endl << "Name of file to write logging information as the program runs. If not specified will write to status.log<group>, or status.log if no group is specified" << std::endl;
		std::cout << std::endl << "no_train (opt):" << std::endl << "Set to 1 if the training part should be skipped (useful in debugging - default 0)" << std::endl;
		std::cout << std::endl << "start_energy (opt - se only)"  << std::endl << "Set to starting energy if want to start training part way through (single energy only -default 0)" << std::endl;
		std::cout << std::endl;
		exit(1);
	}

    std::string input_filename = argv[1];	//List (one per line): id, smiles_or_inchi, group
	std::string feature_filename = argv[2]; //List of features, line-spaced
	std::string config_filename = argv[3];	//Parameter configuration
	std::string peakfile_dir = argv[4];		//Directory containing the peak files for each molecule (in format <id>.txt)

	//Cross validation groups to process
	int min_group = 0, max_group = 9;
	if( argc >= 6 ){ 
		min_group = atoi(argv[5]);
		max_group = min_group;
	}
	std::string status_filename("status.log");	//Status file to write to
	if( argc >= 7 ) status_filename = argv[6];					
	else if( argc >= 6 ) status_filename = "status.log" + std::string(argv[5]);	//status.log<group>

	int no_train = 0;
	int start_energy = 0;
	if( argc >= 8 ) no_train = atoi(argv[7]);
	if( argc >= 9 ) start_energy = atoi(argv[8]);

	if( mpi_rank == MASTER ){
		//Create the tmp_data directory if it doesn't exist
		if( !boost::filesystem::exists("tmp_data") )
			boost::filesystem::create_directory("tmp_data");
		if( !boost::filesystem::exists("tmp_data/enumerated_output") )
			boost::filesystem::create_directory("tmp_data/enumerated_output");
		if( !boost::filesystem::exists("tmp_data/predicted_output") )
			boost::filesystem::create_directory("tmp_data/predicted_output");
		//Delete the status file if it already exists
		if( boost::filesystem::exists( status_filename ) )
			boost::filesystem::remove_all( status_filename );
	}

	if( mpi_rank == MASTER ) std::cout << "Parsing input file...";
	std::vector<MolData> data;
	parseInputFile( data, input_filename, mpi_rank, mpi_nump );
	if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;

	if( mpi_rank == MASTER ) std::cout << "Initialising Feature Calculator..";
	FeatureCalculator fc( feature_filename );
	if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;

	if( mpi_rank == MASTER ) std::cout << "Initialising Parameter Configuration..";
	config_t cfg;
	initConfig( cfg, config_filename );
	if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;

	if( mpi_rank == MASTER ) std::cout << "Computing fragmentation graphs and features..";
	std::vector<MolData>::iterator mit = data.begin();
	for( ; mit != data.end(); ++mit ){
		mit->setIonizationMode(cfg.ionization_mode == NEGATIVE_IONIZATION_MODE);
		//If we're not training, only load the ones we'll be testing
		if( (mit->getGroup() >= min_group && mit->getGroup() <= max_group) || !no_train )	
			mit->computeFragmentGraphAndReplaceMolsWithFVs(cfg.fg_depth, &fc);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;

	if( mpi_rank == MASTER ) std::cout << "Loading spectra..";
	for( mit = data.begin(); mit != data.end(); ++mit ){ 
		if( (mit->getGroup() >= min_group && mit->getGroup() <= max_group) || !no_train ){
			std::string spec_file = peakfile_dir + "/" + mit->getId() + ".txt";
			mit->readInSpectraFromFile( spec_file );
			mit->removePeaksWithNoFragment( cfg.abs_mass_tol, cfg.ppm_mass_tol );
		}
	}
	if( cfg.interpolate_spectra ) 
		for( mit = data.begin(); mit != data.end(); ++mit ) 
			mit->createInterpolatedSpectra( cfg );
	MPI_Barrier(MPI_COMM_WORLD);  
	if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;

	for( int group = min_group; group <= max_group; group++ ){
		if( mpi_rank == MASTER ) std::cout << "Running EM to train parameters for Group " << group;	

		std::string param_filename = "tmp_data/param_output";
		param_filename += boost::lexical_cast<std::string>(group);
		param_filename += ".log";

		if( cfg.use_single_energy_cfm )
			trainSingleEnergyCFM( param_filename, cfg, fc, status_filename, group, data, start_energy, no_train);
		else if(!no_train)
			trainCombinedEnergyCFM( param_filename, cfg, fc, status_filename, group, data);


		MPI_Barrier(MPI_COMM_WORLD);   	//Wait for all threads
		if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;
	
		if( mpi_rank == MASTER ) std::cout << "Generating Peak Predictions for Group " << group << "..." << std::endl;
		Param param( param_filename );
		for( mit = data.begin(); mit != data.end(); ++mit ){ 
			if( mit->getGroup() != group ) continue;
			
			//Full enumeration spectrum
			std::string enum_filename = "tmp_data/enumerated_output/" + mit->getId() + ".txt";
			mit->writeFullEnumerationSpectrumToFile( enum_filename );

			//Predicted spectrum
			mit->computePredictedSpectra( param, cfg, true );
			std::string spectra_filename = "tmp_data/predicted_output/" + mit->getId() + ".txt";
			mit->writePredictedSpectraToFile(spectra_filename);
		}
		if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;

		MPI_Barrier(MPI_COMM_WORLD);   	//Wait for all threads
	}

	MPI_Barrier(MPI_COMM_WORLD);   	//Wait for all threads
	MPI_Finalize();

	return(0);    
}

void trainCombinedEnergyCFM( std::string &param_filename, config_t &cfg, FeatureCalculator &fc, std::string &status_filename, int group, std::vector<MolData> &data){

	int mpi_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );

	//Run EM multiple times with random restarts, taking the final one with the best Q
	double prev_Q = -1000000000;
	for( int repeat = 0; repeat < cfg.num_em_restarts; repeat++ ){
		EM em( &cfg, &fc, status_filename );
		double Q = em.run( data, group );
		std::string repeat_filename = param_filename + boost::lexical_cast<std::string>(repeat);
		if( mpi_rank == MASTER ) em.writeParamsToFile( repeat_filename );
		if( Q > prev_Q ){ 
			if( mpi_rank == MASTER ){ 
				std::cout << "Found better Q!" << std::endl;
				em.writeParamsToFile(param_filename);
			}
			prev_Q = Q;
		}
	}
}

void trainSingleEnergyCFM( std::string &param_filename, config_t &cfg, FeatureCalculator &fc, std::string &status_filename, int group, std::vector<MolData> &data, int start_energy, int no_train){
	
	int mpi_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
	Param *final_params;

	for( int energy = 0; energy < cfg.spectrum_depths.size(); energy++ ){

		std::string eparam_filename = param_filename + boost::lexical_cast<std::string>(energy) + "e";

		config_t se_cfg;
		initSingleEnergyConfig( se_cfg, cfg, energy );

		if( energy >= start_energy && !no_train ){
			//Run EM multiple times with random restarts, taking the final one with the best Q
			double prev_Q = -1000000000;
			for( int repeat = 0; repeat < se_cfg.num_em_restarts; repeat++ ){
				EM em( &se_cfg, &fc, status_filename );
				double Q = em.run( data, group );
				std::string repeat_filename = eparam_filename + boost::lexical_cast<std::string>(repeat);
				if( mpi_rank == MASTER ) em.writeParamsToFile( repeat_filename );
				if( Q > prev_Q ){ 
					if( mpi_rank == MASTER ){ 
						std::cout << "Found better Q!" << std::endl;
						em.writeParamsToFile(eparam_filename);
					}
					prev_Q = Q;
				}
			}
		}

		if(energy == 0 ) final_params = new Param(eparam_filename);
		else{
			Param eparam(eparam_filename);
			final_params->appendNextEnergyParams( eparam, energy );
		}
	}
	if( mpi_rank == MASTER ) final_params->saveToFile(param_filename);
	delete final_params;
}



void parseInputFile(std::vector<MolData> &data, std::string &input_filename, int mpi_rank, int mpi_nump ){

	std::string line, smiles_or_inchi, id;
	std::ifstream ifs ( input_filename.c_str() , std::ifstream::in );
	int group, num_mols = 0;

	//Get the first line - the number of input molecules
	if( ifs.good() ){ 
		getline( ifs, line );
		num_mols = atoi(line.c_str());
	}
	else{
		std::cout << "Could not open input file " << input_filename << std::endl;
	}

	//Now get all the molecules
	int i = 0;
	while( ifs.good() && i < num_mols){
		i++;

		getline( ifs, line );
		if( line.size() < 3 ) continue;

		std::stringstream ss(line);
		ss >> id >> smiles_or_inchi >> group;

		//Split the data between processors. Only load in data for this
		//processor
		if( (i % mpi_nump) == mpi_rank )
			data.push_back( MolData( id, smiles_or_inchi, group) );
	}

}
