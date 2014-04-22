/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# em_test.cpp
#
# Description: Test code for EM.cpp
#
# Author: Felicity Allen
# Created: August 2013
#########################################################################*/

#include "mpi.h"

#include "em_test.h"

EMTestSelfProduction::EMTestSelfProduction(){
	description = "Test integrated CE-CFM EM ability to learn single spectrum";
}

void EMTestSelfProduction::runTest(){
	
	bool pass = true; 
	double intensity_tol = 2.0; //For intensity in range 0-100
    
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump != 1 ){
		std::cout << "Error: Test intended for one processor only" << std::endl;
		passed = false;
		return;
	}

	//Config
	config_t cfg;
	std::string param_cfg_file = "tests/test_data/example_param_config.txt";
	initConfig( cfg, param_cfg_file );
	cfg.lambda = 0.0000001;

	for( int config_state = 0; config_state <= 3; config_state++ ){

		if( config_state == 1 ){
			std::cout << "Testing Depth 1-2-3 Configuration" << std::endl;
			cfg.spectrum_depths[0] = 1; cfg.spectrum_depths[1] = 2;  cfg.spectrum_depths[2] = 3;
			cfg.model_depth = 3;
			initDerivedConfig(cfg);
		}else std::cout << "Testing Depth 2-4-6 Configuration" << std::endl;

		//Feature Calculator
		std::string feature_cfg_file = "tests/test_data/example_feature_config_withquadratic.txt";
		FeatureCalculator fc( feature_cfg_file );

		//Prepare some simple data
		std::vector<MolData> data;
		std::string id = "TestMol", smiles = "NCCCN";
		std::string spec_file = "tests/test_data/example_spectra.txt";
		data.push_back( MolData( id, smiles, 0 ) );
		data[0].setIonizationMode(false);
		data[0].computeFragmentGraphAndReplaceMolsWithFVs(cfg.fg_depth, &fc);
		data[0].readInSpectraFromFile( spec_file );

		//Run EM
		std::string status_file = "tmp_status_file.log";
		EM em( &cfg, &fc, status_file );
		double Q = em.run( data, 1 );
		std::string param_filename = "tmp_param_output.log";
		em.writeParamsToFile( param_filename );

		//Predict the output spectra
		Param param( param_filename );
		data[0].computePredictedSpectra( param, cfg );

		//Compare the original and predicted spectra - should be able to overfit
		//very close to the actual values since training on same (and only same) mol
		for( unsigned int energy = 0; energy < data[0].getNumSpectra(); energy++ ){
			const Spectrum *orig_spec = data[0].getSpectrum(energy);
			const Spectrum *predicted_spec = data[0].getPredictedSpectrum(energy);
			Spectrum::const_iterator ito = orig_spec->begin();
			for( ; ito != orig_spec->end(); ++ito ){
			
				//Find a peak in the predicted spectrum with the same mass
				Spectrum::const_iterator itp = predicted_spec->begin();
				int found = 0;
				for( ; itp != predicted_spec->end(); ++itp ){
					double mass_tol = getMassTol( cfg.abs_mass_tol, cfg.ppm_mass_tol, itp->mass );
					if( fabs(itp->mass - ito->mass ) < mass_tol ){
					
						//Check the intensity values of the matching peaks
						if( fabs(itp->intensity - ito->intensity) > intensity_tol ){
							std::cout << "Mismatch in predicted peak intensity for mass " << ito->mass << ": Expecting " << ito->intensity << " but found " << itp->intensity << std::endl;
							pass = false;
						}
						found = 1; 
						break;
					}
				}
				if( !found ){
					std::cout << "Could not find matching predicted peak for mass " << ito->mass << std::endl;
					pass = false;
				}
			}
		}
	}

	passed = pass;

}

EMTestSingleEnergySelfProduction::EMTestSingleEnergySelfProduction(){
	description = "Test integrated SE-CFM EM ability to learn single spectrum";
}

void EMTestSingleEnergySelfProduction::runTest(){
	
	bool pass = true; 
	double intensity_tol = 2.0; //For intensity in range 0-100
    
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump != 1 ){
		std::cout << "Error: Test intended for one processor only" << std::endl;
		passed = false;
		return;
	}

	//Config
	config_t orig_cfg;
	std::string param_cfg_file = "tests/test_data/example_param_config.txt";
	initConfig( orig_cfg, param_cfg_file );
	orig_cfg.lambda = 0.0000001;
	orig_cfg.use_single_energy_cfm = 1;
	orig_cfg.spectrum_depths[1] = 2;
	orig_cfg.spectrum_depths[2] = 2;

	//Feature Calculator
	std::string feature_cfg_file = "tests/test_data/example_feature_config_withquadratic.txt";
	FeatureCalculator fc( feature_cfg_file );

	//Prepare some simple data
	std::vector<MolData> data;
	std::string id = "TestMol", smiles = "NCCCN";
	std::string spec_file = "tests/test_data/example_spectra.txt";
	data.push_back( MolData( id, smiles, 0 ) );
	data[0].setIonizationMode(false);
	data[0].computeFragmentGraphAndReplaceMolsWithFVs(orig_cfg.fg_depth, &fc);
	data[0].readInSpectraFromFile( spec_file );

	Param *final_params;
	for( int energy = 0; energy < 3; energy++ ){

		config_t cfg;
		initSingleEnergyConfig( cfg, orig_cfg, energy );

		//Run EM
		std::string status_file = "tmp_status_file.log";
		EM em( &cfg, &fc, status_file );
		double Q = em.run( data, 1 );
		std::string param_filename = "tmp_param_output.log";
		em.writeParamsToFile( param_filename );

		if(energy == 0 ) final_params = new Param(param_filename);
		else{ 
			Param eparam(param_filename);
			final_params->appendNextEnergyParams( eparam, energy );
		}
	}

	//Predict the output spectra
	data[0].computePredictedSpectra( *final_params, orig_cfg );

	//Compare the original and predicted spectra - should be able to overfit
	//very close to the actual values since training on same (and only same) mol
	for( unsigned int energy = 0; energy < data[0].getNumSpectra(); energy++ ){
		const Spectrum *orig_spec = data[0].getSpectrum(energy);
		const Spectrum *predicted_spec = data[0].getPredictedSpectrum(energy);
		Spectrum::const_iterator ito = orig_spec->begin();
		for( ; ito != orig_spec->end(); ++ito ){
			
			//Find a peak in the predicted spectrum with the same mass
			Spectrum::const_iterator itp = predicted_spec->begin();
			int found = 0;
			for( ; itp != predicted_spec->end(); ++itp ){
				double mass_tol = getMassTol( orig_cfg.abs_mass_tol, orig_cfg.ppm_mass_tol, itp->mass );
				if( fabs(itp->mass - ito->mass ) < mass_tol ){
					
					//Check the intensity values of the matching peaks
					if( fabs(itp->intensity - ito->intensity) > intensity_tol ){
						std::cout << "Mismatch in predicted peak intensity for mass " << ito->mass << ": Expecting " << ito->intensity << " but found " << itp->intensity << std::endl;
						pass = false;
					}
					found = 1; 
					break;
				}
			}
			if( !found ){
				std::cout << "Could not find matching predicted peak for mass " << ito->mass << std::endl;
				pass = false;
			}
		}
	}
	passed = pass;
}

EMTestMultiProcessor::EMTestMultiProcessor(){
	description = "Test EM running on multiple processors";
}

void EMTestMultiProcessor::runTest(){
	
	bool pass = true;
	double tol = 1e-6;
	
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump <= 2 ){
		std::cout << "Error: Test intended for multiple processors" << std::endl;
		passed = false;
		return;
	}

	//Initialisation
	config_t cfg;
	std::string cfg_file = "tests/test_data/example_param_config.txt";
	initConfig( cfg, cfg_file );
	std::string fc_file = "tests/test_data/example_feature_config_withquadratic.txt";
	FeatureCalculator fc( fc_file );


	std::string id1 = "TestMol1", id2 = "TestMol2", id3 = "TestMol3";
	std::string smiles1 = "NCCCN", smiles2 = "N=CC(OC)CN", smiles3 = "N(CCCC)CCCN";
	std::string spec_file = "tests/test_data/example_spectra.txt";

	//First Run three molecules on the Master Only
	std::cout << "Running on master only..." << std::endl;
	std::vector<MolData> data;
	if( mpi_rank == MASTER ){
		data.push_back( MolData( id1, smiles1, 0 ) );
		data.push_back( MolData( id2, smiles2, 0 ) );
		data.push_back( MolData( id3, smiles3, 0 ) );
		for( int i = 0 ; i < 3; i++ ){
			data[i].computeFragmentGraphAndReplaceMolsWithFVs(cfg.fg_depth, &fc);
			data[i].readInSpectraFromFile( spec_file );
		}
	}
	std::string status_file = "tmp_status_file.log";
	EM em1( &cfg, &fc, status_file );
	double all_master_Q = em1.run( data, 1, true );
	std::cout << "Q = " << all_master_Q << std::endl;

	//Now run the same molecules on 3 separate processes
	std::cout << "Running on 3 processes..." << std::endl;
	data.clear();
	if( mpi_rank == MASTER ) data.push_back( MolData( id1, smiles1, 0 ) );
	else if( mpi_rank == 1 ) data.push_back( MolData( id2, smiles2, 0 ) );
	else if( mpi_rank == 2 ) data.push_back( MolData( id3, smiles3, 0 ) );
	data[0].computeFragmentGraph(cfg.fg_depth);
	data[0].computeFeatureVectors(&fc);
	data[0].readInSpectraFromFile( spec_file );
	EM em2( &cfg, &fc, status_file );
	double separated_Q = em2.run( data, 1, true );
	std::cout << "Q = " << separated_Q << std::endl;

	//Check that the output Q values are the same
	std::cout << "Q1: " << all_master_Q << " Q2: " << separated_Q << std::endl;
	if( fabs( all_master_Q - separated_Q ) > tol ){
		std::cout << "Mismatch in Q values between running on one vs many processes" << std::endl;
		pass = false;
	}

	passed = pass;

}
