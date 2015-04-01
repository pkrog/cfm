/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# main.cpp
#
# Description: 	Predict the MS/MS spectra for a given structure using a
#				pre-trained CFM model.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Config.h"
#include "Param.h"
#include "Features.h"
#include "MolData.h"

#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/filesystem.hpp>

int main(int argc, char *argv[]);

#include <iostream>
#include <fstream>
#include <string>

class SpectrumPredictionException: public std::exception{
private:
    std::string message_;
public:
	SpectrumPredictionException(const std::string& message) throw() : message_(message) {};
	virtual const char* what() const throw(){
		std::cout << message_ << std::endl;
		return message_.c_str();
	}
	~SpectrumPredictionException() throw() {};
};

class FileException: public std::exception{
private:
    std::string message_;
public:
	FileException(const std::string& message) throw() : message_(message) {};
	virtual const char* what() const throw(){
		std::cout << message_ << std::endl;
		return message_.c_str();
	}
	~FileException() throw() {};
};


void parseInputFile(std::vector<MolData> &data, std::string &input_filename );

int main(int argc, char *argv[])
{
	bool to_stdout = true;
	int do_annotate = 0;
	int apply_postprocessing = 1;
	std::string output_filename;
	std::string param_filename = "param_output.log";
	std::string config_filename = "param_config.txt";
	double prob_thresh_for_prune = 0.001;

	if (argc != 6 && argc != 2 && argc != 5 && argc != 3 && argc != 7 && argc != 8)
	{
		std::cout << std::endl << std::endl;
		std::cout << std::endl << "Usage: cfm-predict.exe <input_smiles_or_inchi> <prob_thresh_for_prune> <param_filename> <config_filename> <include_annotations> <output_filename> <apply_post_processing>" << std::endl << std::endl << std::endl;
		std::cout << std::endl << "input_smiles_or_inchi_or_file:" << std::endl << "The smiles or inchi string of the structure whose spectra you want to predict, or a .txt file containing a list of <id smiles> pairs, one per line." << std::endl;
		std::cout << std::endl << "prob_thresh_for_prune (opt):" << std::endl << "The probability below which to prune unlikely fragmentations (default 0.001)" << std::endl;
		std::cout << std::endl << "param_filename (opt):" << std::endl << "The filename where the parameters of a trained cfm model can be found (if not given, assumes param_output.log in current directory)" << std::endl;
		std::cout << std::endl << "config_filename (opt):" << std::endl << "The filename where the configuration parameters of the cfm model can be found (if not given, assumes param_config.txt in current directory)" << std::endl;
		std::cout << std::endl << "include_annotations (opt):" << std::endl << "Whether to include fragment information in the output spectra (0 = NO (default), 1 = YES ). Note: ignored for msp/mgf output." << std::endl;
		std::cout << std::endl << "output_filename_or_dir (opt):" << std::endl << "The filename of the output spectra file to write to (if not given, prints to stdout), OR directory if multiple smiles inputs are given (else current directory) OR msp or mgf file." << std::endl;
		std::cout << std::endl << "apply_postprocessing (opt):" << std::endl << "Whether or not to post-process predicted spectra to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks (whichever comes first) (0 = OFF, 1 = ON (default))." << std::endl;
		exit(1);
	}
	
	std::string input_smiles_or_inchi = argv[1];
	if( argc >= 3 ) prob_thresh_for_prune = atof(argv[2]);
	if( argc >= 5 ){
		param_filename = argv[3];
		config_filename = argv[4];
	}
	if( argc >= 6 ) do_annotate = atoi(argv[5]);
	if( argc >= 7){
		output_filename = argv[6];
		to_stdout = false;
	}
	if( argc == 8 ) apply_postprocessing = atoi(argv[7]);

	//Initialise model configuration
	config_t cfg;
	if( !boost::filesystem::exists( config_filename ) ){
		std::cout << "Could not find file: " <<  config_filename << std::endl;
		throw FileException("Could not find file: " + config_filename);
	}
	initConfig( cfg, config_filename );

	//Read in the parameters
	if( !boost::filesystem::exists( param_filename ) ){
		std::cout << "Could not find file: " <<  param_filename << std::endl;
		throw FileException("Could not find file: " + param_filename);
	}
	Param param( param_filename );

	//Check for mgf or msp output - and setup in exists
	int output_mode = NO_OUTPUT_MODE;
	std::ostream *out;
	std::ofstream of; std::streambuf * buf;
	if( !to_stdout && output_filename.substr( output_filename.size() - 4 ).compare(".msp") == 0 ){
		output_mode = MSP_OUTPUT_MODE;
		of.open( output_filename.c_str() );
		if( !of.is_open() ){ std::cout << "Error: Could not open output msp file " << output_filename << std::endl; throw FileException("Could not open output msp file " + output_filename); }
		buf = of.rdbuf();
		out = new std::ostream(buf);
	}
	else if( !to_stdout && output_filename.substr( output_filename.size() - 4 ).compare(".mgf") == 0 ){
		output_mode = MGF_OUTPUT_MODE;
		of.open( output_filename.c_str() );
		if( !of.is_open() ){ std::cout << "Error: Could not open output mgf file " << output_filename << std::endl; throw FileException("Could not open output mgf file " + output_filename); }
		buf = of.rdbuf();
		out = new std::ostream(buf);
	}

	//Check for batch input - if found, read in inchis and set up output directory, mgf or msp
	std::vector<MolData> data;
	bool batch_run = false;
	std::string output_dir_str = "";
	if( input_smiles_or_inchi.substr( input_smiles_or_inchi.size() - 4 ).compare(".txt") == 0 ){
		parseInputFile( data, input_smiles_or_inchi );
		batch_run = true;
		if( !to_stdout && output_mode == NO_OUTPUT_MODE ){
			if(  output_filename != "." && !boost::filesystem::exists(output_filename) )
				boost::filesystem::create_directory(output_filename);
			output_dir_str = output_filename + "/";
		}
		to_stdout = false;
	}
	else data.push_back( MolData( "NullId", input_smiles_or_inchi.c_str() ) );

	std::vector<MolData>::iterator it = data.begin();
	for( ; it != data.end(); ++it ){
		//Create the MolData structure with the input
		it->setIonizationMode(cfg.ionization_mode == NEGATIVE_IONIZATION_MODE);	
		try{
			//Calculate the pruned FragmentGraph
			it->computeLikelyFragmentGraphAndSetThetas(&param, &cfg, prob_thresh_for_prune, do_annotate);

			//Predict the spectra (and post-process, use existing thetas)
			it->computePredictedSpectra( param, cfg, apply_postprocessing, true );
		}
		catch( RDKit::MolSanitizeException e ){
			std::cout << "Could not sanitize input: " << it->getSmilesOrInchi() << std::endl;
			if(!batch_run) throw SpectrumPredictionException("RDKit could not sanitize input: " + it->getSmilesOrInchi());
			continue;
		}
		catch( RDKit::SmilesParseException pe ){
			std::cout << "Could not parse input: " << it->getSmilesOrInchi() << std::endl;
			if(!batch_run) throw SpectrumPredictionException("RDKit could not parse input: " + it->getSmilesOrInchi());
			continue;
		}
		catch( FragmentGraphGenerationException ){
			std::cout << "Could not compute fragmentation graph for input: " << it->getSmilesOrInchi() << std::endl;
			if(!batch_run) throw SpectrumPredictionException("Could not compute fragmentation graph for input: " + it->getSmilesOrInchi());
			continue;
		}

		//Set up the output stream (if not already set up)
		if( output_mode == NO_OUTPUT_MODE ){
			if( !to_stdout ) {
				if( batch_run ) output_filename = output_dir_str + it->getId() + ".log";
				of.open( output_filename.c_str() );
				if( !of.is_open() ){ std::cout << "Error: Could not open output file " << output_filename << std::endl; throw FileException("Could not open output file " + output_filename ); }
				buf = of.rdbuf();
			} else buf = std::cout.rdbuf();
			out = new std::ostream(buf);
		}

		//Write the spectra to output
		if( output_mode == NO_OUTPUT_MODE ){
			it->outputSpectra(*out, "Predicted", do_annotate);
			*out << std::endl;
			if( do_annotate ) it->getFragmentGraph()->writeFragmentsOnly(*out);
		}
		else if( output_mode == MSP_OUTPUT_MODE ) it->writePredictedSpectraToMspFileStream( *out );
		else if( output_mode == MGF_OUTPUT_MODE ) it->writePredictedSpectraToMgfFileStream( *out );

		if( output_mode == NO_OUTPUT_MODE ){
			if( !to_stdout ) of.close();
			delete out;
		}
	}
	if( output_mode != NO_OUTPUT_MODE ) delete out;
	return(0);    
}


void parseInputFile(std::vector<MolData> &data, std::string &input_filename ){

	std::string line, smiles_or_inchi, id;
	std::ifstream ifs ( input_filename.c_str() , std::ifstream::in );

	//Get the first line - the number of input molecules
	if( !ifs.good() ){
		std::cout << "Could not open input file: " << input_filename << std::endl;
		throw FileException("Could not open input file: " + input_filename);
	}

	//Now get all the molecules
	int i = 0;
	while( ifs.good()){
		i++;

		getline( ifs, line );
		if( line.size() < 3 ) continue;

		std::stringstream ss(line);
		ss >> id >> smiles_or_inchi;

		data.push_back( MolData( id.c_str(), smiles_or_inchi.c_str() ) );
	}

	std::cout << "Read " << data.size() << " molecules from input file." << std::endl;

}
