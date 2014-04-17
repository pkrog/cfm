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

int main(int argc, char *argv[]);

#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char *argv[])
{
	bool to_stdout = true;
	int do_annotate = 0;
	std::string output_filename;
	std::string param_filename = "param_output.log";
	std::string config_filename = "param_config.txt";
	double prob_thresh_for_prune = 0.001;

	if (argc != 6 && argc != 2 && argc != 5 && argc != 3 && argc != 7)
	{
		std::cout << std::endl << std::endl;
		std::cout << std::endl << "Usage: cfm-predict.exe <input_smiles_or_inchi> <prob_thresh_for_prune> <param_filename> <config_filename> <include_annotations> <output_filename>" << std::endl << std::endl << std::endl;
		std::cout << std::endl << "input_smiles_or_inchi:" << std::endl << "The smiles or inchi string of the structure whose spectra you want to predict" << std::endl;
		std::cout << std::endl << "prob_thresh_for_prune (opt):" << std::endl << "The probability below which to prune unlikely fragmentations (default 0.001)" << std::endl;
		std::cout << std::endl << "param_filename (opt):" << std::endl << "The filename where the parameters of a trained cfm model can be found (if not given, assumes param_output.log in current directory)" << std::endl;
		std::cout << std::endl << "config_filename (opt):" << std::endl << "The filename where the configuration parameters of the cfm model can be found (if not given, assumes param_config.txt in current directory)" << std::endl;
		std::cout << std::endl << "include_annotations (opt):" << std::endl << "Whether to include fragment information in the output spectra (0 = NO (default), 1 = YES )" << std::endl;
		std::cout << std::endl << "output_filename (opt):" << std::endl << "The filename of the output spectra file to write to (if not given, prints to stdout)" << std::endl;
		exit(1);
	}
	
	std::string input_smiles_or_inchi = argv[1];
	if( argc >= 3 ) prob_thresh_for_prune = atof(argv[2]);
	if( argc >= 5 ){
		param_filename = argv[3];
		config_filename = argv[4];
	}
	if( argc >= 6 ) do_annotate = atoi(argv[5]);
	if( argc == 7){
		output_filename = argv[6];
		to_stdout = false;
	}

	//Initialise model configuration
	config_t cfg;
	initConfig( cfg, config_filename );

	//Read in the parameters
	Param param( param_filename );

	//Create the MolData structure with the input
 	MolData moldata( "NullId", input_smiles_or_inchi.c_str() );

	//Calculate the pruned FragmentGraph
	moldata.computeLikelyFragmentGraphAndSetThetas(&param, &cfg, prob_thresh_for_prune, do_annotate);

	//Predict the spectra (and post-process, use existing thetas)
	moldata.computePredictedSpectra( param, cfg, true, true );

	//Set up the output stream
	std::streambuf * buf;
	std::ofstream of;
	if( !to_stdout ) {
		of.open(output_filename.c_str());
		if( !of.is_open() ) 
			std::cout << "Warning: Trouble opening output file" << std::endl;
		buf = of.rdbuf();
	} else buf = std::cout.rdbuf();
	std::ostream out(buf);

	//Write the spectra to output
	moldata.outputSpectra(out, "Predicted", do_annotate);
	out << std::endl;
	if( do_annotate ) moldata.getFragmentGraph()->writeFragmentsOnly(out);
	if( !to_stdout ) of.close();

	return(0);    
}
