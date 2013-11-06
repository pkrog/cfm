/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# main.cpp
#
# Description: 	Identify the most likely candidate structure for a given
#				spectrum from a list of candidates.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

int main(int argc, char *argv[]);

#include "FragmentGraphGenerator.h"
#include "MolData.h"
#include "Identifier.h"
#include "Comparators.h"

#include <iostream>
#include <fstream>
#include <string>

void readInCandidates( std::vector<Candidate> &candidates, std::string &candidate_file);
void reportResults( std::vector<Candidate> &candidates, std::ostream &out );

int main(int argc, char *argv[])
{
	bool to_stdout = true;
	std::string output_filename;
	std::string param_filename = "param_output.log";
	std::string config_filename = "param_config.txt";
	int num_highest = -1;
	double abs_mass_tol = 0.01, ppm_mass_tol = 10.0;
	double prob_thresh_for_prune = 0.001;
	std::string score_type = "Jaccard";

	if (argc < 3 || argc > 11 )
	{
		std::cout << std::endl << std::endl;
		std::cout << std::endl << "Usage: cfm-id.exe <spectrum_file> <candidate_file> <num_highest> <ppm_mass_tol> <abs_mass_tol> <param_filename> <config_filename> <output_filename>" << std::endl << std::endl << std::endl;
		std::cout << std::endl << "spectrum_file:" << std::endl << "The filename where the input spectra can be found as a list of peaks 'mass intensity' delimited by lines, with either 'low','med' and 'high' lines beginning spectra of different energy levels, or 'energy0', 'energy1', etc. ";
		std::cout << "e.g." << std::endl << "energy0" << std::endl << "65.02 40.0" << std::endl << "86.11 60.0" << std::endl << "energy1" << std::endl << "65.02 100.0 ... etc" << std::endl;
		std::cout << std::endl << "candidate_file:" << std::endl << "The filename where the input list of candidate structures can be found - line separated 'id smiles_or_inchi' pairs." << std::endl;
		std::cout << std::endl << "num_highest (opt):" << std::endl << "The number of (ranked) candidates to return or -1 for all (if not given, returns all in ranked order)" << std::endl;
		std::cout << std::endl << "ppm_mass_tol (opt):" << std::endl << "The mass tolerance in ppm to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs (if not given defaults to 10ppm)" << std::endl;
		std::cout << std::endl << "abs_mass_tol (opt):" << std::endl << "The mass tolerance in abs Da to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs ( if not given defaults to 0.01Da)" << std::endl;
		std::cout << std::endl << "prob_thresh_for_prune (opt):" << std::endl << "The probabiltiy threshold at which to prune unlikely fragmnetations (default 0.001)" << std::endl;
		std::cout << std::endl << "param_filename (opt):" << std::endl << "The filename where the parameters of a trained cfm model can be found (if not given, assumes param_output.log in current directory)" << std::endl;
		std::cout << std::endl << "config_filename (opt):" << std::endl << "The filename where the configuration parameters of the cfm model can be found (if not given, assumes param_config.txt in current directory)" << std::endl;
		std::cout << std::endl << "score_type (opt):" << std::endl << "The type of scoring function to use when comparing spectra. Options: Jaccard (default), DotProduct" << std::endl;
		std::cout << std::endl << "output_filename (opt):" << std::endl << "The filename of the output file to write to (if not given, prints to stdout)" << std::endl;
		exit(1);
	}

	std::string spectrum_file = argv[1];
	std::string candidate_file = argv[2];
	if( argc > 3 ) num_highest = atoi(argv[3]);
	if( argc > 4 ) ppm_mass_tol = atof(argv[4]);
	if( argc > 5 ) abs_mass_tol = atof(argv[5]);
	if( argc > 6 ) prob_thresh_for_prune = atof(argv[6]);
	if( argc > 7 ) param_filename = argv[7];
	if( argc > 8 ) config_filename = argv[8];
	if( argc > 9 ) score_type = argv[9];
	if( argc > 10 ){
		output_filename = argv[10];
		to_stdout = false;
	}

	//Read in the input spectrum
	MolData targetData( "Unknown", "Unknown" );
	targetData.readInSpectraFromFile( spectrum_file );

	//Initialise model configuration
	config_t cfg;
	initConfig( cfg, config_filename );

	//Read in the parameters
	Param param( param_filename );

	//Fetch the list of candidates
	std::vector<Candidate> candidates;
	readInCandidates( candidates, candidate_file );

	//Set up the comparator and identifier
	Comparator *cmp;
	if( score_type == "DotProduct" ){ 
		std::cout << "Using DotProduct score function" << std::endl;
		cmp = new DotProduct( ppm_mass_tol, abs_mass_tol );
	}else if( score_type == "Jaccard" ){ 
		std::cout << "Using Jaccard score function" << std::endl;
		cmp = new Jaccard( ppm_mass_tol, abs_mass_tol );
	}else{ 
		std::cout << "Unknown comparator type: " << score_type << std::endl;
		exit(1);
	}
	Identifier identifier( &param, &cfg, cmp, prob_thresh_for_prune );

	//Rank the candidates
	identifier.rankCandidatesForSpecMatch( candidates, targetData.getSpectra() );

	//Keep only the num_highest
	if( num_highest > 0 && num_highest < (int)candidates.size() ) candidates.resize( num_highest );

	//Set up the output (to file or stdout)
	std::streambuf * buf;
	std::ofstream of;
	if( !to_stdout ) {
		of.open(output_filename.c_str());
		buf = of.rdbuf();
	} else buf = std::cout.rdbuf();
	std::ostream out(buf);

	//Report the results
	reportResults( candidates, out );
	
	delete cmp;

	return(0);    
}

void readInCandidates( std::vector<Candidate> &candidates, std::string &candidate_file){

	std::string line, smiles_or_inchi, id;
	std::ifstream ifs ( candidate_file.c_str() , std::ifstream::in );

	if( !ifs.good() ){ 
		std::cout << "Could not open input file " << candidate_file << std::endl;
		return;
	}

	while( ifs.good() ){

		getline( ifs, line );
		if( line.size() < 3 ) continue;

		std::stringstream ss(line);
		ss >> id >> smiles_or_inchi;

		candidates.push_back( Candidate( id, smiles_or_inchi ) );
	}

}

void reportResults( std::vector<Candidate> &candidates, std::ostream &out ){

	std::vector<Candidate>::iterator it = candidates.begin();
	out << std::setprecision(8);
	for( int rank = 1; it != candidates.end(); ++it, rank++ ){
		out << rank << " ";
		out << it->getScore() << " ";
		out << *it->getId() << " ";
		out << *it->getSmilesOrInchi() << std::endl;
	}

}