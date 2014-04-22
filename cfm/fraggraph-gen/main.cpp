/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# main.c
#
# Description: 	Search for valid fragment configurations using backtrack.
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

#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char *argv[])
{
	std::string output_filename;
	std::string smiles_or_inchi; 
	int max_depth;

	bool verbose = false;
	bool to_stdout = true;
	bool fullgraph = true;
    
	if (argc != 4 && argc != 5 && argc != 6)
	{
		std::cout << std::endl << "Usage:" << std::endl;
		std::cout << "fraggraph-gen.exe <smiles or inchi string> <max_depth> <mode(+/-)> ";
		std::cout << "<opt: fullgraph (default) or fragonly> <opt: output_filename (else stdout)>" << std::endl << std::endl ;
		exit(1);
	}
    smiles_or_inchi = argv[1];
    max_depth = atoi(argv[2]);
	int negative_mode = -1;
	if( argv[3][0] == '+' ) negative_mode = 0;
	else if(argv[3][0] == '-' ) negative_mode = 1;
	else{ 
		std::cout << "Unknown setting for ionization mode: expecting + or -" << std::endl;
		exit(1);
	}

	if( argc > 4 ){
		std::string type = argv[4];
		if( type == "fullgraph" ) fullgraph = true;
		else if( type == "fragonly" ) fullgraph = false;
		else std::cout << "Unknown setting for output type: " << type << " (expecting fullgraph or fragonly)" << std::endl;
	}
	if( argc > 5 ){ 
		output_filename = argv[5];
		to_stdout = false;
	}
	
	//Set up the output (to file or stdout)
	std::streambuf * buf;
	std::ofstream of;
	if( !to_stdout ) {
		of.open(output_filename.c_str());
		buf = of.rdbuf();
	} else buf = std::cout.rdbuf();
	std::ostream out(buf);

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, negative_mode);
	FragmentGraph *graph = gg.createNewGraph();
	gg.compute( *startNode, max_depth, -1 );
	delete startNode;

	//Write to output
	if( fullgraph ) graph->writeFullGraph( out );
	else graph->writeFragmentsOnly( out );
	delete graph;

	return(0);    
}
