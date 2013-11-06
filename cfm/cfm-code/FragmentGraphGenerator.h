/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentGraphGenerator.h
#
# Description: 	FragmentGraphGenerator class for generating a fragment tree.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __FRAG_GEN_H__
#define __FRAG_GEN_H__

#include "FragmentGraph.h"
#include "FragmentTreeNode.h"
#include "Features.h"
#include "Param.h"
#include "Config.h"

//Class for generating a FragmentGraph given a starting node in the
//enumerated fragmentation tree and a depth
class FragmentGraphGenerator{
public:
	//Constructor
	FragmentGraphGenerator() : verbose(0), mols_to_fv(false) {};
	FragmentGraphGenerator(int a_verbose) : verbose(a_verbose), mols_to_fv(false) {};

	//Constructor to use if the transition molecules are to be replaced by a feature vector as
	//soon as they are created (for less memory usage).
	FragmentGraphGenerator(FeatureCalculator *a_fc ) :
		fc(a_fc), verbose(0), mols_to_fv(true) {};

	//Start a graph. Compute can then add to this graph, but it is the caller's 
	//responsibility to delete it
	FragmentGraph *createNewGraph();

	//Create the starting node from a smiles or inchi string - responsibility of caller to delete
	FragmentTreeNode *createStartNode( std::string &smiles_or_inchi );

	//Compute a FragmentGraph starting at the given node and computing to the depth given.
	//The output will be appended to the current_graph
	void compute( FragmentTreeNode &startnode, int remaining_depth, int parentid );

protected:

	FeatureCalculator *fc;
	FragmentGraph *current_graph;
	bool mols_to_fv;
	int verbose;

private:

	//Record of previous computations so we know to what depth each fragment has been computed
	std::map<int, int> id_depth_computed_cache;	
	//Helper function - check if the fragment has already been computed to at least this depth
	int alreadyComputed(int id, int remaining_depth);
	//Helper function - count the number of spare electron pairs
	int countExtraElectronPairs( RDKit::RWMol *rwmol );

	//Helper functions - used to create labels on atoms and bonds, 
	//that will be used in Feature Calculations and can't be computed once
	//the molecule is broken
	void labelGasteigers( RDKit::RWMol *rwmol );
	void labelAromatics( RDKit::RWMol *rwmol );
	void labelOriginalMasses( RDKit::RWMol *rwmol );
	void initialiseRoots( RDKit::RWMol *rwmol );
};

//Class to be used if a graph is to be pruned as it's created 
//removing anything with probability below a given threshold
class LikelyFragmentGraphGenerator : public FragmentGraphGenerator {
public:
	//Constructor
	LikelyFragmentGraphGenerator(Param *a_param, config_t *a_cfg, double a_prob_thresh ) :
	  cfg( a_cfg ), param(a_param), log_prob_thresh( log(a_prob_thresh) ) { fc = new FeatureCalculator(*param->getFeatureNames());  mols_to_fv = true; };

	//Start a graph. Compute can then add to this graph, but it is the caller's 
	//responsibility to delete it
	FragmentGraph *createNewGraph();

	//Compute a FragmentGraph starting at the given node and computing to the depth given.
	//The output will be appended to the current_graph
	void compute( FragmentTreeNode &startnode, int remaining_depth, int parentid, double parent_log_prob );

private:
	Param *param;
	config_t *cfg;
	double log_prob_thresh;

	//Record of previous computations so we know to what probability each fragment has been computed at
	std::map<int, double> id_prob_computed_cache;	
	int alreadyComputed(int id, double prob_offset);

};

#endif // __FRAG_GEN_H__