/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentGraphGenerator.cpp
#
# Description: 	FragmentGraphGenerator class for generating a fragment tree.
#				Also contains Break and FragmentTreeNode classes, for use in
#				this process.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "FragmentGraphGenerator.h"
#include "MILP.h"

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/PeriodicTable.h>
#include <INCHI-API/inchi.h>


//Start a graph. Compute can then add to this graph, but it is the caller's 
//responsibility to delete it
FragmentGraph *FragmentGraphGenerator::createNewGraph(){
	current_graph = new FragmentGraph();
	id_depth_computed_cache.clear();	//The graph is empty, so clear all computation records
	return current_graph;
}

//Start a graph. Compute can then add to this graph, but it is the caller's 
//responsibility to delete it
FragmentGraph *LikelyFragmentGraphGenerator::createNewGraph(){
	current_graph = new FragmentGraph();
	id_prob_computed_cache.clear();	//The graph is empty, so clear all computation records
	return current_graph;
}

//Create the starting node from a smiles or inchi string - responsibility of caller to delete
FragmentTreeNode *FragmentGraphGenerator::createStartNode( std::string &smiles_or_inchi ){
	
	//Create the RDKit mol - this will be the ion
	RDKit::RWMol *rwmol;
	if( smiles_or_inchi.substr(0,6) == "InChI=" ){
		RDKit::ExtraInchiReturnValues rv;
		rwmol = RDKit::InchiToMol( smiles_or_inchi, rv);
	}else
		rwmol = RDKit::SmilesToMol( smiles_or_inchi ); 

	//This is dirty, but for some reason RDKit doesn't throw the exception...
	if(!rwmol) throw RDKit::SmilesParseException("Error occurred - assuming Smiles Parse  Exception");

	//Compute and label the initial aromatic rings, Gasteiger Charges and Masses (with Hydrogens)
	RDKit::MolOps::removeStereochemistry( *rwmol );
	labelAromatics( rwmol );
	labelGasteigers( rwmol );
	initialiseRoots( rwmol );
	labelOriginalMasses( rwmol );

	//Kekulize and count the excess electron pairs
	RDKit::MolOps::Kekulize( *rwmol );
	int num_ep = countExtraElectronPairs( rwmol );

	return( new FragmentTreeNode(romol_ptr_t(rwmol), num_ep, 0) );
}

void FragmentGraphGenerator::initialiseRoots( RDKit::RWMol *rwmol ){
	RDKit::ROMol::AtomIterator ai;
	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai ){
		(*ai)->setProp("Root", 0 );
		(*ai)->setProp("OtherRoot", 0 );
	}
}

void FragmentGraphGenerator::labelGasteigers( RDKit::RWMol *rwmol ){
	// For each atom, will store the result in prop "OrigGasteigerCharge"
	RDKit::computeGasteigerCharges(rwmol); 
	RDKit::ROMol::AtomIterator ai;
	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai ){
		double gc;
		(*ai)->getProp<double>("_GasteigerCharge", gc);
		(*ai)->setProp("OrigGasteigerCharge", gc );
	}

}

void FragmentGraphGenerator::labelAromatics( RDKit::RWMol *rwmol ){	
	
	//Set bond aromaticity information
	RDKit::ROMol::BondIterator bi;
	for( bi = rwmol->beginBonds(); bi != rwmol->endBonds(); ++bi ){
		int aromatic = (*bi)->getIsAromatic();
		(*bi)->setProp("InAromaticRing", aromatic);
		(*bi)->setProp("InDblAromaticRing", 0);
	}

	//Check for any double-aromatic systems
	RDKit::MolOps::findSSSR( *rwmol );	
	RDKit::RingInfo *rinfo = rwmol->getRingInfo();
	std::vector<int> double_aromatic_idxs;
	for( unsigned int i = 0; i < rwmol->getNumBonds(); i++ ){
		if( rinfo->numBondRings(i) <= 1 ) continue;
		RDKit::Bond *bond = rwmol->getBondWithIdx(i);
		if( bond->getIsAromatic() ) 
			double_aromatic_idxs.push_back(i);
	}
	
	//If any are found, label all the bonds within them
	if(double_aromatic_idxs.size() == 0 ) return;
		
	//Consider each ring...
	RDKit::RingInfo::VECT_INT_VECT brings = rinfo->bondRings();
	RDKit::RingInfo::VECT_INT_VECT::iterator bit = brings.begin();
	for( ; bit != brings.end(); ++bit ){

		//Check for a double aromatic bond within the ring
		bool hasDblArom = false;
		RDKit::RingInfo::INT_VECT::iterator it;
		for( it = bit->begin(); it != bit->end(); ++it ){
			std::vector<int>::iterator ii = double_aromatic_idxs.begin();
			for( ; ii != double_aromatic_idxs.end(); ++ii )
				if( *ii == *it ) hasDblArom = true;
			if(hasDblArom) break;
		}

		//If one exists, label all bonds in the ring
		if( !hasDblArom ) continue;
		for( it = bit->begin(); it != bit->end(); ++it ){
			RDKit::Bond *bond = rwmol->getBondWithIdx(*it);
			bond->setProp("InDblAromaticRing", 1);
		}
		
	}
}

void FragmentGraphGenerator::labelOriginalMasses( RDKit::RWMol *rwmol ){
	RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
	RDKit::ROMol::AtomIterator ai;
	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai ){
		double mass = 0.0;
		std::string symbol = (*ai)->getSymbol();
		mass += pt->getMostCommonIsotopeMass(symbol);
		mass += (*ai)->getTotalNumHs()*pt->getMostCommonIsotopeMass("H");
		(*ai)->setProp("OriginalMass", mass);
	}
}

int FragmentGraphGenerator::countExtraElectronPairs( RDKit::RWMol *rwmol ){

	//Compute the total number of bond electrons
	double total_bond_es = 0;
	RDKit::ROMol::AtomIterator ai;
	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai ){
	
		int Hs = (*ai)->getTotalNumHs();
		double charge = (*ai)->getFormalCharge();
		if( charge != 0 ){ 
			std::cout << "Warning: expecting uncharged atoms. This probably won't work..." << std::endl;
			return 0;
		}
		int implied_val = (*ai)->getExplicitValence() + (*ai)->getImplicitValence();
		total_bond_es += (implied_val - Hs);

	}

	//Adjust to figure out how many are extras (non-single bond electrons)
	//and divide by 2 to count pairs
	double extra_bond_eps = (total_bond_es - rwmol->getNumBonds()*2)/2;
	return (int)extra_bond_eps;

}

//Helper function - check if the fragment has already been computed to at least this depth
int FragmentGraphGenerator::alreadyComputed(int id, int remaining_depth){
	if( id_depth_computed_cache.find(id) == id_depth_computed_cache.end()	//Not found
		|| id_depth_computed_cache[id] < remaining_depth){					//Or computed previously to less depth
		id_depth_computed_cache[id] = remaining_depth;
		return 0;
	}
	return 1;
}

//Compute a FragmentGraph starting at the given node and computing to the depth given.
//The output will be appended to the current_graph
void FragmentGraphGenerator::compute( FragmentTreeNode &node, int remaining_depth, int parentid ){
	
	//Add the node to the graph, and return a fragment id
	int id = -1;
	if( mols_to_fv ) 
		id = current_graph->addToGraphAndReplaceMolWithFV( node.ion, node.nl, parentid, fc );
	else 
		id = current_graph->addToGraph( node.ion, node.nl, parentid );

	//Only compute to the desired depth
	if( remaining_depth <= 0 ) return;

	//If the node was already in the graph at sufficient depth, skip any further computation
	if( alreadyComputed(id, remaining_depth) ){ 
		if(verbose)	std::cout << "Node already computed: Skipping" << std::endl;
		return;
	}

	//Generate Breaks
	std::vector<Break> breaks;
	node.generateBreaks( breaks );

	//Iterate over the possible breaks
	std::vector<Break>::iterator it = breaks.begin();
	for( ; it != breaks.end(); ++it ){

		node.applyBreak( *it );
		node.generateChildrenOfBreak( *it );

		//Recur over children
		std::vector<FragmentTreeNode>::iterator itt = node.children.begin();
		for( ; itt != node.children.end(); ++itt ){
			compute( *itt, remaining_depth-1, id );
		}

		//Undo and remove children
		node.undoBreak( *it );
		node.children = std::vector<FragmentTreeNode>();
	}
}


//Compute a FragmentGraph starting at the given node and computing to the depth given.
//The output will be appended to the current_graph
void LikelyFragmentGraphGenerator::compute( FragmentTreeNode &node,  int remaining_depth, int parentid, double parent_log_prob ){
	
	//Add the node to the graph, and return a fragment id: note, no mols or fv will be set,
	//but the precomputed theta value will be used instead
	int id = -1;
	id = current_graph->addToGraphWithThetas( node.ion, node.nl, node.getAllTmpThetas(), parentid );

	//Reached max depth? 
	if( remaining_depth <= 0 ) return;

	//If we've already run the fragmentation on this fragment with an equal or higher
	//probability offset, we don't need to run again, unless it is persisting
	if( alreadyComputed( id, parent_log_prob )) return;

	//Generate Children 
	std::vector<Break> breaks;
	node.generateBreaks( breaks );
	std::vector<Break>::iterator it = breaks.begin();
	for( ; it != breaks.end(); ++it ){
		node.applyBreak( *it );
		node.generateChildrenOfBreak( *it );
		node.undoBreak( *it );
	} 

	//Compute child thetas
	std::vector<FragmentTreeNode>::iterator itt = node.children.begin();
	for( ; itt != node.children.end(); ++itt ){ 
		Transition tmp_t( -1, -1, itt->nl, itt->ion ); 
		FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
		for( int engy = param->getNumEnergyLevels()-1; engy >= 0; engy-- )
			itt->setTmpTheta( param->computeTheta(*fv, engy), engy );
		delete fv;
	}

	//Compute child probabilities (including persistence) - for all energy levels
	std::vector<double> denom( param->getNumEnergyLevels() );
	for(int i = 0; i < denom.size(); i++ ) denom[i] = 0.0;
	for( itt = node.children.begin(); itt != node.children.end(); ++itt ){
		for(int energy = 0; energy < denom.size(); energy++ )
			denom[energy] = logAdd( denom[energy], itt->getTmpTheta(energy) );
	}

	//Add and recur over likely children if above threshold for any energy level
	double max_child_prob;
	for( itt = node.children.begin(); itt != node.children.end(); ++itt ){
		max_child_prob = log_prob_thresh - 10.0;
		for(int energy = 0; energy < denom.size(); energy++ ){		
			double child_log_prob = itt->getTmpTheta(energy) - denom[energy] + parent_log_prob;
			if( child_log_prob > max_child_prob ) max_child_prob = child_log_prob;
		}
		if( max_child_prob >= log_prob_thresh )
			compute( *itt, remaining_depth-1, id, max_child_prob );
	}

	//Clear the children
	node.children = std::vector<FragmentTreeNode>();
}


//Helper function - check if the fragment has already been computed with at least this probability offset
int LikelyFragmentGraphGenerator::alreadyComputed(int id, double prob_offset){
	if( id_prob_computed_cache.find(id) == id_prob_computed_cache.end()
		|| id_prob_computed_cache[id] < prob_offset){					//Or computed previously with lower prob
		id_prob_computed_cache[id] = prob_offset;
		return 0;
	}
	return 1;
}
