/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentTreeNode.cpp
#
# Description: 	Contains Break and FragmentTreeNode classes, for use in
#				the fragment generation process.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "FragmentTreeNode.h"
#include "MILP.h"

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/RingInfo.h>

//Constructor for a Non-Ring Break
Break::Break( int a_bond_idx ){
	bond_idxs[0] = a_bond_idx;
	ring_break = 0;
}

//Constructor for a Ring Break
Break::Break( int a_bond1_idx, int a_bond2_idx, int a_ring_idx ){
	bond_idxs[0] = a_bond1_idx;
	bond_idxs[1] = a_bond2_idx;
	ring_break = 1;
	ring_idx = a_ring_idx;
}

void FragmentTreeNode::generateChildrenOfBreak( Break &brk ){
	
	int f0_max_e, f1_max_e = 0;
	std::vector<int> f0_output_bmax, f1_output_bmax;
	
	int isringbrk = brk.isRingBreak();
	int brk_ringidx = isringbrk*brk.getRingIdx() - (1-isringbrk);

	//Compute the max electron assignment for F0
	MILP f0_solver(ion.get(), 0, brk_ringidx);
	f0_max_e = f0_solver.runSolver(f0_output_bmax);

	//Compute the max electron assignment for F1
	MILP f1_solver(ion.get(), 1, brk_ringidx);
	f1_max_e = f1_solver.runSolver(f1_output_bmax);

	//Combine the electron allocations of the two fragments
	for( unsigned int i = 0; i < f0_output_bmax.size(); i++ )
		f0_output_bmax[i] += f1_output_bmax[i];

	//Determine the limits for how many electron pairs can be allocated
	int total_free_epairs = ion_free_epairs + 1 + brk.isRingBreak();
	int min_f0 = std::max(total_free_epairs - f1_max_e, 0 );
	int max_f0 = std::min(total_free_epairs, f0_max_e );

	//Iterate through the possible solutions, adding them to the node as children
	for( int e_f0 = min_f0; e_f0 <= max_f0; e_f0++ ){
		addBothChildren( e_f0, total_free_epairs, f0_output_bmax, brk);	//Child with each side as ion
	}
}

void FragmentTreeNode::addBothChildren( int e_f0, int e_to_allocate, std::vector<int> &output_bmax, Break &brk ){

	RDKit::RWMol rwmol = *ion;	//Copy the ion
	//Set the correct bond orders on the fragments
	int broken, fragidx;
	int numbonds = rwmol.getNumBonds();
	int remaining_e[2] = {e_f0, e_to_allocate - e_f0};
	int allocated_e[2] = {e_f0, e_to_allocate - e_f0};
	std::vector<std::pair<int,int> > broken_specs;
	for( int i = 0; i < numbonds; i++ ){
		RDKit::Bond *bond = rwmol.getBondWithIdx(i);
		bond->getProp("Broken", broken);
		if( broken ){ 
			//Save the bond to remove at the end, removing it now causes problems with the iteration/indexing
			broken_specs.push_back( std::pair<int,int>(bond->getBeginAtomIdx(), bond->getEndAtomIdx()) );
			continue;
		}
		bond->setIsAromatic(false);
		bond->getBeginAtom()->getProp("FragIdx", fragidx);
		if( remaining_e[fragidx] > 0){
			int num_to_add = std::min( output_bmax[i], remaining_e[fragidx] );
			bond->setBondType( RDKit::Bond::BondType(1 + num_to_add ) );
			remaining_e[fragidx] -= num_to_add;
		}
		else bond->setBondType( RDKit::Bond::SINGLE );
	}

	//Remove the broken bonds
	std::vector<std::pair<int,int> >::iterator itr = broken_specs.begin();
	for( ; itr != broken_specs.end(); ++itr )
		rwmol.removeBond(itr->first, itr->second);

	//Fill in the Hydrogens
	int orig_val;
	for( unsigned int i = 0; i < rwmol.getNumAtoms(); i++ ){
		RDKit::Atom *atom = rwmol.getAtomWithIdx(i);
		atom->setNoImplicit(true);
		atom->setNumExplicitHs( 0 );
		atom->setFormalCharge( 0 );
		atom->setIsAromatic(false);
		int val = atom->calcExplicitValence();	//Valence excluding hydrogens
		atom->getProp( "OrigValence", orig_val );
		atom->setNumExplicitHs( orig_val - val );
		atom->calcExplicitValence();
	}

	for( int charge_frag = 0; charge_frag <= 1; charge_frag++ ){

		//Set the charge location: (non-ring) N, S, O, then C+
		int charge_idx = findChargeLocation( rwmol, charge_frag );
		if( charge_idx < 0 ) continue;
		RDKit::Atom *atom = rwmol.getAtomWithIdx(charge_idx);
		atom->setFormalCharge(1);
		atom->setNumExplicitHs( atom->getTotalNumHs() + 1 );
		atom->calcExplicitValence();

		//Separate the ion and nl mols and add the child node
		std::vector< boost::shared_ptr< RDKit::ROMol > > mols = RDKit::MolOps::getMolFrags( rwmol );
		if( mols.size() == 2 ){
			int ion_idx = 1 - RDKit::MolOps::getFormalCharge( *mols[0].get() );
		
			//Record some of the properties of the break in the neutral loss
			//e.g. ring break? aromaticity etc.
			romol_ptr_t nl = mols[1-ion_idx];
			nl.get()->setProp("IsRingBreak",brk.isRingBreak());
			RDKit::Bond *bond = ion->getBondWithIdx( brk.getBondIdx() );
			int is_arom, is_dbl_arom;
			bond->getProp("InAromaticRing", is_arom);
			bond->getProp("InDblAromaticRing", is_dbl_arom);
			nl.get()->setProp("IsAromaticRingBreak", is_arom);
			nl.get()->setProp("IsAromaticDblRingBreak", is_dbl_arom);
		
			//Add the child node
			children.push_back( FragmentTreeNode( mols[ion_idx], nl, allocated_e[charge_frag], depth+1));
		}
		//Undo the charge
		atom->setFormalCharge(0);
		atom->setNumExplicitHs( atom->getTotalNumHs() - 1 );
		atom->calcExplicitValence();
	}
}
	
	
int FragmentTreeNode::findChargeLocation( RDKit::RWMol &rwmol, int charge_frag ){
	
	int charge_loc = -1;
	int fragidx;
	unsigned int numrings;
	int Nidx = -1, Sidx = -1, Oidx = -1, Cidx = -1, anyCidx = -1;
	for( unsigned int i = 0; i < rwmol.getNumAtoms(); i++ ){
		RDKit::Atom *atom = rwmol.getAtomWithIdx(i);
		atom->getProp("FragIdx", fragidx);
		if( fragidx != charge_frag ) continue;
		if( atom->getSymbol() == "C" ) anyCidx = i;
		atom->getProp("NumUnbrokenRings", numrings);		
		if( numrings > 0 ) continue;
		if( atom->getSymbol() == "N" && Nidx < 0 ) Nidx = i;
		else if( atom->getSymbol() == "S" && Sidx < 0 ) Sidx = i;
		else if( atom->getSymbol() == "O" && Oidx < 0 ) Oidx = i;
		else if( atom->getSymbol() == "C" && Cidx < 0 ) Cidx = i;
	}
	if( Nidx >= 0 ) charge_loc = Nidx;		 //Non-Ring Nitrogen
	else if( Sidx >= 0 ) charge_loc = Sidx;  //Non-Ring Sulfur
	else if( Oidx >= 0 ) charge_loc = Oidx;  //Non-Ring Oxygen
	else if( Cidx >= 0 ) charge_loc = Cidx;  //Non-Ring Carbon
	else if( anyCidx >= 0 ) charge_loc = anyCidx; //Any Carbon
	
	return charge_loc;

}

void FragmentTreeNode::generateBreaks(std::vector<Break> &breaks){

	//Populate the Ring Info for the ion and set the NumUnbrokenRings and OrigValence properties
	RDKit::MolOps::findSSSR( *ion.get() );	
	RDKit::RingInfo *rinfo = ion.get()->getRingInfo();
	RDKit::ROMol::AtomIterator ai; 
	for( ai = ion.get()->beginAtoms(); ai != ion.get()->endAtoms(); ++ai ){
		(*ai)->setProp("NumUnbrokenRings", rinfo->numAtomRings((*ai)->getIdx()) );
		(*ai)->setProp("OrigValence", (*ai)->getExplicitValence() + (*ai)->getImplicitValence() - (*ai)->getFormalCharge() );
		(*ai)->setProp("Root", 0);
		(*ai)->setProp("OtherRoot", 0);
	}

	//Generate Non-Ring Breaks
	for( unsigned int bidx = 0; bidx < ion.get()->getNumBonds(); bidx++ ){
		RDKit::Bond *bond = ion.get()->getBondWithIdx(bidx);
		bond->setProp("Broken",0);
		bond->setProp("NumUnbrokenRings", rinfo->numBondRings(bidx));
		if( rinfo->numBondRings(bidx) == 0 )
			breaks.push_back( Break(bidx) );
	}
	
	//Ring Breaks
	RDKit::RingInfo::VECT_INT_VECT brings = rinfo->bondRings();
	RDKit::RingInfo::VECT_INT_VECT::iterator bit = brings.begin();
	for( int ringidx = 0; bit != brings.end(); ++bit, ringidx++ ){

		//All pairs of bonds within the ring, that don't
		//belong to any other ring
		RDKit::RingInfo::INT_VECT::iterator it1, it2;
		for( it1 = bit->begin(); it1 != bit->end(); ++it1 ){
			if( rinfo->numBondRings(*it1) != 1 ) continue;
			for( it2 = it1 + 1; it2 != bit->end(); ++it2 ){
				if( rinfo->numBondRings(*it2) != 1 ) continue;
				breaks.push_back( Break( *it1, *it2, ringidx) );
			}
		}
	}
}

void FragmentTreeNode::applyBreak(Break &brk){

	//Label the broken bond and root atoms
	RDKit::Bond *broken_bond = ion.get()->getBondWithIdx( brk.getBondIdx() );
	broken_bond->setProp( "Broken", 1 );
	broken_bond->getBeginAtom()->setProp("Root",1);
	broken_bond->getEndAtom()->setProp("Root",1);
	if( brk.isRingBreak() ){
		RDKit::Bond *broken_bond2 = ion.get()->getBondWithIdx( brk.getSecondBondIdx() );
		broken_bond2->setProp( "Broken", 1 );
		broken_bond2->getBeginAtom()->setProp("OtherRoot",1);
		broken_bond2->getEndAtom()->setProp("OtherRoot",1);
	}

	//Label fragment idxs
	for( unsigned int i = 0; i < ion.get()->getNumAtoms(); i++ ){
		RDKit::Atom *atom = ion.get()->getAtomWithIdx(i);
		atom->setProp( "FragIdx", 0 );
	}
	allocatedCtdToFragment( ion.get(), broken_bond->getBeginAtom() );

	//Decrement the NumUnbrokenRings setting for atoms and bonds in a broken ring
	if( brk.isRingBreak() ){
		int ringidx = brk.getRingIdx();
		RDKit::RingInfo *rinfo = ion.get()->getRingInfo();
		RDKit::RingInfo::INT_VECT::iterator it;		

		RDKit::RingInfo::VECT_INT_VECT arings = rinfo->atomRings();
		for( it = arings[ringidx].begin(); it != arings[ringidx].end(); ++it ){
			unsigned int numrings;
			RDKit::Atom *atom = ion.get()->getAtomWithIdx(*it);
			atom->getProp("NumUnbrokenRings", numrings);
			atom->setProp("NumUnbrokenRings", numrings-1);
		}

		RDKit::RingInfo::VECT_INT_VECT brings = rinfo->bondRings();
		for( it = brings[ringidx].begin(); it != brings[ringidx].end(); ++it ){
			unsigned int numrings;
			RDKit::Bond *bond = ion.get()->getBondWithIdx(*it);
			bond->getProp("NumUnbrokenRings", numrings);
			bond->setProp("NumUnbrokenRings", numrings-1);
		}
	}
}

void FragmentTreeNode::undoBreak(Break &brk){

	//Label the broken bond and root atoms
	RDKit::Bond *broken_bond = ion.get()->getBondWithIdx( brk.getBondIdx() );
	broken_bond->setProp( "Broken", 0 );
	broken_bond->getBeginAtom()->setProp("Root",0);
	broken_bond->getEndAtom()->setProp("Root",0);

	if( brk.isRingBreak() ){
		RDKit::Bond *broken_bond2 = ion.get()->getBondWithIdx( brk.getSecondBondIdx() );
		broken_bond2->setProp( "Broken", 0 );
		broken_bond2->getBeginAtom()->setProp("OtherRoot",0);
		broken_bond2->getEndAtom()->setProp("OtherRoot",0);
	}

	//Increment the NumUnbrokenRings setting for atoms and bonds in a broken ring
	if( brk.isRingBreak() ){
		int ringidx = brk.getRingIdx();
		RDKit::RingInfo *rinfo = ion.get()->getRingInfo();
		RDKit::RingInfo::INT_VECT::iterator it;		

		RDKit::RingInfo::VECT_INT_VECT arings = rinfo->atomRings();
		for( it = arings[ringidx].begin(); it != arings[ringidx].end(); ++it ){
			unsigned int numrings;
			RDKit::Atom *atom = ion.get()->getAtomWithIdx(*it);
			atom->getProp("NumUnbrokenRings", numrings);
			atom->setProp("NumUnbrokenRings", numrings+1);
		}

		RDKit::RingInfo::VECT_INT_VECT brings = rinfo->bondRings();
		for( it = brings[ringidx].begin(); it != brings[ringidx].end(); ++it ){
			unsigned int numrings;
			RDKit::Bond *bond = ion.get()->getBondWithIdx(*it);
			bond->getProp("NumUnbrokenRings", numrings);
			bond->setProp("NumUnbrokenRings", numrings+1);
		}
	}
}

void FragmentTreeNode::allocatedCtdToFragment( RDKit::ROMol *romol, RDKit::Atom *atom ){
	
	int broken, fragidx;
	atom->setProp( "FragIdx", 1 );
	RDKit::ROMol::ADJ_ITER_PAIR itp = romol->getAtomNeighbors( atom );
	for( ; itp.first != itp.second; ++itp.first ){
		RDKit::Atom *nbr_atom = romol->getAtomWithIdx(*itp.first);
		nbr_atom->getProp("FragIdx",fragidx);
		RDKit::Bond *bond = romol->getBondBetweenAtoms( atom->getIdx(), nbr_atom->getIdx() );
		bond->getProp("Broken", broken);
		if(!fragidx && !broken) allocatedCtdToFragment( romol, nbr_atom );
	}
}
