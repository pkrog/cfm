/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# milp.c
#
# Description: 	Functions for running the milp solver looking for
#				valid assignment of bonds and hydrogens.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "MILP.h"
#include "lp_lib.h"
#include <GraphMol/RingInfo.h>
#include <GraphMol/ROMol.h>


MILP::MILP( RDKit::ROMol *a_mol, int a_fragmentidx, int a_broken_ringidx  ){
	mol = a_mol;
	fragmentidx = a_fragmentidx;
	broken_ringidx = a_broken_ringidx;
}

MILP::MILP( RDKit::ROMol *a_mol, int a_fragmentidx ){
	mol = a_mol;
	fragmentidx = a_fragmentidx;
	broken_ringidx = -1;
}

int MILP::runSolver( std::vector<int> &output_bmax ){
   //Uses lp_solve: based on demonstration code provided.
  unsigned int numunbroken;
  int broken, fragidx, origval;
  int Ncol, *colno = NULL, i, ret = 0, output_max_e=0;
  REAL *row = NULL;
  lprec *lp;

  // Variables are: num additional electron pairs
  // added per bond-> i.e. 0 = single, 1 = double, 2 = triple
  int num_bonds = mol->getNumBonds();
  int num_atoms = mol->getNumAtoms();
  Ncol = num_bonds;
  lp = make_lp(0, num_bonds);
  if(lp == NULL)
    ret = 1; /* couldn't construct a new model... */

  if(ret == 0) {
    colno = (int *) malloc(num_bonds * sizeof(*colno));
    row = (REAL *) malloc(num_bonds * sizeof(*row));
    if((colno == NULL) || (row == NULL))
      ret = 2;
  }

  //Bond Constraints
  if(ret == 0) {
    set_add_rowmode(lp, TRUE);

	//All bonds must be at most TRIPLE bonds ( <= 2 ) 
	//except broken bonds ( <= 0 ) and ring bonds ( <= 1 )
	for( i = 0; i < num_bonds && ret == 0; i++ ){
		set_int(lp, i+1, TRUE); //sets variable to integer

		RDKit::Bond *bond = mol->getBondWithIdx(i);
		int limit = 0;	  //bonds that are broken or in the other fragment are limited to 0
		
		bond->getProp("Broken", broken);
		bond->getBeginAtom()->getProp("FragIdx", fragidx);
		if( !broken && fragidx == fragmentidx ){
			limit = 2;
			bond->getProp("NumUnbrokenRings", numunbroken);
			if( numunbroken > 0 ) limit = 1;
		}
		colno[0] = i + 1; //variable idx i.e. bond
		row[0] = 1;		  //multiplier

		if(!add_constraintex(lp, 1, row, colno, LE, limit))
		  ret = 3;
	}

	//Add atom valence constraints to neighbouring bonds
	for( i = 0; i < num_atoms && ret == 0; i++ ){
		RDKit::Atom *atom = mol->getAtomWithIdx(i);
		atom->getProp("FragIdx", fragidx);
		atom->getProp("OrigValence", origval);
		if( fragidx != fragmentidx ) continue;
		

		RDKit::ROMol::OBOND_ITER_PAIR ip = mol->getAtomBonds( atom );
		RDKit::ROMol::OEDGE_ITER it = ip.first;
		int num_unbroken = 0, j = 0;
		for( ; it != ip.second; ++it ){
			RDKit::Bond *cbond =(*mol)[*it].get();
			cbond->getProp("Broken", broken);
			colno[j] = cbond->getIdx() + 1; //variable idx i.e. bond
			row[j++] = 1;		  //multiplier
			if( !broken ) num_unbroken++;
		}
		if(!add_constraintex(lp, j, row, colno, LE, origval - num_unbroken))
		  ret = 3;
	}

	//Add constraints for neighbouring ring bonds (can't have two double in a row)
	RDKit::RingInfo *rinfo = mol->getRingInfo();
	RDKit::RingInfo::VECT_INT_VECT brings = rinfo->bondRings();
	RDKit::RingInfo::VECT_INT_VECT::iterator bit = brings.begin();
	for( int ringidx = 0; bit != brings.end() && ret == 0; ++bit, ringidx++ ){

		if( ringidx == broken_ringidx ) continue;
		row[0] = 1; row[1] = 1;
		
		//Create a vector of flags indicating bonds included in the ring
		std::vector<int> ring_bond_flags(num_bonds);
		for( int i = 0; i < num_bonds;i++ ) ring_bond_flags[i] = 0;
		RDKit::RingInfo::INT_VECT::iterator it;
		for( it = bit->begin(); it != bit->end(); ++it ) ring_bond_flags[*it] = 1;
		
		//Traverse around the ring, creating the constraints
		RDKit::Bond *bond = mol->getBondWithIdx(*(bit->begin())); //Starting Bond	
		RDKit::Bond *start_bond = bond, *prev_bond = bond;
		RDKit::Atom *atom = bond->getBeginAtom();
		while( ret == 0 && (bond = getNextBondInRing(bond, atom, ring_bond_flags)) != start_bond ){
			colno[0] = prev_bond->getIdx() + 1;
			colno[1] = bond->getIdx() + 1;
			if(!add_constraintex(lp, 2, row, colno, LE, 1)) ret = 3;
			atom = bond->getOtherAtom( atom );
			prev_bond = bond;
		}
	}
  }

  if(ret == 0) {
	// set the objective function = sum bond electrons
	set_add_rowmode(lp, FALSE);
	for( int j = 0; j < num_bonds; j++ ){
		colno[j] = j+1; 
		row[j] = 1;	
	}
	if(!set_obj_fnex(lp, num_bonds, row, colno))
		ret = 4;
  }

  //Run optimization
  if(ret == 0) {
    set_maxim(lp);
    set_verbose(lp, IMPORTANT);
    ret = solve(lp);
	if(ret == OPTIMAL) ret = 0;
    else ret = 5;
  }

  //Extract Results
  if(ret == 0) {
	output_bmax.resize( Ncol );
	get_variables(lp, row);
	for( int j = 0; j < Ncol; j++) output_bmax[j] = (int)row[j];
	output_max_e = (int)get_objective(lp);
  }

  // Free allocated memory
  if(row != NULL) free(row);
  if(colno != NULL) free(colno);
  if(lp != NULL) delete_lp(lp);
  
  //Return maximum non-single bond electron pairs that can 
  //be allocated to this fragment.
  return output_max_e;
}

//Helper function - allows traversal of a ring one bond at a time
RDKit::Bond *MILP::getNextBondInRing( RDKit::Bond *bond, RDKit::Atom *atom, std::vector<int> &ring_bond_flags){

	RDKit::ROMol::OBOND_ITER_PAIR ip;
	ip = atom->getOwningMol().getAtomBonds( atom );
	RDKit::ROMol::OEDGE_ITER it = ip.first;
	for( ; it != ip.second; ++it ){
		int idx = (*mol)[*it].get()->getIdx();	//There must be a better way...
		if( ring_bond_flags[idx] && idx != bond->getIdx() )
			return atom->getOwningMol().getBondWithIdx(idx);
	}
	return NULL;	//Error@!
}