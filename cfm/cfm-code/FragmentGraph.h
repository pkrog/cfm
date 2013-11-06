/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentGraph.h
#
# Description: 	FragmentGraph class for holding the results of a generated
#				fragment graph.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __FRAGTREE_H__
#define __FRAGTREE_H__

#include <GraphMol/ROMol.h>
#include <vector>

#include "Util.h"
#include "Features.h"

typedef std::vector<std::vector<int> > tmap_t;

//Class for storing the base fragment state for our model
class Fragment{

public:
	//Constructor, store the ion smiles and a reduced smiles since the ion is
	//not needed and takes more space.
	Fragment( std::string &a_ion_smiles, std::string &a_reduced_smiles, int an_id, double a_mass ) : 
		id( an_id ), ion_smiles(a_ion_smiles), reduced_smiles( a_reduced_smiles ), mass( a_mass ){};

	//Access Functions
	double getMass() const { return mass; };
	int getId() const { return id; };
	void setId( int an_id ){ id = an_id;};
	const std::string *getReducedSmiles() const {return &reduced_smiles; };
	const std::string *getIonSmiles() const {return &ion_smiles; };
	void clearSmiles(){reduced_smiles = std::string(); ion_smiles = std::string();};

private:
	int id;
	std::string reduced_smiles;	//Reduced version of the smiles string (just backbone)
	std::string ion_smiles;		//Full ion smiles (for writing out if called for)
	double mass;
};

//Class for storing possible transitions between fragments
class Transition{

public:
	//Basic constructor
	Transition( int a_from_id, int a_to_id, RootedROMolPtr &a_nl, RootedROMolPtr &an_ion  ) :
		from_id( a_from_id), to_id( a_to_id ), nl( a_nl ), ion( an_ion ) {}; 

	//Alternative constructor that finds the root atoms and sets the
	//root pointers appropriately
	Transition( int a_from_id, int a_to_id, romol_ptr_t &a_nl, romol_ptr_t &an_ion  );
	
	//Access Functions	
	int getFromId() const { return from_id; };
	int getToId() const { return to_id; };
	const RootedROMolPtr *getNeutralLoss() const { return &nl; };
	void deleteNeutralLoss(){ nl.mol.reset(); nl = RootedROMolPtr(); };
	const RootedROMolPtr *getIon() const { return &ion; };
	void deleteIon(){ ion.mol.reset(); ion = RootedROMolPtr(); };
	FeatureVector *getTmpFV() const { return tmp_fv; };
	void setTmpFV( FeatureVector *an_fv_ptr ){ tmp_fv = an_fv_ptr; };
	const std::vector<double> *getTmpThetas() const { return &tmp_thetas; };
	void setTmpThetas( const std::vector<double> *a_thetas ){tmp_thetas = *a_thetas;}; 

private:
	int from_id;
	int to_id;
	RootedROMolPtr nl;
	RootedROMolPtr ion;	//We store the ion on the transition to
						//allow for different roots - the fragment stores
						//only an unrooted shared pointer.
	FeatureVector *tmp_fv;	//Temporary storage for the feature vector pointer
							//while we compute the fragment graph (don't use this
							//directly - it will be moved up into the MolData)
	std::vector<double> tmp_thetas;	
							//Temporary storage for the theta values
							//while we compute the likely fragment graph
							//(as above, don't use directly, will be moved)
};

class FragmentGraph{
public:
	FragmentGraph(){};

	//Add a fragment node to the graph (should be the only way to modify the graph)
	//	-- Add a fragment, or return an id, if it already exists
	//  -- Add a transition to this fragment, based on the provided parent fragment id
	//  -- Update the relevant tmaps
	// Note: If parentid < 0, assumes starting ion, so adds extra H to mass, and
	//       doesn't add transition.
	int addToGraph( romol_ptr_t ion, romol_ptr_t nl, int parentid );
	
	//As for previous function, but delete the mols in the transition and compute and store a feature vector instead
	int addToGraphAndReplaceMolWithFV( romol_ptr_t ion, romol_ptr_t nl, int parentid, FeatureCalculator *fc );

	//As for previous function, but don't store the mols in the transition and insert the pre-computed thetas instead
	int addToGraphWithThetas(romol_ptr_t ion, romol_ptr_t nl, const std::vector<double> *thetas, int parentid );

	//Write the Fragments only to file (formerly the backtrack output - without extra details)
	void writeFragmentsOnly( std::ostream &out );

	//Write the FragmentGraph to file (formerly the transition output - without feature details)
	void writeFullGraph( std::ostream &out );

	//Access functions
	unsigned int getNumTransitions() const { return transitions.size(); };
	unsigned int getNumFragments() const { return fragments.size(); };
	const Transition *getTransitionAtIdx( int index ) const { return &(transitions[index]); };
	const Fragment *getFragmentAtIdx( int index ) const { return &(fragments[index]); };
	const tmap_t *getFromIdTMap() const { return &from_id_tmap; };
	const tmap_t *getToIdTMap() const { return &to_id_tmap; };

	//Functions to delete shared romol pointers for ion and nl mols
	void deleteMolsForTransitionAtIdx( int index );
	void clearAllSmiles();

private:
	std::vector<Fragment> fragments;
	std::vector<Transition> transitions;
	tmap_t from_id_tmap;	//Mapping between from_id and transitions with that from_id
	tmap_t to_id_tmap;		//Mapping between to_id and transitions with that to_id

	//Mapping from rounded mass to list of fragment ids, 
	//to enable fast check for existing fragments
	std::map<double, std::vector<int> > frag_mass_lookup;  

	//Find the id for an existing fragment that matches the input ion and mass
	//or create a new fragment in the case where no such fragment is found
	int addFragmentOrFetchExistingId( romol_ptr_t ion, double mass );

	//Determine if the two fragments match - assumes the masses have already
	//been checked to be roughly the same, now check reduced structure.
	bool areMatching( RDKit::ROMol *f1_reduced_ion, RDKit::ROMol *f2_reduced_ion );

	//Make all bonds single and remove any charge - used to compare structure alone
	void reduceMol( RDKit::RWMol &rwmol );

	//Find the id for an existing transition that matches the input ids
	//or -1 in the case where no such transition is found
	int findMatchingTransition( int from_id, int to_id );
};


#endif // __FRAGTREE_H__
