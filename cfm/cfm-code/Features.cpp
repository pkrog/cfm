/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features.cpp
#
# Description: 	Code for computing features for fragmentations.
#
#				Assume that we have a config file that lists the feature
#				vectors to compute (line separated text).
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>
#include <cmath>
#include <queue>

#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/PeriodicTable.h>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include "Features.h"

typedef std::pair<std::string,std::string> symbol_pair_t;

const boost::ptr_vector<Feature> &FeatureCalculator::featureCogs(){
	
	static boost::ptr_vector<Feature> cogs;
	static bool initialised = false;

	if (!initialised) {
		cogs.push_back( new BreakAtomPair() );
		cogs.push_back( new GasteigerCharges() );
		cogs.push_back( new HydrogenMovement() );
		cogs.push_back( new IonRootPairs() );
		cogs.push_back( new IonRootTriples() );
		cogs.push_back( new NLRootPairs() );
		cogs.push_back( new NLRootTriples() );
		cogs.push_back( new RingFeatures() );
		cogs.push_back( new QuadraticFeatures() );
		initialised = true;
	}
	return cogs;
}

FeatureCalculator::FeatureCalculator( std::string &config_filename ){
		
	//Read the config file into a set containing the names of the features 
	std::ifstream ifs ( config_filename.c_str() , std::ifstream::in );
	
	if(!ifs.good()){
		std::cout << "Trouble opening feature config file: " << config_filename << std::endl;
	}
	while( ifs.good() ){ 
			
		std::string name;
		getline( ifs, name );
		boost::trim(name);
		if( name.length() <= 1 ) break;
		configureFeature( name );		
		
	}

	if( used_feature_idxs.size() == 0 ){
		std::cout << "Error reading config file, no features found" << std::endl;
		throw( InvalidConfigException() );
	}
}

FeatureCalculator::FeatureCalculator( std::vector<std::string> &feature_list ){

	//Find the relevant feature cog for this name
	std::vector<std::string>::iterator itname = feature_list.begin();
	for( ; itname != feature_list.end(); ++itname )
		configureFeature( *itname );	

	if( used_feature_idxs.size() == 0 )
		std::cout << "Warning: No features found in feature list" << std::endl;

}

std::vector<std::string> FeatureCalculator::getFeatureNames(){

	std::vector<std::string> names;
	std::vector<int>::iterator it = used_feature_idxs.begin();
	for( ; it != used_feature_idxs.end(); ++it ){
		const Feature *cog = &(featureCogs()[*it]);
		names.push_back( cog->getName() );
	}
	return names;

}

const std::vector<std::string> FeatureCalculator::getValidFeatureNames(){
	
	static std::vector<std::string> output;
	static bool initialised = false;
	if( initialised ) return output;	

	boost::ptr_vector<Feature>::const_iterator it = featureCogs().begin();
	for( ; it != featureCogs().end(); ++it )
		output.push_back( it->getName() );
	
	initialised = true;
	return output;
}

void FeatureCalculator::configureFeature( std::string &name ){

	//Find the relevant feature cog for this name
	boost::ptr_vector<Feature>::const_iterator it = featureCogs().begin();
	bool found = false;
	for( int idx = 0; it != featureCogs().end(); ++it, idx++ ){
		if( it->getName() == name ){ 
			used_feature_idxs.push_back(idx);
			found = true;
			break;
		}
	}
	if( !found ){
		std::cout << "Unrecognised feature: " << name << std::endl;
		throw( InvalidConfigException() );
	}

}

unsigned int FeatureCalculator::getNumFeatures(){

	unsigned int count = 1;	//Bias
	int quadratic = 0;
	std::vector<int>::iterator it = used_feature_idxs.begin();
	for( ; it != used_feature_idxs.end(); ++it ){
		count += featureCogs()[*it].getSize();
		if(  featureCogs()[*it].getName() == "QuadraticFeatures" )
			quadratic = 1;
	}
	if(quadratic) count += (count-1)*(count-2)/2;
	return count;
}


FeatureVector *FeatureCalculator::computeFV( const RootedROMolPtr *ion, const RootedROMolPtr *nl ){

	FeatureVector *fv = new FeatureVector();
	
	//Add the Bias Feature
	fv->addFeature(1.0);

	//Compute all other features
	std::vector<int>::iterator it = used_feature_idxs.begin();
	for( ; it != used_feature_idxs.end(); ++it )
		featureCogs()[*it].compute( *fv, ion, nl );

	return fv;
}

void FeatureVector::addFeature( double value ){
	if( value != 0.0 ) fv.push_back( fv_idx++ );
	else fv_idx++;
}

void FeatureVector::addFeatureAtIdx( double value, unsigned int idx ){
	if( fv_idx <= idx ) fv_idx = idx + 1;
	if( value != 0.0 ) fv.push_back( idx );
}

//Helper functions for multiple features
const std::vector<std::string> &Feature::OKsymbols() {

	static std::vector<std::string> x;
	static bool initialised = false;

	if (!initialised) {
		x.push_back("Br");
		x.push_back("C");
		x.push_back("Cl");
		x.push_back("F");
		x.push_back("I");
		x.push_back("N");
		x.push_back("O");
		x.push_back("P");
		x.push_back("S");
		x.push_back("Se");
		x.push_back("Si");

		initialised = true;
	}
	return x;
}

const std::vector<std::string> &Feature::OKSymbolsLess() {

	static std::vector<std::string> x;
	static bool initialised = false;

	if (!initialised) {
		x.push_back("C");
		x.push_back("N");
		x.push_back("O");
		x.push_back("P");
		x.push_back("S");
		x.push_back("X");	//For all other

		initialised = true;
	}
	return x;
}

void Feature::replaceUncommonWithX( std::string &symbol ) const{

	//Replace uncommon symbols with X
	if( symbol != "C" && symbol != "N" && symbol != "O" &&
		symbol != "P" && symbol != "S" )
		symbol = "X";
}

//*************************
//FEATURE IMPLEMENTATIONS:
//*************************

void BreakAtomPair::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{
	
	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	std::vector<symbol_pair_t> pairs;

	//Ion Symbol(s)
	std::string irootsymbol, iotherrootsymbol;
	irootsymbol = ion->root->getSymbol();
	replaceUncommonWithX(irootsymbol);
	if( ring_break ){
		iotherrootsymbol = ion->other_root->getSymbol();
		replaceUncommonWithX(iotherrootsymbol);
	}

	//Neutral Loss Symbol(s)
	std::string nlrootsymbol, nlotherrootsymbol;
	nlrootsymbol = nl->root->getSymbol();
	replaceUncommonWithX(nlrootsymbol);
	if( ring_break ){
		nlotherrootsymbol = nl->other_root->getSymbol();
		replaceUncommonWithX(nlotherrootsymbol);
	}

	//Pairs
	pairs.push_back( symbol_pair_t( irootsymbol, nlrootsymbol ));
	if( ring_break ) 
		pairs.push_back( symbol_pair_t( iotherrootsymbol, nlotherrootsymbol ));

	//Iterate through all combinations of atom pairs, appending
	//a feature for each; 1 if it matches, 0 otherwise.
	//Note: the order matters here, ion first then nl
	std::vector<std::string>::const_iterator it1, it2;
	const std::vector<std::string> *ok_symbols = &OKSymbolsLess();
	for( it1 = ok_symbols->begin(); it1 != ok_symbols->end(); ++it1 ){
		for( it2 = ok_symbols->begin(); it2 != ok_symbols->end(); ++it2 ){
			symbol_pair_t sp = symbol_pair_t(*it1, *it2);
			double nonringf = 0.0, ringf = 0.0;
			if( sp == *pairs.begin() ){ 
				nonringf = !ring_break;
				ringf = ring_break;
			}
			if( sp == *pairs.rbegin() ) ringf = ring_break;
			fv.addFeature(nonringf);
			fv.addFeature(ringf);
		} 
	}

}

void RootPathFeature::computeRootPaths(std::vector<path_t> &paths, const RootedROMolPtr *mol, int len, bool ring_break) const{

	path_t path_so_far;
	addPathsFromAtom( paths, mol->root, mol->mol, mol->root, path_so_far, len );
	if( ring_break ) addPathsFromAtom( paths, mol->other_root, mol->mol, mol->other_root, path_so_far, len );

}

void RootPathFeature::addPathsFromAtom( std::vector<path_t> &paths, const RDKit::Atom *atom, const romol_ptr_t mol, const RDKit::Atom *prev_atom, path_t &path_so_far, int len ) const{

	//Add the current symbol
	std::string symbol = atom->getSymbol();
	replaceUncommonWithX(symbol);
	path_so_far.push_back( symbol );

	//Iterate until len is reached, then add the path
	if( len > 1 ){
		RDKit::ROMol::ADJ_ITER_PAIR itp = mol.get()->getAtomNeighbors( atom );
		for( ; itp.first != itp.second; ++itp.first ){
			RDKit::Atom *nbr_atom = mol.get()->getAtomWithIdx(*itp.first);	
			if( nbr_atom != prev_atom ) addPathsFromAtom( paths, nbr_atom, mol, atom, path_so_far, len-1 );
		}
	}
	else paths.push_back( path_so_far );
	
	//Remove the latest symbol
	path_so_far.pop_back();
}

void RootPathFeature::addRootPairFeatures(FeatureVector &fv, std::vector<path_t> &paths, int ring_break) const{

	//Add a feature indicating that there are no pairs
	fv.addFeature( (double)(paths.size() == 0) );

	//Iterate through all combinations of atom pairs, adding a count
	//of the number of times each is seen;
	//Note: the order matters here, root atom first then other
	std::vector<std::string>::const_iterator it1, it2;
	const std::vector<std::string> *ok_symbols = &OKSymbolsLess();
	for( it1 = ok_symbols->begin(); it1 != ok_symbols->end(); ++it1 ){
		for( it2 = ok_symbols->begin(); it2 != ok_symbols->end(); ++it2 ){
			
			//Count how many of each possible symbol pair we have
			path_t sp; sp.push_back(*it1); sp.push_back(*it2);
			std::vector<path_t>::iterator it3 = paths.begin();
			double count = 0.0, ring_count = 0.0;	
			for( ; it3 != paths.end(); ++it3 ){
				if( !ring_break && sp[0] == (*it3)[0] && sp[1] == (*it3)[1])	
					count += 1.0;
				if( ring_break && sp[0] == (*it3)[0] && sp[1] == (*it3)[1])	
					ring_count += 1.0;
			}

			//First feature indicates at least 1
			//Second feature indicates more than 1
			//Non-Ring
			if( count > 0.0 ) fv.addFeature( 1.0 );
			else fv.addFeature( 0.0 );
			if( count > 1.0 ) fv.addFeature( 1.0 );
			else fv.addFeature( 0.0 );
			//Ring
			if( ring_count > 0.0 ) fv.addFeature( 1.0 );
			else fv.addFeature( 0.0 );
			if( ring_count > 1.0 ) fv.addFeature( 1.0 );
			else fv.addFeature( 0.0 );
		} 
	}
}

void RootPathFeature::addRootTripleFeatures(FeatureVector &fv, std::vector<path_t> &paths, int ring_break ) const{
	
	//Add a feature indicating that there are no triples
	fv.addFeature( (double)(paths.size() == 0) );

	//Iterate through all combinations of atom triples, adding a count
	//of the number of times each is seen;
	//Note: the order matters here, root atom first then the next in the path, etc.
	std::vector<std::string>::const_iterator it1, it2, it3;
	const std::vector<std::string> *ok_symbols = &OKSymbolsLess();
	for( it1 = ok_symbols->begin(); it1 != ok_symbols->end(); ++it1 ){
		for( it2 = ok_symbols->begin(); it2 != ok_symbols->end(); ++it2 ){
			for( it3 = ok_symbols->begin(); it3 != ok_symbols->end(); ++it3 ){
				//Count how many of each possible symbol pair we have
				path_t trp; trp.push_back(*it1); trp.push_back(*it2); trp.push_back(*it3);
				std::vector<path_t>::iterator it4 = paths.begin();
				double count = 0.0, ring_count = 0.0;
				for( ; it4 != paths.end(); ++it4 ){
					if( trp[0] == (*it4)[0] && trp[1] == (*it4)[1] && trp[2] == (*it4)[2] ){
						if( !ring_break ) count += 1.0;
						else ring_count += 1.0;
					}
				}

				//First feature indicates at least 1
				//Second feature indicates more than 1
				//Non-Ring
				if( count > 0.0 ) fv.addFeature( 1.0 );
				else fv.addFeature( 0.0 );
				if( count > 1.0 ) fv.addFeature( 1.0 );
				else fv.addFeature( 0.0 );
				//Ring
				if( ring_count > 0.0 ) fv.addFeature( 1.0 );
				else fv.addFeature( 0.0 );
				if( ring_count > 1.0 ) fv.addFeature( 1.0 );
				else fv.addFeature( 0.0 );
			}
		} 
	}

}

void IonRootPairs::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const{

	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	std::vector<path_t> paths;
	computeRootPaths( paths, ion, 2, ring_break );
	addRootPairFeatures( fv, paths, ring_break );
}

void NLRootPairs::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{

	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	std::vector<path_t> paths;
	computeRootPaths( paths, nl, 2, ring_break );
	addRootPairFeatures( fv, paths, ring_break);
}

void IonRootTriples::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{

	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	std::vector<path_t> paths;
	computeRootPaths( paths, ion, 3, ring_break );
	addRootTripleFeatures( fv, paths, ring_break);
}

void NLRootTriples::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{

	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	std::vector<path_t> paths;
	computeRootPaths( paths, nl, 3, ring_break );
	addRootTripleFeatures( fv, paths, ring_break);
}

void GasteigerCharges::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{
	
	//Collect the charges from the root atoms
	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	typedef std::pair<double, double> gasteiger_t;
	std::vector<gasteiger_t> gasteigers;
	
	//Ion
	double icharge, iothercharge;
	ion->root->getProp<double>("OrigGasteigerCharge", icharge);
	if( ring_break ) ion->other_root->getProp<double>("OrigGasteigerCharge", iothercharge);

	//Neutral Loss
	double nlcharge, nlothercharge;
	nl->root->getProp<double>("OrigGasteigerCharge", nlcharge);
	if( ring_break ) nl->other_root->getProp<double>("OrigGasteigerCharge", nlothercharge);

	//Collate the charges
	gasteigers.push_back( gasteiger_t( icharge, nlcharge ) );
	if( ring_break ) gasteigers.push_back( gasteiger_t( iothercharge, nlothercharge ) );

	if( !ring_break ){ 
				
		//Then there are 6 x 6 = 36 possible configurations of the charge
		// - Allocate one bit to each
		int gc_ion = discretizeGasteigerCharge( gasteigers[0].first );
		int gc_nl = discretizeGasteigerCharge( gasteigers[0].second );
		for( int i = 0; i <= 5; i++ ){
			for( int j = 0; j <= 5; j++ ){
				if( i == gc_ion && j == gc_nl ) fv.addFeature( 1.0 );
				else fv.addFeature( 0.0 );	
			}
		}
		//Ring break charges
		for( int i = 0; i < 36; i++ ) fv.addFeature( 0.0 );	
	}
	else{
		//Non-Ring charges 
		for( int i = 0; i < 36; i++ ) fv.addFeature( 0.0 );	

		//Ring Charges - set bit if either breaks fit the rule
		int gc0_ion = discretizeGasteigerCharge( gasteigers[0].first );
		int gc0_nl = discretizeGasteigerCharge( gasteigers[0].second );
		int gc1_ion = discretizeGasteigerCharge( gasteigers[1].first );
		int gc1_nl = discretizeGasteigerCharge( gasteigers[1].second );

		for( int i = 0; i <= 5; i++ ){
			for( int j = 0; j <= 5; j++ ){
				if( i == gc0_ion && j == gc0_nl ) fv.addFeature( 1.0 );
				else if( i == gc1_ion && j == gc1_nl ) fv.addFeature( 1.0 );
				else fv.addFeature( 0.0 );	
			}
		}
	}
}

int GasteigerCharges::discretizeGasteigerCharge( double gc ) const{

	//Discretize the Gasteiger Charges into 6 levels:
	// x < -0.5
	// -0.5 <= x < -0.1
	// -0.1 <= x < 0
	// 0 <= x < 0.1
	// 0.1 <= x <= 0.5
	// 0.5 <= x

	if( gc < -0.5 ) return 0;
	else if( gc < -0.1 ) return 1;
	else if( gc < 0 ) return 2;
	else if( gc < 0.1 ) return 3;
	else if( gc < 0.5 ) return 4;
	else return 5;
}

void HydrogenMovement::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl  ) const{

	double h_movement = 0.0;

	//Compute the mass difference in the ion
	RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
	RDKit::ROMol::AtomIterator ai;
	for( ai = ion->mol.get()->beginAtoms(); ai != ion->mol.get()->endAtoms(); ++ai ){
		double orig_mass, mass = 0.0;
		std::string symbol = (*ai)->getSymbol();
		mass += pt->getMostCommonIsotopeMass(symbol);
		mass += (*ai)->getTotalNumHs()*pt->getMostCommonIsotopeMass("H");
		(*ai)->getProp<double>("OriginalMass", orig_mass);
		h_movement += (mass - orig_mass);
	}

	//Binary on/off indicating whether a particular transfer occurred
	for( double h = -4.0; h <= 4.0; h += 1.0 ){	
		if( fabs( h - h_movement ) < 0.5 ) fv.addFeature( 1.0 );
		else fv.addFeature(0.0);
	}
	//Catch-all for all other hydrogen transfers
	if( fabs( h_movement ) > 4.0 ) fv.addFeature( 1.0 );
	else fv.addFeature( 0.0 );
}

void RingFeatures::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const{
	
	//Non-Ring break
	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	if( !ring_break){
		for( int i = 0; i < 12; i++ ) fv.addFeature(0);
	}
	else{
	//Ring Break
		//Find a broken bond, and fetch the aromaticity labels from it
		int is_arom, is_dbl_arom;
		nl->mol.get()->getProp( "IsAromaticRingBreak", is_arom );
		nl->mol.get()->getProp( "IsAromaticDblRingBreak", is_dbl_arom );

		int root_dist_ion = calcRootDistance( ion );
		int root_dist_nl = calcRootDistance( nl );

		int break_dist = std::min( root_dist_ion, root_dist_nl) + 1;
		int ring_size = root_dist_nl + root_dist_ion + 2;

		fv.addFeature(!is_arom);
		fv.addFeature(is_arom);
		fv.addFeature(is_dbl_arom);
		for( int dist = 1; dist <=3; dist++ )
			fv.addFeature( break_dist == dist );
		fv.addFeature( break_dist > 3 );
		for( int rsize = 3; rsize <= 6; rsize++ )
			fv.addFeature( rsize == ring_size );
		fv.addFeature( ring_size > 6 );
	}

}

int RingFeatures::calcRootDistance(const RootedROMolPtr *mol) const{
	
	//Compute the distance from the root to the other root (breadth-first-search)
	typedef std::pair<RDKit::Atom *,int> item_t;
	std::set<int> seen_idxs;
	std::queue<item_t> queue;
	queue.push( item_t( mol->root, 0 ) );
	while( queue.size() > 0 ){
		item_t item = queue.front();
		RDKit::Atom *at = item.first;
		int depth = item.second;
		seen_idxs.insert(at->getIdx());

		//Have we found the root?
		if( at == mol->other_root ) return depth;

		//If not, keep looking
		RDKit::ROMol::ADJ_ITER_PAIR itp = mol->mol.get()->getAtomNeighbors( at );
		for( ; itp.first != itp.second; ++itp.first ){
			RDKit::Atom *nbr_atom = mol->mol.get()->getAtomWithIdx(*itp.first);
			if( seen_idxs.find( nbr_atom->getIdx() ) == seen_idxs.end() ) 
				queue.push( item_t(nbr_atom, depth+1) );
		}
		queue.pop();
	}
	return -1; //OtherRoot not found
}


void QuadraticFeatures::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{

	//Compute quadratic feature indexes for all existing features
	int n = fv.getTotalLength();
	std::vector<feature_t>::const_iterator it1, it2;
	std::vector<int> quadratic_indexes;
	for( it1 = fv.getFeatureBegin() + 1; it1 != fv.getFeatureEnd(); ++it1 ){
		//Due to symmetry and not wanting to include square features,
		//or bias features, only the lower left triangle of the 
		//feature x feature matrix is included (minus bias row/col)
		//- the offset gives the index for the first used feature in each row.
		int offset = n + (*it1 - 2)*(*it1 - 1)/2;	
		for( it2 = fv.getFeatureBegin() + 1; it2 != it1; ++it2 )
			quadratic_indexes.push_back( offset + *it2 - 1 );
	}
	//Add the features 
	//Note: modifying the feature vector in the above loop causes problems...
	std::vector<int>::iterator it = quadratic_indexes.begin();
	for( ; it != quadratic_indexes.end(); ++it )
		fv.addFeatureAtIdx(1.0, *it);

	//Update the feature length (without modifying any features).
	int total_num_features = n + (n-1)*(n-2)/2;
	fv.addFeatureAtIdx(0.0, total_num_features-1); 
}