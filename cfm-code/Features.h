/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Features.h
#
# Description: 	Code for computing features for fragmentations.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __FEATURE_H__
#define __FEATURE_H__

#include "Util.h"

#include <boost/ptr_container/ptr_vector.hpp>

#include <string>
#include <vector>

struct input_file_t;

//Exception to throw when the input feature configuration file is invalid 
class InvalidConfigException: public std::exception{

	virtual const char* what() const throw(){
		return "Invalid Configuration File";
	}
};

//Structure to hold a sparse computed feature vector
typedef unsigned int feature_t;
class FeatureVector{
public:
	FeatureVector(){fv_idx = 0;};
	void addFeature( double value );
	void addFeatureAtIdx( double value, unsigned int idx );
	unsigned int const getTotalLength() const {return fv_idx;};
	feature_t const getFeature( int idx ) const { return fv[idx];};
	std::vector<feature_t>::const_iterator getFeatureBegin() const { return fv.begin(); };
	std::vector<feature_t>::const_iterator getFeatureEnd() const { return fv.end(); };
	unsigned int const getNumSetFeatures() const { return fv.size();};

private:
	std::vector<feature_t> fv;
	unsigned int fv_idx;
};

//Base class to compute a feature - all features should inherit from this
class Feature{

public:
	virtual void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const = 0;
	unsigned int getSize() const { return size; };
	std::string getName() const { return name; };
	virtual ~Feature(){};

protected:
	unsigned int size;	
	std::string name;
	static const std::vector<std::string> &OKsymbols();
	static const std::vector<std::string> &OKSymbolsLess();
	void replaceUncommonWithX( std::string &symbol ) const;
};

//Class to compute a feature vector
class FeatureCalculator{

public:
	//Constructor: Initialise the calculator using a config file listing features
	FeatureCalculator( std::string &config_filename );

	//Constructor: Initialise the calculator using a list of feature names
	FeatureCalculator( std::vector<std::string> &feature_list );
	
	//Compute the expected number of total features
	unsigned int getNumFeatures();
	
	//Retrieve the list of feature names being used
	std::vector<std::string> getFeatureNames();

	//Retrieve a list of valid feature names (for testing)
	static const std::vector<std::string> getValidFeatureNames();

	//Compute the feature vector for the input ion and nl (with labeled Root atoms)
	// - NB: responsibility of caller to delete.
	FeatureVector *computeFV( const RootedROMolPtr *ion, const RootedROMolPtr *nl );

private:
	//List of feature classes ready to be used
	static const boost::ptr_vector<Feature> &featureCogs();

	//Indexes of feature classes that are selected for use
	std::vector<int> used_feature_idxs;

	//Helper function - Configure feature for use
	void configureFeature( std::string &name );
};

//*************************
//FEATURE IMPLEMENTATIONS:
//*************************

class BreakAtomPair : public Feature {
public:
	BreakAtomPair(){ size = 72; name = "BreakAtomPair"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};

class RootPathFeature : public Feature {
protected:
	typedef std::vector<std::string> path_t;
	void computeRootPaths(std::vector<path_t> &paths, const RootedROMolPtr *mol, int len, bool ring_break) const;
	void addRootPairFeatures(FeatureVector &fv, std::vector<path_t> &paths, int ring_break) const;
	void addRootTripleFeatures(FeatureVector &fv, std::vector<path_t> &paths, int ring_break) const;
private:
	void addPathsFromAtom( std::vector<path_t> &paths, const RDKit::Atom *atom, const romol_ptr_t mol, const RDKit::Atom *prev_atom, path_t &path_so_far, int len ) const;
};

class IonRootPairs : public RootPathFeature {
public:
	IonRootPairs(){ size = 145; name = "IonRootPairs"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLRootPairs : public RootPathFeature {
public:
	NLRootPairs(){ size = 145; name = "NLRootPairs"; }; 
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class IonRootTriples : public RootPathFeature {
public:
	IonRootTriples(){ size = 865; name = "IonRootTriples"; }; 
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLRootTriples : public RootPathFeature {
public:
	NLRootTriples(){ size = 865; name = "NLRootTriples"; }; 
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class GasteigerCharges : public RootPathFeature {
public:
	GasteigerCharges(){ size = 72; name = "GasteigerCharges"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl)  const;
private:
	int discretizeGasteigerCharge( double gc ) const;
};

class HydrogenMovement : public Feature {
public:
	HydrogenMovement(){ size = 10; name = "HydrogenMovement"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class RingFeatures : public Feature {
public:
	RingFeatures(){ size = 12; name = "RingFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
private:
	//Helper function - compute the distance between two root 
	//atoms in a molecule (assumes ring break)
	int calcRootDistance(const RootedROMolPtr *mol)  const;
};

class QuadraticFeatures : public Feature {
public:
	QuadraticFeatures(){ size = 0; name = "QuadraticFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};


#endif // __FEATURE_H__