/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# MolData.h
#
# Description: 	Class to hold the input data belonging to a molecule:
#					 - An ID and smiles/inchi
#					 - (optional) A cross-validation group identifier
#					 - (optional) A computed fragmentation graph
#					 - (optional) A computed set of features corresponding to that graph
#					 - (optional) A computed set of theta values for that graph
#				     - (optional) A computed set of transition probabilities using those thetas.
#					 - (optional) A set of spectra
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __MOLDATA_H__
#define __MOLDATA_H__

#include "Features.h"
#include "FragmentGraph.h"
#include "FragmentGraphGenerator.h"
#include "Param.h"
#include "Config.h"
#include "Message.h"

typedef std::pair<int,double> annotation_t; //<Fragment Id, Score>

static bool sort_annotations_by_score(const annotation_t &u, const annotation_t &v){
   return u.second > v.second;
}

class Peak{
public:
	Peak(){};
	Peak(double a_mass, double an_intensity) : 
	  mass(a_mass), intensity(an_intensity) {};
	double mass;
	double intensity;
	std::vector<annotation_t> annotations; 
};

typedef std::vector<Peak> Spectrum;

static bool sort_peaks_by_intensity(const Peak &u, const Peak &v){
   return u.intensity > v.intensity;
}

static bool sort_peaks_by_mass(const Peak &u, const Peak &v){
   return u.mass < v.mass;
}

class MolData{
public:
	MolData( std::string &an_id, std::string &an_smiles_or_inchi, int a_group ) :
	  id(an_id), smiles_or_inchi(an_smiles_or_inchi), group(a_group), pos_graph_computed(0), pos_ev_graph_computed(0), neg_graph_computed(0), neg_ev_graph_computed(0) {}; 
	MolData( std::string &an_id, std::string &an_smiles_or_inchi ) :
	  id(an_id), smiles_or_inchi(an_smiles_or_inchi), group(0), pos_graph_computed(0), pos_ev_graph_computed(0), neg_graph_computed(0), neg_ev_graph_computed(0) {}; 
	MolData( const char *an_id, const char *an_smiles_or_inchi ) :
	  id(an_id), smiles_or_inchi(an_smiles_or_inchi), group(0), pos_graph_computed(0), pos_ev_graph_computed(0), neg_graph_computed(0), neg_ev_graph_computed(0) {}; 
		
	//Access functions
	const FragmentGraph *getFragmentGraph(){if(using_negative_mode) return neg_fg; else return pos_fg;};
	const EvidenceFragmentGraph *getEvidenceFragmentGraph(){if(using_negative_mode) return neg_ev_fg; else return pos_ev_fg;};
	const Spectrum *getSpectrum(int energy){return &(spectra[energy]);};
	const std::vector<Spectrum> *getSpectra(){return &spectra;};
	const Spectrum *getPredictedSpectrum(int energy){return &(predicted_spectra[energy]);};
	unsigned int getNumSpectra(){return spectra.size();};
	const FeatureVector *getFeatureVectorForIdx(int index){ return (*fvs)[index]; };
	double getThetaForIdx(int energy, int index)
	{ return (*thetas)[energy][index];};
	double getLogTransitionProbForIdx(int energy, int index)
	{ return (*log_probs)[energy][index];};
	double getLogPersistenceProbForIdx(int energy, int index)
	{ return (*log_probs)[energy][(*fg)->getNumTransitions() + index];};
	int getGroup(){ return group; };
	std::string getId(){ return id; };
	std::string getSmilesOrInchi(){ return smiles_or_inchi; };

	void readInSpectraFromFile( const std::string &filename, bool readToPredicted = false );
	void removePeaksWithNoFragment( double abs_tol, double ppm_tol );
	void writePredictedSpectraToFile( std::string &filename );
	void writeFullEnumerationSpectrumToFile( std::string &filename );
	void outputSpectra( std::ostream &out, const char*spec_type, bool do_annotate = false );
	void createInterpolatedSpectra( config_t &cfg );
	
	//More memory efficient alternative to calling computeFragmentGraph and 
	//then computeFeatureVectors with deleteMols = true
	void computeFragmentGraphAndReplaceMolsWithFVs( int depth, FeatureCalculator *fc, bool retain_smiles = false);

	//Replaces computeFragmentGraph, computeFeatureVectors and computeTransitionThetas
	//below (delteMols = true), pruning according to prob_thresh_for_prune value.
	void computeLikelyFragmentGraphAndSetThetas( Param *param, config_t *cfg, double prob_thresh_for_prune, bool retain_smiles = false );

	//Note that the following should be called in this order
	//since each one assumes all previous have already been called.
	void computeFragmentGraph( int depth );
	void computeFeatureVectors( FeatureCalculator *fc, bool deleteMols = false );
	void computeTransitionThetas( Param &param );
	void computeTransitionProbabilities();
	void computePredictedSpectra( Param &param, config_t &cfg, bool postprocess = false, bool use_existing_thetas = false );
	
	//Function to compute a much reduced fragment graph containing only those
	//fragmentations as actually occur in the spectra, based on a computed set of beliefs
	//thresholding inclusion in the graph by the provided belief_thresh value (log domain)
	void computeEvidenceFragmentGraph( beliefs_t *beliefs, double log_belief_thresh, config_t *cfg );
	void annotatePeaks(double abs_tol, double ppm_tol, bool prune_deadends = true);
	void setIonizationMode(bool is_negative_mode);
	~MolData();

protected:	//These features are protected rather than private for access during tests.
	int group;
	std::string id;
	std::string smiles_or_inchi;

	//Fragment Graphs
	FragmentGraph **fg; //Working fg -> set by setIonizationMode()
	FragmentGraph *pos_fg;
	FragmentGraph *neg_fg;
	EvidenceFragmentGraph **ev_fg; //Working ev_fg -> set by setIonizationMode()
	EvidenceFragmentGraph *pos_ev_fg;
	EvidenceFragmentGraph *neg_ev_fg;
	
	//Flags indicating when the respective fragment graphs have been computed
	int *graph_computed;
	int *ev_graph_computed;
	int pos_graph_computed;
	int pos_ev_graph_computed;
	int neg_graph_computed;
	int neg_ev_graph_computed;
	
	//Spectra
	std::vector<Spectrum> spectra;
	std::vector<Spectrum> predicted_spectra;
	
	//Model values
	std::vector<FeatureVector *> *fvs;
	std::vector<std::vector<double> > *thetas;
	std::vector<std::vector<double> > *log_probs;
	std::vector<FeatureVector *> pos_fvs;
	std::vector<std::vector<double> > pos_thetas;
	std::vector<std::vector<double> > pos_log_probs;
	std::vector<FeatureVector *> neg_fvs;
	std::vector<std::vector<double> > neg_thetas;
	std::vector<std::vector<double> > neg_log_probs;

	//Ionization mode used in current computations
	bool using_negative_mode;
	
	void computeGraphWithGenerator( FragmentGraphGenerator &fgen, int depth );
	void interpolateSpectra( Spectrum &output, Spectrum &lower_spec, Spectrum &higher_spec, double ratio, config_t &cfg );
	void postprocessPredictedSpectra();

	void computePredictedSingleEnergySpectra( Param &param, config_t &cfg, bool postprocess, bool use_existing_thetas );
	void translatePeaksFromMsgToSpectra( Spectrum &out_spec, Message *msg );
	void computeFragmentEvidenceValues(std::vector<double> &evidence, int frag_idx, const beliefs_t *beliefs, config_t *cfg );
	static void postprocessSpectrum( Spectrum &spectrum );
	static void normalizeAndSortSpectrum( Spectrum &spectrum );
	static void sortAndNormalizeAnnotations(Spectrum &spectrum);

};

#endif // __MOLDATA_H__
