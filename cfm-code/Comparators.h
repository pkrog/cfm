/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Comparators.h
#
# Description: 	Functions for comparing spectra and returning a score.
#				So far, just includes the dot product as described in Stein 1994.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __COMPARATORS_H__
#define __COMPARATORS_H__

#include "MolData.h"

typedef std::pair<Peak, Peak> peak_pair_t;

//Base class to compute a spectrum comparison score
class Comparator{
public:
	Comparator( double a_ppm_tol, double a_abs_tol ) : ppm_tol(a_ppm_tol), abs_tol( a_abs_tol ) {};
	virtual double computeScore( const Spectrum *measured, const Spectrum *predicted ) const = 0;
protected:
	double ppm_tol;
	double abs_tol;

	void getMatchingPeakPairs( std::vector<peak_pair_t> &peak_pairs,  const Spectrum *p, const Spectrum *q ) const;
};

class DotProduct : public Comparator {
public:
	DotProduct( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const;
private:
	double getTotalPeakSum( const Spectrum *spectrum ) const;
	double getAdjustedIntensity( double intensity, double mass ) const;
};

class Precision : public Comparator {
public:
	Precision( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const;
private:
};

class WeightedRecall : public Comparator {
public:
	WeightedRecall( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const;
private:
};

class Jaccard : public Comparator {
public:
	Jaccard( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const;
private:
};


#endif // __COMPARATORS_H__