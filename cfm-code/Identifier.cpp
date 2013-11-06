/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Identifier.cpp
#
# Description: 	Identifier class for ranking candidate structures according
#				to their match with a set of target spectra.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Identifier.h"

#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

static bool sort_candidates(const Candidate &u, const Candidate &v){
	return u.getScore() > v.getScore();
}

//Ranks the list of candidates according to the match between their predicted spectra and the target
void Identifier::rankCandidatesForSpecMatch( std::vector<Candidate> &candidates, const std::vector<Spectrum> *target_spectra ){

	//Compute the scores for each candidate
	std::vector<Candidate>::iterator it = candidates.begin();
	for( ; it != candidates.end(); ++it ){
	
		double score = 0.0;
		try{


			//Create the MolData structure with the input
			MolData moldata( it->getId()->c_str(), it->getSmilesOrInchi()->c_str() );
	
			//Calculate the pruned FragmentGraph
			moldata.computeLikelyFragmentGraphAndSetThetas(param, cfg, prob_thresh_for_prune);

			//Predict the spectra (and post-process, use existing thetas)
			moldata.computePredictedSpectra( *param, *cfg, true, true );

			//Compute the score for the current candidate:
			// - The score is the sum of the comparison scores between all the target and predicted spectra
			for( unsigned int energy = 0; energy < target_spectra->size(); energy++ )
				score += cmp->computeScore( &((*target_spectra)[energy]), moldata.getPredictedSpectrum(energy) );
			
		}
		catch( RDKit::MolSanitizeException e ){
			std::cout << "Could not sanitize " << *it->getSmilesOrInchi() << std::endl;
		}
		catch( RDKit::SmilesParseException pe ){
			std::cout << "Could not parse " << *it->getSmilesOrInchi() << std::endl;		
		}

		it->setScore(score);
	}

	//Sort the candidates by decreasing score
	std::sort( candidates.begin(), candidates.end(), sort_candidates );

}

