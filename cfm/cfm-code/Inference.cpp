/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# inference.cpp
#
# Description: 	Custom implementation of jtree based exact inference for
#				this exact network (for efficiency).
#				
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Inference.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>

#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions.hpp>

//Big log negative value that will effectively result in zero probability
static const double NULL_PROB = -1000000000000.0;
static const int DOWN = 0;
static const int UP = 1;

void createMessage( factor_probs_t &tmp_log_probs, Message &m, Message &prev_m, int direction, int depth, MolData *moldata );
void passMessage(factor_probs_t &tmp_log_probs, int direction, int depth, MolData *moldata, Message &m, int energy );

void initTmpFactorProbSizes( factor_probs_t &tmp_log_probs, unsigned int num_frag, unsigned int num_trans, unsigned int model_depth ){
	
	//Initialise persistence vector sizes
	tmp_log_probs.ps.resize( num_frag );
	for( unsigned int i = 0; i < num_frag; i++ ) 
		tmp_log_probs.ps[i].resize( model_depth );
	
	//Initialise transition vector sizes
	tmp_log_probs.tn.resize( num_trans );
	for( unsigned int i = 0; i < num_trans; i++ ) 
		tmp_log_probs.tn[i].resize( model_depth );
}


void runInferenceDownwardPass( std::vector<Message> &down_msgs, MolData *moldata, config_t &cfg, int to_depth){

	const FragmentGraph *fg = moldata->getFragmentGraph();

	//Initialise the messages
	down_msgs.resize(cfg.model_depth);

	//Create the tmp factor probs
	factor_probs_t tmp_log_probs;
	initTmpFactorProbSizes( tmp_log_probs, fg->getNumFragments(), fg->getNumTransitions(), cfg.model_depth );
	
	//Factor (F0,F1) => Create F1 Message
	int energy = cfg.map_d_to_energy[0];
	const std::vector<int> *tmap = &((*fg->getFromIdTMap())[0]);
	std::vector<int>::const_iterator it = tmap->begin();
	down_msgs[0].reset(fg->getNumFragments());
	down_msgs[0].addToIdx(0, moldata->getLogPersistenceProbForIdx(energy,0));
	for( ; it != tmap->end(); ++it ){
		const Transition *t = fg->getTransitionAtIdx(*it);
		down_msgs[0].addToIdx(t->getToId(), moldata->getLogTransitionProbForIdx(energy,*it));
	}

	//Update Factor (F1,F2) => Create Message F2 => Update Factor (F2,F3) ...etc as per MODEL_DEPTH
	for( int i = 0; i < to_depth-1; i++ ){	
		energy = cfg.map_d_to_energy[i+1];
		passMessage(tmp_log_probs, DOWN, i, moldata, down_msgs[i], energy );
		createMessage( tmp_log_probs, down_msgs[i+1], down_msgs[i], DOWN, i, moldata );
	}

}

void createMessage( factor_probs_t &tmp_log_probs, Message &m, Message &prev_m, int direction, int depth, MolData *moldata ){
	
	const FragmentGraph *fg = moldata->getFragmentGraph();
	m.reset(fg->getNumFragments());
	for( unsigned int id = 0; id < fg->getNumFragments(); id++ ){

		//Going up or down?
		const std::vector<int> *tmap;
		if( direction == DOWN ) tmap = &((*fg->getToIdTMap())[id]);
		else tmap = &((*fg->getFromIdTMap())[id]);
		
		//Marginalize out the upper/lower variable
		double log_sum = -DBL_MAXIMUM;
		if( prev_m.getIdx(id) > -DBL_MAXIMUM ) log_sum = tmp_log_probs.ps[id][depth];

		std::vector<int>::const_iterator itt = tmap->begin();
		for( ; itt != tmap->end(); ++itt ){
			const Transition *t = fg->getTransitionAtIdx(*itt);
			if( (direction == DOWN && prev_m.getIdx(t->getFromId()) > -DBL_MAXIMUM ) ||
			    (direction == UP && prev_m.getIdx(t->getToId()) > -DBL_MAXIMUM )){
				log_sum = logAdd(log_sum, tmp_log_probs.tn[*itt][depth]);
			}
		}
		if( log_sum > -DBL_MAXIMUM ) m.addToIdx(id, log_sum);
	}
}

void passMessage(factor_probs_t &tmp_log_probs, int direction, int depth, MolData *moldata, Message &m, int energy ){

	const FragmentGraph *fg = moldata->getFragmentGraph();
	Message::const_iterator it = m.begin();
	for( ; it != m.end(); ++it ){

		unsigned int idx = it.index();

		//Apply to the persistence term (idx -> idx)
		tmp_log_probs.ps[idx][depth] = moldata->getLogPersistenceProbForIdx(energy,idx) + m.getIdx(idx);
	
		//Apply to all the other transitions applicable for this message element
		const std::vector<int> *tmap;
		if( direction == DOWN ) tmap = &((*fg->getFromIdTMap())[idx]);
		else tmap = &((*fg->getToIdTMap())[idx]);
		std::vector<int>::const_iterator itt = tmap->begin();
		for( ; itt != tmap->end(); ++itt )
			tmp_log_probs.tn[*itt][depth] = moldata->getLogTransitionProbForIdx(energy, *itt) + m.getIdx(idx);
	}
}




