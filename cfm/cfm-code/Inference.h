/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# inference.h
#
# Description: 	Custom implementation of jtree based exact inference for
#				this network (for efficiency).
#				
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __INFER_H__
#define __INFER_H__

#include "MolData.h"
#include "Message.h"
#include "Config.h"

struct factor_probs_t{
	std::vector<std::vector<double> > tn;	//Transition:  Indexed by transition, then depth
	std::vector<std::vector<double> > ps;	//Persistence: Indexed by fragment, then depth
};

void runInferenceDownwardPass( std::vector<Message> &down_msgs, MolData *moldata, config_t &cfg, int to_depth);
void completeInferenceUpwardsPass( 	Message &up_m, std::vector<Message> &up_msgs, std::vector<Message> &down_msgs, std::vector<Message> &spec_msgs, MolData *moldata, config_t &cfg, int to_depth );

#endif // __INFER_H__