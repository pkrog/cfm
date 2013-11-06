/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# em_test.h
#
# Description: 	Test code for EM.cpp
#
# Author: Felicity Allen
# Created: August 2013
#########################################################################*/

#ifndef __EM_TEST_H__
#define __EM_TEST_H__

#include "../test.h"
#include "EM.h"

class EMTestSelfProduction : public Test {
public:
	EMTestSelfProduction();
	void runTest();
};

class EMTestSingleEnergySelfProduction : public Test {
public:
	EMTestSingleEnergySelfProduction();
	void runTest();
};

class EMTestMultiProcessor: public Test {
public:
	EMTestMultiProcessor();
	void runTest();
};

#endif // __EM_TEST_H__