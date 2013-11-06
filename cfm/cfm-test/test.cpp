/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# test.cpp
#
# Description: 	Test code.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "mpi.h"
#include "test.h"

#include <iostream> 
#include <exception>
#include <vector>

#include "tests/message_test.h"
#include "tests/features_test.h"
#include "tests/param_test.h"
#include "tests/ipfp_test.h"
#include "tests/comms_test.h"
#include "tests/em_test.h"
#include "tests/comparator_test.h"

#include <boost/ptr_container/ptr_vector.hpp>

Test::Test(){
	description = "None";
}

void Test::run(){
	try{
		excpt_occurred = false;
		passed = false;
		runTest();
	}
	catch(std::exception e){
		excpt_occurred = true;
		std::cout << e.what();
		throw e;
	}
}

void Test::runTest(){
	std::cout << "Running parent runTest" << std::endl;
	//Do nothing: placeholder function
}

void runTests(boost::ptr_vector<Test> &tests){

	int num_passed = 0;
	int num_failed = 0;
	int num_run = 0;
	int num_exceptions = 0;

	boost::ptr_vector<Test>::iterator it = tests.begin();
	for( ; it != tests.end(); ++it ){
		std::cout << it->description << "...";
		it->run();
		if( it->excpt_occurred ){ 
			num_exceptions++;
			std::cout << "EXCEPTION" << std::endl;
		}else if( it->passed ){ 
			num_passed ++;
			std::cout << "PASS" << std::endl;
		}
		else{ 
			std::cout << "FAIL" << std::endl;
			num_failed++;
		}
		num_run++;
	}

	std::cout << std::endl << "----------------------------------------------" << std::endl;
	std::cout << num_run << " Tests Run" << std::endl;
	std::cout << num_passed << " Passed, ";
	std::cout << num_failed << " Failed, ";
	std::cout << num_exceptions << " Exceptions";
	std::cout << std::endl << "----------------------------------------------" << std::endl;

}

int main(int argc, char *argv[])
{
	boost::ptr_vector<Test> tests;
    MPI_Init( &argc, &argv );
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump == 1 ){
		tests.push_back( new ComparatorsTestRecall() );
		tests.push_back( new ComparatorsTestPrecision() );
		tests.push_back( new ComparatorsTestJaccard() );
		tests.push_back( new MessageTestBasics() );
		tests.push_back( new MessageTestIteration() );
		tests.push_back( new MessageTestCopyAssignment() );
		tests.push_back( new IPFPTestComputeBeliefsSingleEnergy() );
		tests.push_back( new IPFPTestComputeBeliefsConverge() );
		//tests.push_back( new IPFPTestComputeBeliefsInterpolate() );
		//tests.push_back( new IPFPTestComputeBeliefsInterpolateNoDirect() );
		tests.push_back( new IPFPTestComputeBeliefsNonConverge() );
		tests.push_back( new IPFPTestComputeBeliefsSharedMass() );
		tests.push_back( new FeaturesTestInit() );
		tests.push_back( new FeaturesTestBreakAtomPair() );
		tests.push_back( new FeaturesTestRootPairs() );
		tests.push_back( new FeaturesTestRootTriples() );
		tests.push_back( new FeaturesTestGasteigerCharges() );
		tests.push_back( new FeaturesTestHydrogenMovement() );
		tests.push_back( new FeaturesTestRingFeatures() );
		tests.push_back( new FeaturesTestQuadraticFeatures() );
		tests.push_back( new FeaturesTestMetlinExample() );
		tests.push_back( new FeaturesTestLength() );
		tests.push_back( new ParamsTestComputeTransitionThetas() );
		tests.push_back( new ParamsTestComputeAndAccumulateGradient() );
		tests.push_back( new EMTestSelfProduction() );
		tests.push_back( new EMTestSingleEnergySelfProduction() );
	}
	if( mpi_nump > 2 ){
		tests.push_back( new CommsTestSetMasterUsedIdxs() );
		tests.push_back( new CommsTestCollectQInMaster() );
		tests.push_back( new CommsTestCollectGradsInMaster() );
		tests.push_back( new CommsTestBroadcastParams() );
		tests.push_back( new EMTestMultiProcessor() );
	}

	runTests( tests );

	MPI_Finalize();

	return(0);    
}



