/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# comparator_test.cpp
#
# Description: Test code for Comparators.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/

#include "mpi.h"

#include "comparator_test.h"
#include "MolData.h"

#include <boost/filesystem.hpp>


ComparatorsTestRecall::ComparatorsTestRecall(){
	description = "Test computing of recall score";
}

void ComparatorsTestRecall::runTest(){
	
	bool pass = true;
    double tol = 1e-10;

	WeightedRecall cmp(10.0, 0.01);
	std::vector<double> scores(3);

	//Test with exact matches
	MolData moldata("Test1","C");
	std::string specfile = "tests/test_data/test_spec/Test1.txt";
	moldata.readInSpectraFromFile(specfile);
	std::string pspecfile = "tests/test_data/test_pspec/Test1.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i) );
	if( fabs( scores[0] - 100.0 ) > tol || fabs( scores[1] - 100.0 ) > tol || fabs( scores[2] - 100.0 ) > tol ){
		std::cout << "Perfect Matches: Expecting 100.0 scores but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	//Test with no matches
	moldata = MolData("Test2","C");
	specfile = "tests/test_data/test_spec/Test2.txt";
	moldata.readInSpectraFromFile(specfile);
	pspecfile = "tests/test_data/test_pspec/Test2.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i) );
	if( fabs( scores[0] - 0.0 ) > tol || fabs( scores[1] - 0.0 ) > tol || fabs( scores[2] - 0.0 ) > tol ){
		std::cout << "No Match: Expecting 0.0 scores but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	//Test between
	moldata = MolData("Test3","C");
	specfile = "tests/test_data/test_spec/Test3.txt";
	moldata.readInSpectraFromFile(specfile);
	pspecfile = "tests/test_data/test_pspec/Test3.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	double exp_scores[3] = {100.0, 65.0, 50.0};
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i) );
	if( fabs( scores[0] - exp_scores[0] ) > tol || fabs( scores[1] - exp_scores[1] ) > tol || fabs( scores[2] - exp_scores[2] ) > tol ){
		std::cout << "Mixed: Expecting " << exp_scores[0] << " " << exp_scores[1] << " " << exp_scores[2] << " but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	passed = pass;

}

ComparatorsTestPrecision::ComparatorsTestPrecision(){
	description = "Test computing of precision score";
}

void ComparatorsTestPrecision::runTest(){
	
	bool pass = true;
    double tol = 1e-5;

	Precision cmp(10.0, 0.01);
	std::vector<double> scores(3);

	//Test with exact matches
	MolData moldata("Test1","C");
	std::string specfile = "tests/test_data/test_spec/Test1.txt";
	moldata.readInSpectraFromFile(specfile);
	std::string pspecfile = "tests/test_data/test_pspec/Test1.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i) );
	if( fabs( scores[0] - 100.0 ) > tol || fabs( scores[1] - 100.0 ) > tol || fabs( scores[2] - 100.0 ) > tol ){
		std::cout << "Perfect Matches: Expecting 100.0 scores but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	//Test with no matches
	moldata = MolData("Test2","C");
	specfile = "tests/test_data/test_spec/Test2.txt";
	moldata.readInSpectraFromFile(specfile);
	pspecfile = "tests/test_data/test_pspec/Test2.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i) );
	if( fabs( scores[0] - 0.0 ) > tol || fabs( scores[1] - 0.0 ) > tol || fabs( scores[2] - 0.0 ) > tol ){
		std::cout << "No Match: Expecting 0.0 scores but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	//Test between
	moldata = MolData("Test3","C");
	specfile = "tests/test_data/test_spec/Test3.txt";
	moldata.readInSpectraFromFile(specfile);
	pspecfile = "tests/test_data/test_pspec/Test3.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	double exp_scores[3] = {50.0, 57.142857142, 50.0};
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i) );
	if( fabs( scores[0] - exp_scores[0] ) > tol || fabs( scores[1] - exp_scores[1] ) > tol || fabs( scores[2] - exp_scores[2] ) > tol ){
		std::cout << "Mixed: Expecting " << exp_scores[0] << " " << exp_scores[1] << " " << exp_scores[2] << " but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	passed = pass;

}

ComparatorsTestJaccard::ComparatorsTestJaccard(){
	description = "Test computing of jaccard score";
}

void ComparatorsTestJaccard::runTest(){
	
	bool pass = true;
    double tol = 1e-5;

	Jaccard cmp(10.0, 0.01);
	std::vector<double> scores(3);

	//Test with exact matches
	MolData moldata("Test1","C");
	std::string specfile = "tests/test_data/test_spec/Test1.txt";
	moldata.readInSpectraFromFile(specfile);
	std::string pspecfile = "tests/test_data/test_pspec/Test1.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i) );
	if( fabs( scores[0] - 1.0 ) > tol || fabs( scores[1] - 1.0 ) > tol || fabs( scores[2] - 1.0 ) > tol ){
		std::cout << "Perfect Matches: Expecting 1.0 scores but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	//Test with no matches
	moldata = MolData("Test2","C");
	specfile = "tests/test_data/test_spec/Test2.txt";
	moldata.readInSpectraFromFile(specfile);
	pspecfile = "tests/test_data/test_pspec/Test2.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i) );
	if( fabs( scores[0] - 0.0 ) > tol || fabs( scores[1] - 0.0 ) > tol || fabs( scores[2] - 0.0 ) > tol ){
		std::cout << "No Match: Expecting 0.0 scores but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	//Test between
	moldata = MolData("Test3","C");
	specfile = "tests/test_data/test_spec/Test3.txt";
	moldata.readInSpectraFromFile(specfile);
	pspecfile = "tests/test_data/test_pspec/Test3.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	double exp_scores[3] = {0.666666667, 0.61538461538461542, 0.5};
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i) );
	if( fabs( scores[0] - exp_scores[0] ) > tol || fabs( scores[1] - exp_scores[1] ) > tol || fabs( scores[2] - exp_scores[2] ) > tol ){
		std::cout << "Mixed: Expecting " << exp_scores[0] << " " << exp_scores[1] << " " << exp_scores[2] << " but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	passed = pass;

}
