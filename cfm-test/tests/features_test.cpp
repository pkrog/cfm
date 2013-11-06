/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features_test.cpp
#
# Description: Test code for features.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/

#include "features_test.h"
#include "Util.h"
#include "MolData.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <boost/filesystem.hpp>

void initMolProps( romol_ptr_t &mol ){
	RDKit::ROMol::AtomIterator ai; 
	for( ai = mol.get()->beginAtoms(); ai != mol.get()->endAtoms(); ++ai ){	
		(*ai)->setProp("Root",0);
		(*ai)->setProp("OtherRoot",0);
	}
	mol.get()->setProp("IsRingBreak", 0);
}

FeaturesTestInit::FeaturesTestInit(){
	description = "Test initialisation of feature class";
}

void FeaturesTestInit::runTest(){
	
	bool pass = true;
    
	//Valid Config
	std::string config_filename = "tests/test_data/valid_feature_config.txt";
	FeatureCalculator *fc = new FeatureCalculator( config_filename );
	std::vector<std::string> fnames = fc->getFeatureNames();
	if( fnames.size() != 2 ){
		std::cout << "Unexpected number of feature names" << std::endl;
		pass = false;
	}
	else{
		if( fnames[0] != "BreakAtomPair" ||
			fnames[1] != "RingFeatures" ){
			std::cout << "Unexpected feature names" << std::endl;
			pass = false;		
		}
		if( fc->getNumFeatures() != 85 ){
			std::cout << "Unexpected feature count";
			pass = false;			
		}
	}
	delete fc;

	//Invalid Config
	bool except_thrown = false;
	try{
		std::string config_filename = "tests/test_data/invalid_feature_config.txt";
		FeatureCalculator *fc = new FeatureCalculator( config_filename );
	}
	catch( InvalidConfigException e ){
		except_thrown = true;
	}
	if( !except_thrown ){
		std::cout << "Allowed invalid feature configuration" << std::endl;
		pass = false;		
	}
	passed = pass;

}

FeaturesTestBreakAtomPair::FeaturesTestBreakAtomPair(){
	description = "Test BreakAtomPair feature";
}

void FeaturesTestBreakAtomPair::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("BreakAtomPair");
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	//Non-Ring Break C-N
	RDKit::Atom *null_atom = NULL;
	romol_ptr_t ion( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	initMolProps(ion);
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom);
	romol_ptr_t nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("N")) );
	initMolProps(nl);
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom);

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 3 ){
			std::cout << "Unexpected value for non-ring C-N root pair" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Non-Ring Break X-C
	ion = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("B")) );
	initMolProps(ion);
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), null_atom );
	nl = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	initMolProps(nl);
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), null_atom );

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 61 ){
			std::cout << "Unexpected value for non-ring X-C root pair" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Ring Break double C-N
	ion = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("CC")));
	initMolProps(ion);
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(1) );
	nl = romol_ptr_t(RDKit::SmilesToMol("NN"));
	initMolProps(nl);
	nl.get()->setProp("IsRingBreak", 1);
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(1) );

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 4 ){
			std::cout << "Unexpected value for ring double C-N root pair" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Ring Break C-N, X-X
	ion = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("CB")));
	initMolProps(ion);
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(1) );
	nl = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("NB")));
	initMolProps(nl);
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(1) );
	nl.get()->setProp("IsRingBreak", 1);

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 3 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 4 ){
			std::cout << "Unexpected value for ring C-N, X-X root pair" << std::endl;
			pass = false;			
		}
		if( fv->getFeature(2) != 72 ){
			std::cout << "Unexpected value for ring C-N, X-X root pair" << std::endl;
			pass = false;			
		}
	}
	delete fv;
	delete fc;
	passed = pass;
}

FeaturesTestRootPairs::FeaturesTestRootPairs(){
	description = "Test IonRootPairs and NLRootPairs feature";
}

void FeaturesTestRootPairs::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("IonRootPairs"); //NLRootPairs should work the same
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	//Non-ring "C-N,C-N,C-X"
	RDKit::Atom *null_atom = NULL;
	romol_ptr_t ion( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C(N)(B)N")) );
	romol_ptr_t nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	initMolProps(ion);
	initMolProps(nl);
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom );
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom );

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 4 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		bool ok = true;
		if( fv->getFeature(1) != 6 ) ok = false;
		if( fv->getFeature(2) != 7 ) ok = false;
		if( fv->getFeature(3) != 22 ) ok = false;
		if(!ok){
			std::cout << "Unexpected value for non-ring C-N,C-N,C-X" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Non-ring "X-X,X-N"
	ion = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("B(B)N")) );
	initMolProps(ion);
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), null_atom );
	
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 3 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		bool ok = true;
		if( fv->getFeature(1) != 126 ) ok = false;	
		if( fv->getFeature(2) != 142 ) ok = false;	
		if(!ok){
			std::cout << "Unexpected value for non-ring X-X,X-N" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Ring "C-N,C-N,X-X,C-X,X-N"
	ion = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C(N)(B)NNBB")) );
	initMolProps(ion);	
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(5) );
	nl.get()->setProp("IsRingBreak", 1);

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 6 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		bool ok = true;
		if( fv->getFeature(1) != 8 ) ok = false;
		if( fv->getFeature(2) != 9 ) ok = false;
		if( fv->getFeature(3) != 24 ) ok = false;
		if( fv->getFeature(4) != 128 ) ok = false;	
		if( fv->getFeature(5) != 144 ) ok = false;	
		if(!ok){
			std::cout << "Unexpected value for ring C-N,C-N,X-X,C-X,X-N" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Empty pairs (set feature indicating no pairs)
	ion = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	initMolProps(ion);	
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), null_atom );
	nl.get()->setProp("IsRingBreak", 0);

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 1 ){ 
			std::cout << "Missing feature flagging no pairs" << std::endl;
			pass = false;
		}
	}
	delete fv;

	delete fc;
	passed = pass;

}

FeaturesTestRootTriples::FeaturesTestRootTriples(){
	description = "Test IonRootTriples and NLRootTriples feature";
}

void FeaturesTestRootTriples::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("NLRootTriples"); //IonRootTriples should work the same
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	//Non-ring "C-C-N,C-C-N,C-X-C"
	romol_ptr_t ion( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	romol_ptr_t nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C(BC)C(N)N")) );
	initMolProps(ion);
	initMolProps(nl);
	RDKit::Atom *null_atom = NULL;
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom );
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom ); 

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 4 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		bool ok = true;
		if( fv->getFeature(1) != 6 ) ok = false;
		if( fv->getFeature(2) != 7 ) ok = false;
		if( fv->getFeature(3) != 122 ) ok = false;
		if(!ok){
			std::cout << "Unexpected value for non-ring C-C-N,C-C-N,C-X-C" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Non-ring "N-X-X"
	nl = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("NBB")) );
	initMolProps(nl);	
	rtd_nl = RootedROMolPtr(nl, nl.get()->getAtomWithIdx(0), null_atom );

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 286 ){
			std::cout << "Unexpected value for non-ring N-X-X" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Ring "C-C-N,C-C-N,C-X-C,N-X-X"
	ion = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	nl = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C(BC)C(N)NBBN")) );
	initMolProps(ion);
	initMolProps(nl);	
	rtd_nl = RootedROMolPtr(nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(8) );
	nl.get()->setProp("IsRingBreak", 1);

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 5 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		bool ok = true;
		if( fv->getFeature(1) != 8 ) ok = false;
		if( fv->getFeature(2) != 9 ) ok = false;
		if( fv->getFeature(3) != 124 ) ok = false;
		if( fv->getFeature(4) != 288 ) ok = false;	
		if(!ok){
			std::cout << "Unexpected value for non-ring C-C-N,C-C-N,C-X-C,N-X-X" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Empty triples (set feature indicating no pairs)
	nl = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	initMolProps(nl);
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), null_atom );
	
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 1 ){ 
			std::cout << "Missing feature flagging no triples" << std::endl;
			pass = false;
		}
	}
	delete fv;
	delete fc;
	passed = pass;
}

FeaturesTestGasteigerCharges::FeaturesTestGasteigerCharges(){
	description = "Test Gasteiger Charges feature";
}

void FeaturesTestGasteigerCharges::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("GasteigerCharges"); //IonRootTriples should work the same
	FeatureCalculator *fc = new FeatureCalculator( fnames );
	RDKit::Atom *null_atom = NULL;

	//Non-ring (-0.33042488,-0.00652530)
	romol_ptr_t ion( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	initMolProps(ion);
	ion.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge", -0.33042488 );
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom );
	romol_ptr_t nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	initMolProps(nl);	
	nl.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge", -0.00652530 );
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom );
	nl.get()->setProp("IsRingBreak",0);

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 9 ){		
			std::cout << "Unexpected value for ion gasteiger charge (non-ring) " << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Ring Retain Order (-0.01,0.9), (0.05,-0.4)
	ion = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("CC")) );
	initMolProps(ion);	
	ion.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge", -0.01 );
	ion.get()->getAtomWithIdx(1)->setProp<double>("OrigGasteigerCharge", 0.05 );
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(1) );	
	nl = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("CC")) );
	initMolProps(nl);
	nl.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge", 0.9 );
	nl.get()->getAtomWithIdx(1)->setProp<double>("OrigGasteigerCharge", -0.4 );
	nl.get()->setProp("IsRingBreak",1);	
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(1) );
	

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 3 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 54 ){		
			std::cout << "Unexpected value for ion gasteiger charge ring " << fv->getFeature(1) << std::endl;
			pass = false;
		}
		if( fv->getFeature(2) != 56 ){		
			std::cout << "Unexpected value for nl gasteiger charge ring " << fv->getFeature(2) << std::endl;
			pass = false;
		}
	}
	delete fv;
	delete fc;
	passed = pass;

}

FeaturesTestHydrogenMovement::FeaturesTestHydrogenMovement(){
	description = "Test Hydrogen Movement feature";
}

void FeaturesTestHydrogenMovement::runTest(){
	
	bool pass = true;
	RDKit::Atom *null_atom = NULL;
	std::vector<std::string> fnames;
	fnames.push_back("HydrogenMovement"); //IonRootTriples should work the same
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	//Positive within range
	romol_ptr_t ion( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	initMolProps(ion);
	double h_movement = 3.00452;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom );
	romol_ptr_t nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	initMolProps(nl);
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom );

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 8 ){		
			std::cout << "Unexpected idx for hydrogen movement" << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Negative within range
	h_movement = -1.115;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 4 ){		
			std::cout << "Unexpected idx for hydrogen movement" << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Zero
	h_movement = 0.0;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 5 ){		
			std::cout << "Unexpected idx for hydrogen movement" << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Out of range positive
	h_movement = 5.0;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 10 ){		
			std::cout << "Unexpected idx for hydrogen movement" << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Out of range negative
	h_movement = -5.0;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 10 ){		
			std::cout << "Unexpected idx for hydrogen movement" << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;
	delete fc;
	passed = pass;

}

FeaturesTestRingFeatures::FeaturesTestRingFeatures(){
	description = "Test Ring features";
}

void FeaturesTestRingFeatures::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("RingFeatures");
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	//1,1,3,6
	romol_ptr_t ion( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("CCC")) );
	romol_ptr_t nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("CCC")) );
	initMolProps(ion);
	initMolProps(nl);
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(2) );
	nl.get()->setProp("IsRingBreak",1);
	nl.get()->setProp("IsAromaticRingBreak",1);
	nl.get()->setProp("IsAromaticDblRingBreak",1);
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(2) );

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 5 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 2 || fv->getFeature(2) != 3 || fv->getFeature(3) != 6 || fv->getFeature(4) != 11 ){		
			std::cout << "Unexpected features for ring break 1,1,3,6 " << std::endl;
			pass = false;
		}
	}
	delete fv;

	//0,0,4,9
	ion = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("CCCC")) );
	nl = romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("CCCCC")) );
	initMolProps(ion);
	initMolProps(nl);
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(3) );
	nl.get()->setProp("IsRingBreak",1);
	nl.get()->setProp("IsAromaticRingBreak",0);
	nl.get()->setProp("IsAromaticDblRingBreak",0);
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(4) );

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 4 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 1 || fv->getFeature(2) != 7 || fv->getFeature(3) != 12 ){		
			std::cout << "Unexpected features for ring break 0,0,4,9 " << std::endl;
			pass = false;
		}
	}
	delete fv;
	delete fc;
	passed = pass;

}


FeaturesTestQuadraticFeatures::FeaturesTestQuadraticFeatures(){
	description = "Test quadratic features";
}

void FeaturesTestQuadraticFeatures::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("BreakAtomPair"); 
	fnames.push_back("HydrogenMovement"); 
	fnames.push_back("QuadraticFeatures"); 
	FeatureCalculator *fc = new FeatureCalculator( fnames );	

	//Simple initial vector with 3 bits set (indexes: 0,3,80 )
	romol_ptr_t ion( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	initMolProps(ion);	
	RDKit::Atom *null_atom = NULL;
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom );
	double h_movement = 3.00452;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	romol_ptr_t nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("N")) );
	initMolProps(nl);	
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom );
	nl.get()->setProp("IsRingBreak",0);

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 4 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;	
	}
	else{
		if( fv->getFeature(0) != 0 || fv->getFeature(1) != 3 || fv->getFeature(2) != 80 ){
			std::cout << "Unexpected singular features" << std::endl;
			pass = false;				
		}
		if( fv->getFeature(3) != 3166 ){
			std::cout << "Unexpected quadratic feature:" << fv->getFeature(3) << std::endl;
			pass = false;				
		}
	}
	int total = fv->getTotalLength();
	if( total != 3404 ){
		std::cout << "Unexpected total feature count: " << total << std::endl;
		pass = false;		
	}
	int numfeatures = fc->getNumFeatures();
	if( numfeatures != 3404 ){
		std::cout << "Unexpected total feature calculation: " << numfeatures << std::endl;
		pass = false;		
	}
	delete fv;
	delete fc;
	passed = pass;
}


FeaturesTestLength::FeaturesTestLength(){
	description = "Test length of features";
}

void FeaturesTestLength::runTest(){

	bool pass = true;
	
	//Create the feature calculator
	const std::vector<std::string> names = FeatureCalculator::getValidFeatureNames();

	//Create some aribitrary input data
	romol_ptr_t ion( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
	initMolProps(ion);		
	RDKit::Atom *null_atom = NULL;
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom );
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 15.0);
	ion.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge", 1.0);
	romol_ptr_t nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("N")) );
	initMolProps(nl);
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom );
	nl.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge",-1.0);
	nl.get()->setProp("IsRingBreak",0);

	//Check all feature lengths
	std::vector<std::string>::const_iterator it = names.begin();
	for( ; it != names.end(); ++it ){
		
		std::vector<std::string> feature_list;
		feature_list.push_back( *it );
		FeatureCalculator fc( feature_list );

		//Compute the feature vector
		FeatureVector *fv = fc.computeFV( &rtd_ion, &rtd_nl );

		//Check the length
		if( fv->getTotalLength() != fc.getNumFeatures() ){ 
			std::cout << "Feature length incorrect for " << *it;
			std::cout << ": expecting " << fc.getNumFeatures();
			std::cout << " but found " << fv->getTotalLength() << std::endl; 
			pass = false;
		}
		delete fv;
	}
	passed = pass;
}

FeaturesTestMetlinExample::FeaturesTestMetlinExample(){
	description = "Test Metlin Example Features";
}

void FeaturesTestMetlinExample::runTest(){

	bool pass = true;
	std::string smiles_Metlin_21361 = "O=C(NC(CCC(N)=O)C(O)=O)C(C)NC(=O)C(N)CC(C)C";

	std::string config_filename = "tests/test_data/example_feature_config.txt";
	FeatureCalculator fc( config_filename );

	//Ingegration Test - compute the fragment tree and transitions
	std::string id = "Metlin_21361";
	MolData mol( id, smiles_Metlin_21361, 0 );
	mol.computeFragmentGraph(1);
	mol.computeFeatureVectors( &fc );

	//Transition 0: Check that the features are as expected
	const FeatureVector *fv = mol.getFeatureVectorForIdx(0);
	if( fv->getNumSetFeatures() != 11 ){
		std::cout << "Unexpected number of non-zero features " << fv->getNumSetFeatures();
		pass = false;
	}
	else{
		if( fv->getFeature(0) != 0 ){ 
			std::cout << "Bias error: " << fv->getFeature(0) <<std::endl;
			pass = false;
		}
		if( fv->getFeature(1) != 5 ){ 
			std::cout << "Incorrect BreakAtomPair: " << fv->getFeature(1) << std::endl;
			pass = false;
		}
		if( fv->getFeature(2) != 98 ){
			std::cout << "Incorrect Gasteiger: " << fv->getFeature(2) << std::endl;
			pass = false;
		}
		if( fv->getFeature(3) != 148 ){
			std::cout << "Incorrect HMovement: " << fv->getFeature(3)<< std::endl;
			pass = false;
		}
		if( fv->getFeature(4) != 156){
			std::cout << "Incorrect IonRootPair C-C: " << fv->getFeature(4) << std::endl;
			pass = false;
		}
		if( fv->getFeature(5) != 160 ){
			std::cout << "Incorrect IonRootPair C-N: " << fv->getFeature(5) << std::endl;
			pass = false;
		}
		if( fv->getFeature(6) != 301 ){
			std::cout << "Incorrect IonRootTriple C-C-C: " << fv->getFeature(6) << std::endl;
			pass = false;
		}
		if( fv->getFeature(7) != 305 ){
			std::cout << "Incorrect IonRootTriple C-C-N: " << fv->getFeature(7) << std::endl;
			pass = false;
		}
		if( fv->getFeature(8) != 325 ){
			std::cout << "Incorrect IonRootTriple C-N-C: " << fv->getFeature(8)<< std::endl;
			pass = false;
		}
		if( fv->getFeature(9) != 1165 ){
			std::cout << "Incorrect NLRootPair (None): " << fv->getFeature(9) << std::endl;
			pass = false;
		}
		if( fv->getFeature(10) != 1310 ){
			std::cout << "Incorrect NLRootTriple (None): " << fv->getFeature(10) << std::endl;
			pass = false;
		}
	}

	passed = pass;
}
