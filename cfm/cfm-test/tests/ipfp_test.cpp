/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# ipfp_test.cpp
#
# Description: Test code for IPFP.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/

#include "ipfp_test.h"
#include "Config.h"
#include "FragmentTreeNode.h"

#include <boost/filesystem.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>

class IPFPTestMolConverge : public MolData {
public:
	IPFPTestMolConverge() : MolData("IPFP Test Mol Converge Case", ""){
	
		setIonizationMode(false);

		//Create a molecule based on what was previously in test_bn_transition_ipfp.txt
		(*fg) = new FragmentGraph();
		romol_ptr_t basic_nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("O=C(O)CNC(=O)C(NC(=O)CN)Cc1ccccc1")) ), basic_nl, -1, -1, false), -1 ); //id = 0
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("NCC(=O)[NH2+]C(=C=O)Cc1ccccc1")) ), basic_nl, -1, -1, false), 0 ); // id = 1, 0 -> 1
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("N=CC(=O)[NH2+]CCc1ccccc1")) ), basic_nl, -1, -1, false), 0 ); //id = 2, 0 -> 2
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("[NH2+]=CCc1ccccc1")) ), basic_nl, -1, -1, false), 2 ); // id = 3, 2 -> 3

		//Set thetas/transition probs to match matlab reference2
		(*thetas).resize( 3 );
		for( int energy = 0; energy < 3; energy++ ){
			(*thetas)[energy].resize( (*fg)->getNumTransitions() );
			(*thetas)[energy][0] = 1.386294361119891;	//0->1
			(*thetas)[energy][1] = 2.456735772821304;	//0->2
			(*thetas)[energy][2] = 0.0;	//2->3
		}
		computeTransitionProbabilities();

		//Set the spectra
		spectra.resize(3);
		spectra[0].push_back( Peak(120.082039, 35.914500 ) );
		spectra[0].push_back( Peak(177.103613, 50.000000 ) );
		spectra[0].push_back( Peak(280.133689, 50.146650) );
		spectra[1].push_back( Peak(120.083820, 100.00000 ) );
		spectra[1].push_back( Peak(177.106284, 33.000500) );
		spectra[2].push_back( Peak(120.081802, 100.00000) );
		for( int i = 0; i <=2; i++ ) 
			MolData::normalizeAndSortSpectrum(spectra[i]);

	}
};

class IPFPTestMolNonConverge : public MolData {
public:
	IPFPTestMolNonConverge() : MolData("IPFP Test Mol Non-Converge Case", ""){

		setIonizationMode(false);

		//Create a molecule based on what was previously in test_bn_transition_ipfp_nonconverge.txt
		(*fg) = new FragmentGraph();
		romol_ptr_t basic_nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("O=C(O)CNC(=O)C(NC(=O)CN)Cc1ccccc1")) ), basic_nl, -1, -1, false), -1 ); //id = 0
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("NCC(=O)[NH2+]C(=C=O)Cc1ccccc1")) ), basic_nl, -1, -1, false), 0 ); // id = 1, 0 -> 1
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("N=CC(=O)[NH2+]CCc1ccccc1")) ), basic_nl, -1, -1, false), 0 ); //id = 2, 0 -> 2
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("[NH2+]=CCc1ccccc1")) ), basic_nl, -1, -1, false), 2 ); // id = 3, 2 -> 3

		//Set thetas/transition probs to match matlab reference2
		(*thetas).resize( 3 );
		for( int energy = 0; energy < 3; energy++ ){
			(*thetas)[energy].resize( (*fg)->getNumTransitions() );
			(*thetas)[energy][0] = 1.386294361119891;	//0->1
			(*thetas)[energy][1] = 2.456735772821304;	//0->2
			(*thetas)[energy][2] = 0.0;	//2->3
		}
		computeTransitionProbabilities();

		//Set the spectra
		spectra.resize(3);
		spectra[0].push_back( Peak(120.082039, 35.914500 ) );
		spectra[0].push_back( Peak(177.103613, 50.000000 ) );
		spectra[0].push_back( Peak(280.133689, 50.146650) );
		spectra[0].push_back( Peak(205.098456, 35.223700) );
		spectra[1].push_back( Peak(120.083820, 100.00000 ) );
		spectra[1].push_back( Peak(177.106284, 33.000500) );
		spectra[2].push_back( Peak(120.081802, 100.00000) );
		spectra[2].push_back( Peak(177.106284, 100.00000) );
		for( int i = 0; i <=2; i++ ) 
			MolData::normalizeAndSortSpectrum(spectra[i]);
	}
};

class IPFPTestMolSharedMass : public MolData {
public:
	IPFPTestMolSharedMass() : MolData("IPFP Test Mol Shared Mass Case", ""){
	
		setIonizationMode(false);

		//Create a molecule based on what was previously in test_bn_transition_ipfp_sharedmass.txt
		(*fg) = new FragmentGraph();
		romol_ptr_t basic_nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("O=C(O)CNC(=O)C(NC(=O)CN)Cc1ccccc1")) ), basic_nl,-1, -1, false), -1 ); //id = 0
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("NCC(=O)[NH2+]C(=C=O)Cc1ccccc1")) ), basic_nl,-1, -1, false), 0 ); // id = 1, 0 -> 1
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("N=CC(=O)[NH2+]CCc1ccccc1")) ), basic_nl,-1, -1, false), 0 ); //id = 2, 0 -> 2
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("[NH2+]=CCc1ccccc1")) ), basic_nl,-1, -1, false), 2 ); // id = 3, 2 -> 3
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("N=CC[NH2+]C(=O)Cc1ccccc1")) ), basic_nl,-1, -1, false), 0 ); //id = 4, 0 -> 4
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("[NH2+]=CCc1ccccc1")) ), basic_nl,-1, -1, false), 4 ); //  4 -> 3		

		//Set thetas/transition probs to match matlab reference2
		(*thetas).resize( 3 );
		for( int energy = 0; energy < 3; energy++ ){
			(*thetas)[energy].resize( (*fg)->getNumTransitions() );
			(*thetas)[energy][0] = 1.386294361119891;	//0->1
			(*thetas)[energy][1] = 2.456735772821304;	//0->2
			(*thetas)[energy][2] = 0.0;	//2->3
			(*thetas)[energy][3] = 5.0;	//0->4
		}
		computeTransitionProbabilities();

		//Set the spectra
		spectra.resize(3);
		spectra[0].push_back( Peak(120.082039, 35.914500 ) );
		spectra[0].push_back( Peak(177.103613, 50.000000 ) );
		spectra[0].push_back( Peak(280.133689, 50.146650) );
		spectra[1].push_back( Peak(120.083820, 100.00000 ) );
		spectra[1].push_back( Peak(177.106284, 33.000500) );
		spectra[2].push_back( Peak(120.081802, 100.00000) );
		for( int i = 0; i <=2; i++ ) 
			MolData::normalizeAndSortSpectrum(spectra[i]);
	}
};

class IPFPTestMolInterpolate : public MolData {
public:
	IPFPTestMolInterpolate() : MolData("IPFP Test Mol Interpolation Case", ""){
	
		setIonizationMode(false);

		//Create a molecule based on what was previously in test_bn_transition_ipfp_interp.txt
		(*fg) = new FragmentGraph();
		romol_ptr_t basic_nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("O=C(O)CNC(=O)C(NC(=O)CN)Cc1ccccc1")) ), basic_nl,-1, -1, false), -1 ); //id = 0
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("NCC(=O)[NH2+]C(=C=O)Cc1ccccc1")) ), basic_nl,-1, -1, false), 0 ); // id = 1, 0 -> 1
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("N=CC(=O)[NH2+]CCc1ccccc1")) ), basic_nl,-1, -1, false), 0 ); //id = 2, 0 -> 2
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("N=CC(=O)[NH2+]CCc1ccccc1")) ), basic_nl,-1, -1, false), 1 ); //1 -> 2
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("[NH2+]=CCc1ccccc1")) ), basic_nl, -1, -1, false),2 ); // id = 3, 2 -> 3	

		//Set thetas/transition probs to match matlab reference2
		(*thetas).resize( 3 );
		for( int energy = 0; energy < 3; energy++ ){
			(*thetas)[energy].resize( (*fg)->getNumTransitions() );
			(*thetas)[energy][0] = 1.386294361119891;	//0->1
			(*thetas)[energy][1] = 2.456735772821304;	//0->2
			(*thetas)[energy][2] = 5.0;	//1->2
			(*thetas)[energy][3] = 0.0;	//2->3
		}
		computeTransitionProbabilities();

		//Set the spectra
		spectra.resize(3);
		spectra[0].push_back( Peak(120.082039, 35.914500 ) );
		spectra[0].push_back( Peak(177.103613, 50.000000 ) );
		spectra[0].push_back( Peak(280.133689, 50.146650) );
		spectra[1].push_back( Peak(120.083820, 100.00000 ) );
		spectra[1].push_back( Peak(177.106284, 33.000500) );
		spectra[2].push_back( Peak(120.081802, 100.00000) );
		for( int i = 0; i <=2; i++ ) 
			MolData::normalizeAndSortSpectrum(spectra[i]);
	}
};

class IPFPTestMolInterpolateNoDirect : public MolData {
public:
	IPFPTestMolInterpolateNoDirect() : MolData("IPFP Test Mol Interpolation with No Direct Case", ""){
	
		setIonizationMode(false);

		//Create a molecule based on what was previously in test_bn_transition_ipfp_interp.txt
		(*fg) = new FragmentGraph();
		romol_ptr_t basic_nl( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("C")) );
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("O=C(O)CNC(=O)C(NC(=O)CN)Cc1ccccc1")) ), basic_nl,-1, -1, false), -1 ); //id = 0
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("NCC(=O)[NH2+]C(=C=O)Cc1ccccc1")) ), basic_nl,-1, -1, false), 0 ); // id = 1, 0 -> 1
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("N=CC(=O)[NH2+]CCc1ccccc1")) ), basic_nl, -1, -1, false),1 ); //id = 2, 1 -> 2
		(*fg)->addToGraph( FragmentTreeNode( romol_ptr_t( static_cast<RDKit::ROMol *>(RDKit::SmilesToMol("[NH2+]=CCc1ccccc1")) ), basic_nl, -1, -1, false),2 ); // id = 3, 2 -> 3	

		//Set thetas/transition probs to match matlab reference2
		(*thetas).resize( 3 );
		for( int energy = 0; energy < 3; energy++ ){
			(*thetas)[energy].resize( (*fg)->getNumTransitions() );
			(*thetas)[energy][0] = 1.386294361119891;	//0->1
			(*thetas)[energy][1] = 5.0;	//1->2
			(*thetas)[energy][2] = 0.0;	//2->3
		}
		computeTransitionProbabilities();

		//Set the spectra
		spectra.resize(3);
		spectra[0].push_back( Peak(120.082039, 35.914500 ) );
		spectra[0].push_back( Peak(177.103613, 50.000000 ) );
		spectra[0].push_back( Peak(280.133689, 50.146650) );
		spectra[1].push_back( Peak(120.083820, 100.00000 ) );
		spectra[1].push_back( Peak(177.106284, 33.000500) );
		spectra[2].push_back( Peak(120.081802, 100.00000) );
		for( int i = 0; i <=2; i++ ) 
			MolData::normalizeAndSortSpectrum(spectra[i]);
	}
};


void printBeliefs( std::vector<double> &beliefs ){
	for( unsigned int i = 0; i < beliefs.size(); i++ ){
		std::cout << exp(beliefs[i]) << std::endl;
	}
}

IPFPTestComputeBeliefsConverge::IPFPTestComputeBeliefsConverge(){
	description = "Test computing of ipfp beliefs converging case";
}

void IPFPTestComputeBeliefsConverge::runTest(){
	
	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 6;
	int tmp_array[3] = {2,4,6};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	double tmp_array2[3] = {0.33,0.33,0.34};
	cfg.spectrum_weights.assign(tmp_array2, tmp_array2+3);
	cfg.interpolate_spectra = 0;
	initDerivedConfig(cfg);
	double tol = cfg.ipfp_converge_thresh;

	//Load some molecule data with transition probabilities set as in matlab reference2
	IPFPTestMolConverge moldata;

	//Run code
	for( int algorithm_to_use = 0; algorithm_to_use <= 2; algorithm_to_use++ ){

		cfg.ipfp_algorithm = algorithm_to_use;
		if( algorithm_to_use == IPFP_ALGORITHM ) std::cout << "with IPFP...";
		else if( algorithm_to_use == GEMA_ALGORITHM ) std::cout << "with GEMA...";
		else std::cout << "with IPFP with oscillatory adjustment...";
		
		IPFP *ipfp = new IPFP( &moldata, &cfg ); 
		beliefs_t *beliefs = ipfp->calculateBeliefs();

		if( ipfp->status != COMPLETE_CONVERGE ){
			std::cout << "IPFP did not converge" << std::endl;
			pass = false;
		}

		//Check results
		//0->0:
		if( fabs(exp(beliefs->ps[0][0]) - 0.407932704207945) > tol ||
			fabs(exp(beliefs->ps[0][1]) - 0.368559651303844) > tol ||
			fabs(exp(beliefs->ps[0][2]) - 0.018441085137235) > tol ||
			fabs(exp(beliefs->ps[0][3]) - 0.0) > tol ||
			fabs(exp(beliefs->ps[0][4]) - 0.0) > tol ||
			fabs(exp(beliefs->ps[0][5]) - 0.0) > tol ){
			std::cout << "Unexpected value for 0->0 persistence belief" << std::endl;
			printBeliefs( beliefs->ps[0] );
			pass = false;
		}

		//0->1:
		if( fabs(exp(beliefs->tn[0][0]) - 0.0) > tol ||
			fabs(exp(beliefs->tn[0][1]) - 0.0) > tol ||
			fabs(exp(beliefs->tn[0][2]) - 0.0) > tol ||
			fabs(exp(beliefs->tn[0][3]) - 0.0) > tol ||
			fabs(exp(beliefs->tn[0][4]) - 0.0) > tol ||
			fabs(exp(beliefs->tn[0][5]) - 0.0) > tol ){
			std::cout << "Unexpected value for 0->1 transition belief" << std::endl;
			printBeliefs( beliefs->tn[0] );		
			pass = false;
		}

		//0->2:
		if( fabs(exp(beliefs->tn[1][0]) - 0.592067295792055) > tol ||
			fabs(exp(beliefs->tn[1][1]) - 0.039373052904101) > tol ||
			fabs(exp(beliefs->tn[1][2]) - 0.350118566166608) > tol ||
			fabs(exp(beliefs->tn[1][3]) - 0.018441085137235) > tol ||
			fabs(exp(beliefs->tn[1][4]) - 0.0) > tol ||
			fabs(exp(beliefs->tn[1][5]) - 0.0) > tol ){
			std::cout << "Unexpected value for 0->2 transition belief" << std::endl;
			printBeliefs( beliefs->tn[1] );
			pass = false;
		}

		//2->2:
		if( fabs(exp(beliefs->ps[2][0]) - 0.0) > tol ||
			fabs(exp(beliefs->ps[2][1]) - 0.328108774200844) > tol ||
			fabs(exp(beliefs->ps[2][2]) - 0.173164830863293) > tol ||
			fabs(exp(beliefs->ps[2][3]) - 0.229682042219429) > tol ||
			fabs(exp(beliefs->ps[2][4]) - 0.082707709118888) > tol ||
			fabs(exp(beliefs->ps[2][5]) - 0.0) > tol ){
			std::cout << "Unexpected value for 2->2 persistence belief" << std::endl;
			printBeliefs( beliefs->ps[2] );
			pass = false;
		}

		//2->3:
		if( fabs(exp(beliefs->tn[2][0]) - 0.0) > tol ||
			fabs(exp(beliefs->tn[2][1]) - 0.263958521591211) > tol ||
			fabs(exp(beliefs->tn[2][2]) - 0.194316996241652) > tol ||
			fabs(exp(beliefs->tn[2][3]) - 0.293601354810473) > tol ||
			fabs(exp(beliefs->tn[2][4]) - 0.165415418237776) > tol ||
			fabs(exp(beliefs->tn[2][5]) - 0.082707709118888) > tol ){
			std::cout << "Unexpected value for 2->3 transition belief" << std::endl;
			printBeliefs( beliefs->tn[2] );
			pass = false;
		}

		//3->3:
		if( fabs(exp(beliefs->ps[3][0]) - 0.0) > tol ||
			fabs(exp(beliefs->ps[3][1]) - 0.0) > tol ||
			fabs(exp(beliefs->ps[3][2]) - 0.263958521591211) > tol ||
			fabs(exp(beliefs->ps[3][3]) - 0.458275517832863) > tol ||
			fabs(exp(beliefs->ps[3][4]) - 0.751876872643336) > tol ||
			fabs(exp(beliefs->ps[3][5]) - 0.917292290881112) > tol ){
			std::cout << "Unexpected value for 3->3 persistence belief" << std::endl;
			printBeliefs( beliefs->ps[3] );
			pass = false;
		}
		delete ipfp;
	}
	passed = pass;

}

IPFPTestComputeBeliefsNonConverge::IPFPTestComputeBeliefsNonConverge(){
	description = "Test computing of ipfp beliefs non-converging case";
}

void IPFPTestComputeBeliefsNonConverge::runTest(){

	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 6;
	int tmp_array[3] = {2,4,6};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	double tmp_array2[3] = {1.0,1.0,1.0}; 
	cfg.spectrum_weights.assign(tmp_array2, tmp_array2+3);
	cfg.interpolate_spectra = 0;
	initDerivedConfig(cfg);

	//Load some test molecule data
	IPFPTestMolNonConverge moldata;

	//Run code for IPFP
	cfg.ipfp_algorithm = IPFP_ALGORITHM;
	IPFP *ipfp = new IPFP( &moldata, &cfg); 
	beliefs_t *beliefs = ipfp->calculateBeliefs();
    
	if( ipfp->status != NON_CONVERGE ){
		std::cout << "IPFP: Incorrect status for non-convergent case - expecting NON_CONVERGE" << std::endl;
		pass = false;
	}
	delete ipfp;
		
	//Run code for GEMA
	cfg.ipfp_algorithm = GEMA_ALGORITHM;
	//ipfp = new IPFP( &moldata, &cfg ); 
	//   beliefs = ipfp->calculateBeliefs();
	//   
	//if( ipfp->status != COMPLETE_CONVERGE ){
	//	std::cout << "GEMA: Incorrect status for non-convergent case - expecting COMPLETE_CONVERGE" << std::endl;
	//	pass = false;
	//}
	//delete ipfp;

	//Run code for IPFP with mod
	cfg.ipfp_algorithm = IPFP_WITH_MOD_ALGORITHM;
	ipfp = new IPFP( &moldata, &cfg ); 
	beliefs = ipfp->calculateBeliefs();
    
	if( ipfp->status != CONVERGE_AFTER_MOD ){
		std::cout << "IPFP with mod: Incorrect status for non-convergent case - expecting CONVERGE_AFTER_MOD" << std::endl;
		pass = false;
	}
	delete ipfp;

	passed = pass;

}


IPFPTestComputeBeliefsSharedMass::IPFPTestComputeBeliefsSharedMass(){
	description = "Test computing of ipfp beliefs, case with shared mass";
}

void IPFPTestComputeBeliefsSharedMass::runTest(){

	bool pass = true;
	double tol = 1e-6;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 6;
	int tmp_array[3] = {2,4,6};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	double tmp_array2[3] = {0.33,0.33,0.34};
	cfg.spectrum_weights.assign(tmp_array2, tmp_array2+3);
	cfg.interpolate_spectra = 0;
	initDerivedConfig(cfg);

	//Load some molecule data
	IPFPTestMolSharedMass moldata;

	//Run code
	cfg.ipfp_algorithm = IPFP_ALGORITHM;
	IPFP *ipfp = new IPFP( &moldata, &cfg ); 
	beliefs_t *beliefs = ipfp->calculateBeliefs();
    
	if( ipfp->status != COMPLETE_CONVERGE ){
		std::cout << "Incorrect status - expecting COMPLETE_CONVERGE" << std::endl;
		pass = false;
	}

	//Check the resulting marginals for fragments 2 and 4 (probability mass should be split
	//proportionally between them according to the theta values set above)
	//std::vector<message_t> *marginals = ipfp->getMostRecentMarginals();
	//double marg2 = exp( (*marginals)[0].message[2] - (*marginals)[0].log_sum );
	//double marg4 = exp( (*marginals)[0].message[4] - (*marginals)[0].log_sum );
	//if( fabs(marg2/marg4 - 6) > tol ){
	//	std::cout << "Expecting ratio of 6 between fragment 2 and 4 but found " << marg2/marg4 << std::endl;
	//	pass = false;	
	//}
	//if( fabs(marg2 + marg4 -  0.3674818) > tol ){
	//	std::cout << "Expecting sum of marginals 2 and 4 to be 0.3674818 but found " << marg2 + marg4 << std::endl;
	//	pass = false;		
	//}
	passed = pass;

}

IPFPTestComputeBeliefsInterpolate::IPFPTestComputeBeliefsInterpolate(){
	description = "Test computing of ipfp beliefs with interpolation";
}

void IPFPTestComputeBeliefsInterpolate::runTest(){

	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 6;
	int tmp_array[3] = {2,4,6};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	double tmp_array2[3] = {0.33,0.33,0.34};
	cfg.spectrum_weights.assign(tmp_array2, tmp_array2+3);
	cfg.interpolate_spectra = 1;
	cfg.intermediate_weights = 0.05;
	initDerivedConfig(cfg);

	//Load some molecule data and create the interpolated spectra
	IPFPTestMolInterpolate moldata;
	moldata.createInterpolatedSpectra(cfg);
		
	//Run code for GEMA
	cfg.ipfp_algorithm = GEMA_ALGORITHM;
	IPFP *ipfp = new IPFP( &moldata, &cfg ); 
    beliefs_t *beliefs = ipfp->calculateBeliefs();
    
	if( ipfp->status != COMPLETE_CONVERGE ){
		std::cout << "GEMA: Incorrect status for interpolated case - expecting COMPLETE_CONVERGE" << std::endl;
		pass = false;
	}
	
	//If the interpolation is working correctly, 0->1 and 1->2 should both be
	//very unlikely.
	double prob_thresh = 0.001;
	std::vector<double>::iterator it = beliefs->tn[0].begin();
	for( ; it != beliefs->tn[0].end(); ++it ){
		if( exp(*it) > prob_thresh ){
			std::cout << "0->1 probability over limit in interpolated case: " << exp(*it) << std::endl;
			pass = false;
		}
	}
	it = beliefs->tn[2].begin();
	for( ; it != beliefs->tn[2].end(); ++it ){
		if( exp(*it) > prob_thresh ){
			std::cout << "1->2 probability over limit in interpolated case: " << exp(*it) << std::endl;
			pass = false;
		}
	}

	delete ipfp;
	passed = pass;

}

IPFPTestComputeBeliefsInterpolateNoDirect::IPFPTestComputeBeliefsInterpolateNoDirect(){
	description = "Test computing of ipfp beliefs with interpolation for the case of no-direct path";
}

void IPFPTestComputeBeliefsInterpolateNoDirect::runTest(){

	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 6;
	int tmp_array[3] = {2,4,6};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	double tmp_array2[3] = {1.0,1.0,1.0};
	cfg.spectrum_weights.assign(tmp_array2, tmp_array2+3);
	cfg.interpolate_spectra = 1;
	cfg.intermediate_weights = 0.2;
	initDerivedConfig(cfg);

	//Load some molecule data
	IPFPTestMolInterpolateNoDirect moldata;
	moldata.createInterpolatedSpectra(cfg);
		
	//Run code for GEMA
	cfg.ipfp_algorithm = GEMA_ALGORITHM;
	IPFP *ipfp = new IPFP( &moldata, &cfg ); 
    beliefs_t *beliefs = ipfp->calculateBeliefs();
    
	if( ipfp->status != COMPLETE_CONVERGE ){
		std::cout << "GEMA: Incorrect status for interpolated case - expecting COMPLETE_CONVERGE" << std::endl;
		pass = false;
	}
	
	//Since there's no other path to 2 and 3, 0->1 and 1->2 should both be likely
	double prob_thresh = 0.01;
	if( exp(beliefs->tn[0][0]) < prob_thresh || exp(beliefs->tn[0][2]) < prob_thresh ){
		std::cout << "0->1 probability under limit in non-direct case: " << std::endl;
		printBeliefs( beliefs->tn[0] );
		pass = false;
	}

	if( exp(beliefs->tn[1][1]) < prob_thresh || exp(beliefs->tn[1][3]) < prob_thresh ){
		std::cout << "1->2 probability under limit in non-direct case: " << std::endl;
		printBeliefs( beliefs->tn[1] );
		pass = false;
	}

	delete ipfp;
	passed = pass;

}

IPFPTestComputeBeliefsSingleEnergy::IPFPTestComputeBeliefsSingleEnergy(){
	description = "Test computing of ipfp beliefs for single energy case";
}

void IPFPTestComputeBeliefsSingleEnergy::runTest(){

	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 6;
	int tmp_array[3] = {2,2,2};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	double tmp_array2[3] = {1.0,1.0,1.0};
	cfg.spectrum_weights.assign(tmp_array2, tmp_array2+3);
	initDerivedConfig(cfg);

	config_t se_cfg;
	initSingleEnergyConfig(se_cfg, cfg, 1);

	//Load some molecule data
	IPFPTestMolConverge moldata;
		
	//Run code for IPFP_OSC
	se_cfg.ipfp_algorithm = IPFP_WITH_MOD_ALGORITHM;
	IPFP *ipfp = new IPFP( &moldata, &se_cfg ); 
    beliefs_t *beliefs = ipfp->calculateBeliefs();
    
	if( ipfp->status != COMPLETE_CONVERGE ){
		std::cout << "IPFP_WITH_MOD_ALGORITHM: Incorrect status for single energy case - expecting COMPLETE_CONVERGE" << std::endl;
		pass = false;
	}

	passed = pass;
}