/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Isotope.cpp
#
# Description:  Class for generation of isotope peaks for a molecular formula.
#				Provides wrapper to emass from:
#		
# A. Rockwood and P. Haimi, "Efficient calculation of accurate masses of isotopic peaks.", 
# Journal of the American Society for Mass Spectrometry, 17:3 p415-9 2006.				
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Isotope.h"
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/MolOps.h>
#include <fstream>

typedef unsigned long ulong;

static const double DUMMY_MASS = -10000000;

#define ISOTOPES \
"X  2" "\n" \
"1  0.9" "\n" \
"2  0.1" "\n" \
"" "\n" \
"H  2" "\n" \
"1.0078246  0.99985" "\n" \
"2.0141021  0.00015" "\n" \
"" "\n" \
"He  2" "\n" \
"3.01603    0.00000138" "\n" \
"4.00260    0.99999862" "\n" \
"" "\n" \
"Li  2" "\n" \
"6.015121   0.075" "\n" \
"7.016003   0.925" "\n" \
"" "\n" \
"Be  1" "\n" \
"9.012182   1.0" "\n" \
"" "\n" \
"B  2" "\n" \
"10.012937  0.199" "\n" \
"11.009305  0.801" "\n" \
"" "\n" \
"C  2" "\n" \
"12.0000000 0.988930" "\n" \
"13.0033554 0.011070" "\n" \
"" "\n" \
"N  2" "\n" \
"14.0030732 0.996337" "\n" \
"15.0001088 0.003663" "\n" \
"" "\n" \
"O  3" "\n" \
"15.9949141 0.997590" "\n" \
"16.9991322 0.000374" "\n" \
"17.9991616 0.002036" "\n" \
"" "\n" \
"F  1" "\n" \
"18.9984032 1.0" "\n" \
"" "\n" \
"Ne  3" "\n" \
"19.992435  0.9048" "\n" \
"20.993843  0.0027" "\n" \
"21.991383  0.0925" "\n" \
"" "\n" \
"Na  1" "\n" \
"22.989767  1.0" "\n" \
"" "\n" \
"Mg  3" "\n" \
"23.985042  0.7899" "\n" \
"24.985837  0.1000" "\n" \
"25.982593  0.1101" "\n" \
"" "\n" \
"Al  1" "\n" \
"26.981539  1.0" "\n" \
"" "\n" \
"Si  3" "\n" \
"27.976927  0.9223" "\n" \
"28.976495  0.0467" "\n" \
"29.973770  0.0310" "\n" \
"" "\n" \
"P  1" "\n" \
"30.973762  1.0" "\n" \
"" "\n" \
"S  4" "\n" \
"31.972070  0.9502" "\n" \
"32.971456  0.0075" "\n" \
"33.967866  0.0421" "\n" \
"35.967080  0.0002" "\n" \
"" "\n" \
"Cl  2" "\n" \
"34.9688531 0.755290" "\n" \
"36.9659034 0.244710" "\n" \
"" "\n" \
"Ar  3" "\n" \
"35.967545  0.00337" "\n" \
"37.962732  0.00063" "\n" \
"39.962384  0.99600" "\n" \
"" "\n" \
"K  3" "\n" \
"38.963707  0.932581" "\n" \
"39.963999  0.000117" "\n" \
"40.961825  0.067302" "\n" \
"" "\n" \
"Ca  6" "\n" \
"39.962591  0.96941" "\n" \
"41.958618  0.00647" "\n" \
"42.958766  0.00135" "\n" \
"43.955480  0.02086" "\n" \
"45.953689  0.00004" "\n" \
"47.952533  0.00187" "\n" \
"" "\n" \
"Sc  1" "\n" \
"44.955910  1.0" "\n" \
"" "\n" \
"Ti  5" "\n" \
"45.952629  0.080" "\n" \
"46.951764  0.073" "\n" \
"47.947947  0.738" "\n" \
"48.947871  0.055" "\n" \
"49.944792  0.054" "\n" \
"" "\n" \
"V  2" "\n" \
"49.947161  0.00250" "\n" \
"50.943962  0.99750" "\n" \
"" "\n" \
"Cr  4" "\n" \
"49.946046  0.04345" "\n" \
"51.940509  0.83790" "\n" \
"52.940651  0.09500" "\n" \
"53.938882  0.02365" "\n" \
"" "\n" \
"Mn  1" "\n" \
"54.938047  1.0" "\n" \
"" "\n" \
"Fe  4" "\n" \
"53.939612  0.0590" "\n" \
"55.934939  0.9172" "\n" \
"56.935396  0.0210" "\n" \
"57.933277  0.0028" "\n" \
"" "\n" \
"Co  1" "\n" \
"58.933198  1.0" "\n" \
"" "\n" \
"Ni  5" "\n" \
"57.935346  0.6827" "\n" \
"59.930788  0.2610" "\n" \
"60.931058  0.0113" "\n" \
"61.928346  0.0359" "\n" \
"63.927968  0.0091" "\n" \
"" "\n" \
"Cu  2" "\n" \
"62.939598  0.6917" "\n" \
"64.927793  0.3083" "\n" \
"" "\n" \
"Zn  5" "\n" \
"63.929145  0.486" "\n" \
"65.926034  0.279" "\n" \
"66.927129  0.041" "\n" \
"67.924846  0.188" "\n" \
"69.925325  0.006" "\n" \
"" "\n" \
"Ga  2" "\n" \
"68.925580  0.60108" "\n" \
"70.924700  0.39892" "\n" \
"" "\n" \
"Ge  5" "\n" \
"69.924250  0.205" "\n" \
"71.922079  0.274" "\n" \
"72.923463  0.078" "\n" \
"73.921177  0.365" "\n" \
"75.921401  0.078" "\n" \
"" "\n" \
"As  1" "\n" \
"74.921594  1.0" "\n" \
"" "\n" \
"Se  6" "\n" \
"73.922475  0.009" "\n" \
"75.919212  0.091" "\n" \
"76.919912  0.076" "\n" \
"77.9190    0.236" "\n" \
"79.916520  0.499" "\n" \
"81.916698  0.089" "\n" \
"" "\n" \
"Br  2" "\n" \
"78.918336  0.5069" "\n" \
"80.916289  0.4931" "\n" \
"" "\n" \
"Kr  6" "\n" \
"77.914     0.0035" "\n" \
"79.916380  0.0225" "\n" \
"81.913482  0.116" "\n" \
"82.914135  0.115" "\n" \
"83.911507  0.570" "\n" \
"85.910616  0.173" "\n" \
"" "\n" \
"Rb  2" "\n" \
"84.911794  0.7217" "\n" \
"86.909187  0.2783" "\n" \
"" "\n" \
"Sr  4" "\n" \
"83.913430  0.0056" "\n" \
"85.909267  0.0986" "\n" \
"86.908884  0.0700" "\n" \
"87.905619  0.8258" "\n" \
"" "\n" \
"Y  1" "\n" \
"88.905849  1.0" "\n" \
"" "\n" \
"Zr  5" "\n" \
"89.904703  0.5145" "\n" \
"90.905644  0.1122" "\n" \
"91.905039  0.1715" "\n" \
"93.906314  0.1738" "\n" \
"95.908275  0.0280" "\n" \
"" "\n" \
"Nb  1" "\n" \
"92.906377  1.0" "\n" \
"" "\n" \
"Mo  7" "\n" \
"91.906808  0.1484" "\n" \
"93.905085  0.0925" "\n" \
"94.905840  0.1592" "\n" \
"95.904678  0.1668" "\n" \
"96.906020  0.0955" "\n" \
"97.905406  0.2413" "\n" \
"99.907477  0.0963" "\n" \
"" "\n" \
"Tc  1" "\n" \
"98.0   1.0" "\n" \
"" "\n" \
"Ru  7" "\n" \
"95.907599  0.0554" "\n" \
"97.905287  0.0186" "\n" \
"98.905939  0.127" "\n" \
"99.904219  0.126" "\n" \
"100.905582  0.171" "\n" \
"101.904348  0.316" "\n" \
"103.905424  0.186" "\n" \
"" "\n" \
"Rh  1" "\n" \
"102.905500  1.0" "\n" \
"" "\n" \
"Pd  6" "\n" \
"101.905634  0.0102" "\n" \
"103.904029  0.1114" "\n" \
"104.905079  0.2233" "\n" \
"105.903478  0.2733" "\n" \
"107.903895  0.2646" "\n" \
"109.905167  0.1172" "\n" \
"" "\n" \
"Ag  2" "\n" \
"106.905092  0.51839" "\n" \
"108.904757  0.48161" "\n" \
"" "\n" \
"Cd  8" "\n" \
"105.906461  0.0125" "\n" \
"107.904176  0.0089" "\n" \
"109.903005  0.1249" "\n" \
"110.904182  0.1280" "\n" \
"111.902758  0.2413" "\n" \
"112.904400  0.1222" "\n" \
"113.903357  0.2873" "\n" \
"115.904754  0.0749" "\n" \
"" "\n" \
"In  2" "\n" \
"112.904061  0.043" "\n" \
"114.903880  0.957" "\n" \
"" "\n" \
"Sn  10" "\n" \
"111.904826  0.0097" "\n" \
"113.902784  0.0065" "\n" \
"114.903348  0.0036" "\n" \
"115.901747  0.1453" "\n" \
"116.902956  0.0768" "\n" \
"117.901609  0.2422" "\n" \
"118.903310  0.0858" "\n" \
"119.902200  0.3259" "\n" \
"121.903440  0.0463" "\n" \
"123.905274  0.0579" "\n" \
"" "\n" \
"Sb  2" "\n" \
"120.903821  0.574" "\n" \
"122.904216  0.426" "\n" \
"" "\n" \
"Te  8" "\n" \
"119.904048  0.00095" "\n" \
"121.903054  0.0259" "\n" \
"122.904271  0.00905" "\n" \
"123.902823  0.0479" "\n" \
"124.904433  0.0712" "\n" \
"125.903314  0.1893" "\n" \
"127.904463  0.3170" "\n" \
"129.906229  0.3387" "\n" \
"" "\n" \
"I  1" "\n" \
"126.904473  1.0" "\n" \
"" "\n" \
"Xe  9" "\n" \
"123.905894  0.0010" "\n" \
"125.904281  0.0009" "\n" \
"127.903531  0.0191" "\n" \
"128.904780  0.264" "\n" \
"129.903509  0.041" "\n" \
"130.905072  0.212" "\n" \
"131.904144  0.269" "\n" \
"133.905395  0.104" "\n" \
"135.907214  0.089" "\n" \
"" "\n" \
"Cs  1" "\n" \
"132.905429  1.0" "\n" \
"" "\n" \
"Ba  7" "\n" \
"129.906282  0.00106" "\n" \
"131.905042  0.00101" "\n" \
"133.904486  0.0242" "\n" \
"134.905665  0.06593" "\n" \
"135.904553  0.0785" "\n" \
"136.905812  0.1123" "\n" \
"137.905232  0.7170" "\n" \
"" "\n" \
"La  2" "\n" \
"137.90711   0.00090" "\n" \
"138.906347  0.99910" "\n" \
"" "\n" \
"Ce  4" "\n" \
"135.907140  0.0019" "\n" \
"137.905985  0.0025" "\n" \
"139.905433  0.8843" "\n" \
"141.909241  0.1113" "\n" \
"" "\n" \
"Pr  1" "\n" \
"140.907647  1.0" "\n" \
"" "\n" \
"Nd  7" "\n" \
"141.907719  0.2713" "\n" \
"142.909810  0.1218" "\n" \
"143.910083  0.2380" "\n" \
"144.912570  0.0830" "\n" \
"145.913113  0.1719" "\n" \
"147.916889  0.0576" "\n" \
"149.920887  0.0564" "\n" \
"" "\n" \
"Pm  1" "\n" \
"145.0  1.0" "\n" \
"" "\n" \
"Sm  7" "\n" \
"143.911998  0.031" "\n" \
"146.914895  0.150" "\n" \
"147.914820  0.113" "\n" \
"148.917181  0.138" "\n" \
"149.917273  0.074" "\n" \
"151.919729  0.267" "\n" \
"153.922206  0.227" "\n" \
"" "\n" \
"Eu  2" "\n" \
"150.919847  0.478" "\n" \
"152.921225  0.522" "\n" \
"" "\n" \
"Gd  7" "\n" \
"151.919786  0.0020" "\n" \
"153.920861  0.0218" "\n" \
"154.922618  0.1480" "\n" \
"155.922118  0.2047" "\n" \
"156.923956  0.1565" "\n" \
"157.924099  0.2484" "\n" \
"159.927049  0.2186" "\n" \
"" "\n" \
"Tb  1" "\n" \
"158.925342  1.0" "\n" \
"" "\n" \
"Dy  7" "\n" \
"155.925277  0.0006" "\n" \
"157.924403  0.0010" "\n" \
"159.925193  0.0234" "\n" \
"160.926930  0.189" "\n" \
"161.926795  0.255" "\n" \
"162.928728  0.249" "\n" \
"163.929171  0.282" "\n" \
"" "\n" \
"Ho  1" "\n" \
"164.930319  1.0" "\n" \
"" "\n" \
"Er  6" "\n" \
"161.928775  0.0014" "\n" \
"163.929198  0.0161" "\n" \
"165.930290  0.336" "\n" \
"166.932046  0.2295" "\n" \
"167.932368  0.268" "\n" \
"169.935461  0.149" "\n" \
"" "\n" \
"Tm  1" "\n" \
"168.934212  1.0" "\n" \
"" "\n" \
"Yb  7" "\n" \
"167.933894  0.0013" "\n" \
"169.934759  0.0305" "\n" \
"170.936323  0.143" "\n" \
"171.936378  0.219" "\n" \
"172.938208  0.1612" "\n" \
"173.938859  0.318" "\n" \
"175.942564  0.127" "\n" \
"" "\n" \
"Lu  2" "\n" \
"174.940770  0.9741" "\n" \
"175.942679  0.0259" "\n" \
"" "\n" \
"Hf  6" "\n" \
"173.940044  0.00162" "\n" \
"175.941406  0.05206" "\n" \
"176.943217  0.18606" "\n" \
"177.943696  0.27297" "\n" \
"178.945812  0.13629" "\n" \
"179.946545  0.35100" "\n" \
"" "\n" \
"Ta  2" "\n" \
"179.947462  0.00012" "\n" \
"180.947992  0.99988" "\n" \
"" "\n" \
"W  5" "\n" \
"179.946701  0.0012" "\n" \
"181.948202  0.263" "\n" \
"182.950220  0.1428" "\n" \
"183.950928  0.307" "\n" \
"185.954357  0.286" "\n" \
"" "\n" \
"Re  2" "\n" \
"184.952951  0.3740" "\n" \
"186.955744  0.6260" "\n" \
"" "\n" \
"Os  7" "\n" \
"183.952488  0.0002" "\n" \
"185.953830  0.0158" "\n" \
"186.955741  0.016" "\n" \
"187.955860  0.133" "\n" \
"188.958137  0.161" "\n" \
"189.958436  0.264" "\n" \
"191.961467  0.410" "\n" \
"" "\n" \
"Ir  2" "\n" \
"190.960584  0.373" "\n" \
"192.962917  0.627" "\n" \
"" "\n" \
"Pt  6" "\n" \
"189.959917  0.0001" "\n" \
"191.961019  0.0079" "\n" \
"193.962655  0.329" "\n" \
"194.964766  0.338" "\n" \
"195.964926  0.253" "\n" \
"197.967869  0.072" "\n" \
"" "\n" \
"Au  1" "\n" \
"196.966543  1.0" "\n" \
"" "\n" \
"Hg  7" "\n" \
"195.965807  0.0015" "\n" \
"197.966743  0.100" "\n" \
"198.968254  0.169" "\n" \
"199.968300  0.231" "\n" \
"200.970277  0.132" "\n" \
"201.970617  0.298" "\n" \
"203.973467  0.0685" "\n" \
"" "\n" \
"Tl  2" "\n" \
"202.972320  0.29524" "\n" \
"204.974401  0.70476" "\n" \
"" "\n" \
"Pb  4" "\n" \
"203.973020  0.014" "\n" \
"205.974440  0.241" "\n" \
"206.975872  0.221" "\n" \
"207.976627  0.524" "\n" \
"" "\n" \
"Bi  1" "\n" \
"208.980374  1.0" "\n" \
"" "\n" \
"Po  1" "\n" \
"209.0  1.0" "\n" \
"" "\n" \
"At  1" "\n" \
"210.0  1.0" "\n" \
"" "\n" \
"Rn  1" "\n" \
"222.0  1.0" "\n" \
"" "\n" \
"Fr  1" "\n" \
"223.0  1.0" "\n" \
"" "\n" \
"Ra  1" "\n" \
"226.025  1.0" "\n" \
"" "\n" \
"Ac  1" "\n" \
"227.028  1.0" "\n" \
"" "\n" \
"Th  1" "\n" \
"232.038054  1.0" "\n" \
"" "\n" \
"Pa  1" "\n" \
"231.0359  1.0" "\n" \
"" "\n" \
"U  3" "\n" \
"234.040946  0.000055" "\n" \
"235.043924  0.00720" "\n" \
"238.050784  0.992745" "\n" \
"" "\n" \
"Np  1" "\n" \
"237.048  1.0" "\n" \
"" "\n" \
"Pu  1" "\n" \
"244.0  1.0" "\n" \
"" "\n" \
"Am  1" "\n" \
"243.0  1.0" "\n" \
"" "\n" \
"Cm  1" "\n" \
"247.0  1.0" "\n" \
"" "\n" \
"Bk  1" "\n" \
"247.0  1.0" "\n" \
"" "\n" \
"Cf  1" "\n" \
"251.0  1.0" "\n" \
"" "\n" \
"Es  1" "\n" \
"252.0  1.0" "\n" \
"" "\n" \
"Fm  1" "\n" \
"257.0  1.0" "\n" \
"" "\n" \
"Md  1" "\n" \
"258.0  1.0" "\n" \
"" "\n" \
"No  1" "\n" \
"259.0  1.0" "\n" \
"" "\n" \
"Lr  1" "\n" \
"260.0  1.0" "\n" \
"" "\n"

void IsotopeCalculator::computeIsotopeSpectrum( Spectrum &output, const romol_ptr_t mol, long charge ){
	
	FormMap fm;
	setFormulaMap( fm, mol );
	
	//Adjust H count for charge
	int mol_q = RDKit::MolOps::getFormalCharge( *mol.get() );
	if( mol_q != charge ){
		if( mol_q == 0 && charge > 0 )
			fm[em["H"]] += 1;
		if( mol_q == 0 && charge < 0 )
			fm[em["H"]] -= 1;
	}

	Pattern result(1);
	Pattern tmp;

    // initialize the result
    result.resize(1);
    result.front().mass = 0.0;
    result.front().rel_area = 1.0;

	calculate(tmp, result, fm, 0, charge);

    if(verbose) print_pattern(result, 10);
	print_to_output(output, result);
	output.normalizeAndSort();

}

//This function sets the formula map (atom_idx -> count) structure for the input molecule
void IsotopeCalculator::setFormulaMap( FormMap &output, const romol_ptr_t mol ){

	if( em.find("H") == em.end() ){
		std::cout << "Could not find symbol H in ElemMap" << std::endl;
		throw IsotopeCalculationException();	
	}
	std::string hstr = "H";
	size_t hidx = em[hstr];
	output[hidx] = 0;

	RDKit::ROMol::AtomIterator ait = mol.get()->beginAtoms();
	for( ; ait!=mol.get()->endAtoms(); ++ait){
		std::string symbol = (*ait)->getSymbol();
		if( em.find(symbol) == em.end() ){
			std::cout << "Could not find symbol " << symbol << " in ElemMap" << std::endl;
			throw IsotopeCalculationException();
		}
		size_t atomidx = em[symbol];
		if( output.find(atomidx) == output.end() ) output[atomidx] = 1;
		else output[atomidx] += 1;

		output[hidx] += (*ait)->getTotalNumHs();
	}
}

//The remainder of the functions are copied directly from emass (with minor mods)
//-------------------------------------------------------------------------------
/*This collective work is Copyright (C)2005 by Perttu Haimi 
Individual portions may be copyright by individual 
contributors, and are included in this collective work with 
permission of the copyright owners.

All rights reserved.

Redistribution and use in source and binary forms, 
with or without modification, are permitted provided 
that the following conditions are met:

    * Redistributions of source code must retain the 
      above copyright notice, this list of conditions 
      and the following disclaimer.
    * Redistributions in binary form must reproduce 
      the above copyright notice, this list of conditions 
      and the following disclaimer in the documentation 
      and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors 
      may be used to endorse or promote products derived 
      from this software without specific prior written 
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND 
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//-------------------------------------------------------------------------------

void IsotopeCalculator::init_data(){
  
  std::istringstream f(ISOTOPES);

  sad.clear();
  em.clear();

  std::string line;
  ulong elemindex = 0;
  int state = 0;
  while(getline(f, line)) {
    std::istringstream ist(line);
    std::string element;
    switch(state) {
    case 0: // new element
      ist >> element;
      em[element] = elemindex;
      sad.push_back(SuperAtomList(1));
      sad.back().reserve(8); // reserve room for 8 superatoms
      elemindex++;
      state = 1;
      break;
    case 1: // isotope
      ipeak p;
      Pattern & idist = sad.back()[0];
      if(ist >> p.mass >> p.rel_area) {
	// fill the gaps in the patterns with zero abundancy peaks
	if(idist.size() > 0) {
	  double prevmass = idist.back().mass;
	  for(int i = 0; i < int(p.mass - prevmass - 0.5); i++) {
	    ipeak filler;
	    filler.mass = DUMMY_MASS;
	    filler.rel_area = 0;
	    idist.push_back(filler);
	  }
	}
	// insert the peak
	idist.push_back(p);                                                        
      } else  
	state = 0; // no more isotope data
      break;
    }
  }  
}


// Merge two patterns to one.
void IsotopeCalculator::convolute_basic(Pattern & h, const Pattern & g, const Pattern & f)
{
  h.clear();
  size_t g_n = g.size();
  size_t f_n = f.size();
  if(g_n == 0 || f_n == 0)
     return;
  for(size_t k = 0; k < g_n + f_n - 1; k++) {
    double sumweight = 0, summass = 0;
    size_t start = k < (f_n - 1) ? 0 : k - f_n + 1; // max(0, k-f_n+1)
    size_t end = k < (g_n - 1) ? k : g_n - 1;       // min(g_n - 1, k)
    for(size_t i = start; i <= end; i++) {
      double weight = g[i].rel_area * f[k - i].rel_area;
      double mass = g[i].mass + f[k - i].mass;
      sumweight += weight;
      summass += weight * mass;
    }
    ipeak p;
    if(sumweight == 0)
      p.mass = DUMMY_MASS;
    else
      p.mass = summass / sumweight;
    p.rel_area = sumweight;
    h.push_back(p);
  }
}

// Prune the small peaks from both sides but
// leave them within the pattern.
void IsotopeCalculator::prune(Pattern & f, double limit)
{
  // prune the front
  Pattern::iterator i = f.begin();
  while(i != f.end()) {
    if((*i).rel_area > limit)
      break;
    i++;
  }
  f.erase(f.begin(), i);

  // prune the end
  while(1) {
    if(f.size() == 0)
      break;
    if(f.back().rel_area > limit)
      break;
    f.pop_back();
  } 
}

void IsotopeCalculator::calculate(Pattern & tmp, Pattern & result, FormMap & fm, double limit, long charge)
{
  for(FormMap::iterator i = fm.begin(); i != fm.end(); i++) {
    size_t atom_index = (*i).first;
    SuperAtomList sal = sad[atom_index];
    ulong n = (*i).second;
    ulong j = 0;
    while(n > 0) {
      size_t sz = sal.size();
      if(j == sz) { 
	sal.resize(sz + 1); // Make new superatom from previous
                            // largest superatom. We are trying to
                            // avoid copying on assignment here.
	convolute_basic(sal[j], sal[j - 1], sal[j - 1]);
	prune(sal[j], limit);
      }
      if(n & 1) { // digit is 1, convolute result
	convolute_basic(tmp , result, sal[j]);
	prune(tmp, limit);
	swap(tmp, result);  // Hopefully the swap implementation
                            // will not copy all elements.
      }
      n >>= 1; 
      j++;
    }
  }
  
  // take charge into account
  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    if(charge > 0)
      (*i).mass = (*i).mass / abs(charge) - MASS_ELECTRON;
    else if (charge < 0)
      (*i).mass = (*i).mass / abs(charge) + MASS_ELECTRON;
  }
}


void IsotopeCalculator::print_pattern(Pattern & result, int digits)
{
  // find the maximum
  double max_area = 0;
  double sum_area = 0;
  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    if(max_area < (*i).rel_area)
      max_area = (*i).rel_area;
    sum_area += (*i).rel_area;
  }
  if(max_area == 0)
    return; // empty pattern

  std::wcout.setf(std::ios::fixed);
  std::wcout.precision(digits);
  double print_limit = pow(10.0, -digits) / 2;
  //wcout.precision(30);
  //double print_limit = 0.000001 / 2;
  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    double mass = (*i).mass;
    double rel_area = (*i).rel_area;
    double val_perc = rel_area / max_area * 100;
    //double val_norm = rel_area / sum_area;
    if(mass != DUMMY_MASS && val_perc >= print_limit)
      std::wcout << mass << L" " << val_perc << std::endl;
    //wcout << mass << L" " << val_perc << L" " << val_norm << endl;
  }
}

void IsotopeCalculator::print_to_output(Spectrum & output, Pattern & result)
{
  // find the maximum
  double max_area = 0;
  double sum_area = 0;
  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    if(max_area < (*i).rel_area)
      max_area = (*i).rel_area;
    sum_area += (*i).rel_area;
  }
  if(max_area == 0)
    return; // empty pattern

  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    double mass = (*i).mass;
    double rel_area = (*i).rel_area;
    double val_perc = rel_area / max_area * 100;
	if(mass != DUMMY_MASS && val_perc >= intensity_thresh)
		output.push_back( Peak(mass, rel_area) );
  }

  if( output.size() == 0 ) throw EmptyIsotopeSpectrumException();
  
}
