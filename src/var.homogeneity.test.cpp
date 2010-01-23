//#=====================================================================================
//#
//#       Filename: var.homogeneity.test.cpp 
//#
//#    Description:  Function to calculate homogenety of variance.
//#
//#        Version:  0.1
//#        Created:  06-Apr-2009
//#       Revision:  none
//#       
//#
//#         Author:  Maksim V. Struchalin
//#        Company:  ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.
//#          Email:  m.struchalin@@erasmusmc.nl
//#
//#=====================================================================================

#include "gtps_container.h"
#include "bartlett_test.h"
#include <cstdlib>  
#include <math.h>
#include <map>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <Rinternals.h>



#include <vector>
#include <vector>
#include <map>
#include <iterator>

#include <R.h>  // to include Rconfig.h 

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("pkg", String)
// replace pkg as appropriate 
#else
#define _(String) (String)
#endif


extern "C" {

enum {AAvsABvsBB, AAvsABandBB, ABvsAAandBB, BBvsAAandAB};



void bartlett_test_R(char *gtps, int * nids,  int *nsnps,
			 						 double * trait, int * is_trait_na,
									 double * chi2, char** analys_type_)
{

int analys_type = AAvsABvsBB;


if(*analys_type_ == "AAvsABvsBB")
	{
	analys_type = AAvsABvsBB;
	}
if(*analys_type_ == "AAvsABandBB")
	{
	analys_type = AAvsABandBB;
	}
if(*analys_type_ == "ABvsAAandBB")
	{
	analys_type = ABvsAAandBB;
	}
if(*analys_type_ == "BBvsAAandAB")
	{
	analys_type = BBvsAAandAB;
	}


gtps_container gwa_data(gtps, *nids, *nsnps);

std::vector<double> NA, AA, AB, BB; // here we will store trait for different genotype group

unsigned counter=0, step=100000;
unsigned counter1 = step;



int genotype;

for(unsigned snp=0 ; snp < *nsnps ; snp++)
	{
	for(unsigned id=0 ; id < *nids ; id++)
		{	
		if(is_trait_na[id] == 1) continue;

		//spread ids trait among genotype group
		genotype = int(gwa_data.get(id+1, snp+1));

		

		switch(genotype)
			{
			case 0:
				{
				//NA.push_back(trait[id]);
				}
			break;
			case 1:
				{
				AA.push_back(trait[id]);
				}
			break;
			case 2:
				{
				AB.push_back(trait[id]);
				}
			break;
			case 3:
				{
				BB.push_back(trait[id]);
				}
			break;
			default:
				{
				Rprintf("error: bartlett_test_R: wrong genortype returned from gtps_container\n");
				return;	
				}
			break;
			}
		}

	
	
//_________________________________________________
//1df conversion:
	//b.insert(b.end(), a.begin(), a.end());
  //enum {AAvsABvsBB, AAvsABandBB, ABvsAAandBB, BBvsAAandAB};

	if(analys_type == AAvsABandBB)
		{
		AB.insert(AB.end(), BB.begin(), BB.end());
		BB.clear();
		}
	if(analys_type == ABvsAAandBB)
		{
		AA.insert(AA.end(), BB.begin(), BB.end());
		BB.clear();
		}
	if(analys_type == BBvsAAandAB)
		{
		AB.insert(AB.end(), AA.begin(), AA.end());
		AA.clear();
		}


//_________________________________________________

	unsigned NA_size = NA.size();
	unsigned AA_size = AA.size();
	unsigned AB_size = AB.size();
	unsigned BB_size = BB.size();

	double *na, *aa, *ab, *bb;


	std::list<my_small_vector> trait_groups;







		if(AA_size > 1)
			{
			aa = new double[AA_size];
			for(unsigned i=0 ; i<AA_size ; i++)
				{
				aa[i] = AA[i];
				}
			trait_groups.push_back(my_small_vector(aa, AA_size));
			}



		if(AB_size > 1)
			{
			ab = new double[AB_size];
			for(unsigned i=0 ; i<AB_size ; i++)
				{
				ab[i] = AB[i];
				}
			trait_groups.push_back(my_small_vector(ab, AB_size));
			}
	

	
		if(BB_size > 1)
			{
			bb = new double[BB_size];
			for(unsigned i=0 ; i<BB_size ; i++)
				{
				bb[i] = BB[i];
				}
			trait_groups.push_back(my_small_vector(bb, BB_size));
			}








		//Start analysis of homogeneity
		if(trait_groups.size() >= 2)
			{
			chi2[snp] = bartlett_test(&trait_groups);	
			}
		else
			{
			chi2[snp] = -1;
			}

//if(chi2[snp] == -1) {std::cout<<"snp="<<snp<<" (chi2=-1), AA_size="<<AA_size<<", AB_size="<<AB_size<<", BB_size="<<BB_size<<"\n";}
//if(chi2[snp] == -1) {std::cout<<"AA.size()="<<AA.size()<<"AB.size()="<<AB.size()<<"BB.size()="<<BB.size()<<"\n";}



if(AA_size > 1) delete[] aa;
if(AB_size > 1) delete[] ab;
if(BB_size > 1) delete[] bb;


AA.clear();
AB.clear();
BB.clear();



counter++;
if(counter == counter1)
	{
	counter1 += step;
	std::cout<<"\t\t\t\t\t\t\t\t\t"<<counter<<" SNPs done\n";	
	}


	}




} // end of function bartlett_test









}// end of extern "C"
