//#=====================================================================================
//#
//#       Filename:  var_meta_plink_C.cpp
//#
//#    Description:  Function for meta analysis ov variance. Read plink format flat files and perform metanalysis.
//#
//#        Version:  1.0
//#        Created:  02-July-2009
//#       Revision:  none
//#
//#
//#         Author:  Maksim V. Struchalin
//#        Company:  ErasmusMC, Epidemiology, The Netherlands.
//#          Email:  m.struchalin@erasmusmc.nl
//#
//#=====================================================================================


//Description:
//
//Function
//"void var_meta_plink_C(const char** filenames, unsigned *file_amount_, const char **output_filename_, unsigned *skip_first_lines_amount_, char **delim_)"
//is wrtitten to meta analyze variance of snps. Input parameteres are 
//filenames - massive of file names where snp variances store. These files came from plink in cases of running it with key --qt-means like in bellow example
//plink --bfile $genotypes --allow-no-sex --pheno $phe_filename --all-pheno --missing-phenotype -999 --qt-means --out out --assoc
//input file format of such file is:
//CHR          SNP  VALUE      G11      G12      G22
//1   rs11497407   GENO      A/A      A/G      G/G
//1   rs11497407 COUNTS        0        2     5699
//1   rs11497407   FREQ        0 0.0003508   0.9996
//1   rs11497407   MEAN       NA  -0.2434 8.542e-05
//1   rs11497407     SD       NA   0.6545        1
//1   rs12565286   GENO      C/C      C/G      G/G
//1   rs12565286 COUNTS        3      459     5239
//1   rs12565286   FREQ 0.0005262  0.08051    0.919
//1   rs12565286   MEAN   0.2318 -0.02798 0.002319
//1   rs12565286     SD   0.1208    0.966    1.003
//
//Function var_meta_plink_C is robust enougph to little deviation from this format: 
//1) order of columns can be arbitrary
//2) there can be any other collumns. unnecessary columns will be ignored
//3) same is with unnecessary rowa. it is ignored if it is unnecessary


//input parameter "file_amount_" is number of files
//output_filename_ - output filename. All metaanalysed snps will be there
//skip_first_lines_amount_ - how many snps should be skipped in input files. 0 by default. In cases somebody wants use differing from plink format
//delim_ - which delim is used, blank space (' ') by default.





#include "dometa.h"
#include <map>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <Rinternals.h>

#include <vector>

#include <R.h>  // to include Rconfig.h

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("pkg", String)
// replace pkg as appropriate
#else
#define _(String) (String)
#endif


#define VARIABLE_TYPE double
#define GENO_TYPES_NUM 3


#include <cmath>

		


const std::string chromosome_column_name = "CHR";
const std::string snpname_column_name = "SNP";
const std::string value_column_name = "VALUE";
const std::string g11_column_name = "G11";
const std::string g12_column_name = "G12";
const std::string g22_column_name = "G22";

const std::string geno_value_name = "GENO";
const std::string counts_value_name = "COUNTS";
const std::string freq_value_name = "FREQ";
const std::string mean_value_name = "MEAN";
const std::string sd_value_name = "SD";



const VARIABLE_TYPE NA_value = -999.999; // which numeric value is ised as NA
const int NA_value_int = -999; // which numeric value is ised as NA
const unsigned precision_output=4;



//all info for a SNP is here
//_________________________________________________________
struct snp_var_data
	{
	snp_var_data() //set all variables to zero
		{
		reset();
		}	
			
	inline void reset() 
		{
		for(int i=0 ; i<GENO_TYPES_NUM ; i++)
			{
			GENO[i] = "NA";
			COUNTS[i] = NA_value_int;
			MEAN[i] = NA_value;
			SD[i]=NA_value;
			Z=NA_value;
			}
		chromosome=NA_value_int;
		snpname="";
		}
	
	std::string snpname;	
	std::string GENO[GENO_TYPES_NUM];
	int COUNTS[GENO_TYPES_NUM];
	VARIABLE_TYPE MEAN[GENO_TYPES_NUM];
	VARIABLE_TYPE SD[GENO_TYPES_NUM];
	int chromosome;
	double Z; //homogeneity test
	double Z_2df; //homogeneity test 2df only
	};
//_________________________________________________________





typedef std::map<std::string, snp_var_data*> Snp_store_type; //use std::string as key because for some reason it doesn't want to work with const char*

bool include_snp(Snp_store_type *, snp_var_data*); //include snp into common storage
snp_var_data* snp_var_meta(snp_var_data* , snp_var_data*); //metaanalysis of two snps
bool is_na(const VARIABLE_TYPE val, const VARIABLE_TYPE na_reference=NA_value); //is numerical value recognized as NA 
void save_snps_data_into_file(Snp_store_type *snps_data, const char *output_filename, char delim); //save all snps into flat file in plink like format
void save_snps_tests_into_file(Snp_store_type *snps_data, const char *output_filename, char delim); //save all snps into flat file in plink like format
std::string double_2_str(VARIABLE_TYPE val, const unsigned precision=precision_output); // convert double to string
bool check(snp_var_data* snp1, snp_var_data* snp2);

bool unify_snp(snp_var_data* snp);
std::string get_uniq_symbols(std::string alleles_snp);

double perform_bartlett_test_for_snp(snp_var_data * snp, bool all_genogroup_only);

double my_median(std::vector<double> * vec);




//Start of function body 
//___________________________________________________________
//Function get names of files where are means and variance for each geno group. Return metanalysed variances and means.
//Function get as many functions as you have. Metanalysis goes file by file. 
//Files must have columns CHR, SNP, VALUE, G11, G12, G22. And each snp has to have parameters
//GENO, COUNTS, MEAN, SD (order isn't important). And don't argue with me!!!

extern "C" {
void var_meta_plink_C(const char** filenames, unsigned *file_amount_, const char **output_filename_, unsigned *skip_first_lines_amount_, char **delim_,
											double* lambdas)
{
unsigned file_amount = *file_amount_;
unsigned skip_first_lines_amount = *skip_first_lines_amount_;
char delim = *delim_[0];

double lambda; //inflation factor


std::string output_filename = output_filename_[0];


		
Snp_store_type snps_data; //all analysed data will be here


Rprintf("Each file has to contain columns with names %s, %s, %s, %s, %s, %s\n", chromosome_column_name.c_str(),
																																								snpname_column_name.c_str(),
																																								value_column_name.c_str(),
																																								g11_column_name.c_str(),
																																								g12_column_name.c_str(),
																																								g22_column_name.c_str());



//read the first file
std::ifstream file;

std::stringstream chromosome_stream, snpname_stream, g11_value_stream, g12_value_stream, g22_value_stream;
std::stringstream num_to_string;
std::string value;
std::string str_from_stream;

std::vector<double> bartlet_tests_vec;//, lambda_vec;





for(unsigned file_num=0 ; file_num < file_amount ; file_num++)
	{
	Rprintf("\nProcessing file \"%s\"...\n", filenames[file_num]);
	file.open(filenames[file_num]);
	if(!file.is_open()){error("Can not open file %s\n", filenames[file_num]);}

	//skip first line
	for(unsigned i=0 ; i<skip_first_lines_amount ; i++) 
		{
		getline(file, str_from_stream);
		if(file.eof()) {error("Tried to skip %i lines in file %s but there is %i at all ", skip_first_lines_amount, filenames[file_num], i);}
		}

	//read header and determine position of our columns

	int CHR_position=-1, VALUE_position=-1, SNP_position=-1, G11_position=-1, G12_position=-1, G22_position=-1; //0 means te first column

	getline(file, str_from_stream);
	std::stringstream line_stream(str_from_stream);

	if(file.eof()) break;

	for(unsigned col=0 ; !line_stream.eof() ; col++ )
		{
		getline(line_stream, str_from_stream, delim);
		if(str_from_stream.size() == 0) {col--; continue;}
	

		if(str_from_stream == chromosome_column_name) {CHR_position=col;}
		else if(str_from_stream == snpname_column_name) {SNP_position=col;}
		else if(str_from_stream == value_column_name) {VALUE_position=col;}
		else if(str_from_stream == g11_column_name) {G11_position=col;}
		else if(str_from_stream == g12_column_name) {G12_position=col;}
		else if(str_from_stream == g22_column_name) {G22_position=col;}
		}

	if(CHR_position==-1) {error("Can not find column \"%s\"\n", chromosome_column_name.c_str());}
	if(SNP_position==-1) {error("Can not find column \"%s\"\n", snpname_column_name.c_str());}
	if(VALUE_position==-1) {error("Can not find column \"%s\"\n", value_column_name.c_str());}
	if(G11_position==-1) {error("Can not find column \"%s\"\n", g11_column_name.c_str());}
	if(G12_position==-1) {error("Can not find column \"%s\"\n", g12_column_name.c_str());}
	if(G22_position==-1) {error("Can not find column \"%s\"\n", g22_column_name.c_str());}



	//start reading data from file
	line_stream.str("");
	line_stream.clear(); //set eof to false
	
	static bool GENO_bool, COUNTS_bool, MEAN_bool, SD_bool;
	GENO_bool = COUNTS_bool = MEAN_bool = SD_bool = false;


	static unsigned snp_number;
	snp_number=0;
	
	static unsigned current_line_number;
	current_line_number=0;

	static unsigned step;
	step = 100000;

	while(1) // one iteration = one snp
		{
		snp_var_data *snp = new snp_var_data;	

		snp_number++;
	
		while(1) //run this cycle untill a snp is ready
			{
			getline(file, str_from_stream); //get next line from file
			
			current_line_number++;
					
			if(file.eof()) break; //exit if EOF
				
			line_stream.str(str_from_stream);	
			
			for(unsigned col=0; !line_stream.eof() ; col++ )  //start spliting current line
				{
				getline(line_stream, str_from_stream, delim);
			  if(str_from_stream.size() == 0) {col--; continue;}
				
				if(str_from_stream == "NA" || str_from_stream == "na" || str_from_stream == "nan" || str_from_stream == "NaN") 
					{
					num_to_string << NA_value;
					str_from_stream = num_to_string.str();
					num_to_string.str("");
					num_to_string.clear();
					}
			
				
				if(col == CHR_position)				 {chromosome_stream << str_from_stream;}
				else if(col == SNP_position) 	 {snpname_stream << str_from_stream;}
				else if(col == VALUE_position) {value = str_from_stream;}
				else if(col == G11_position)   {g11_value_stream << str_from_stream;}
				else if(col == G12_position)   {g12_value_stream << str_from_stream;}
				else if(col == G22_position)   {g22_value_stream << str_from_stream;}
				}
			//Now the line has been splited. Let's try to recognize what we splited and put it into snp info storage	

			
			//check wether we read same snp as in prevoius iteration. If this is new SNP then skip snp with incomplete info
			if(snp->snpname == "")
				{
				snp->snpname = snpname_stream.str();
				chromosome_stream >> snp->chromosome;
				}
			else if(snp->snpname != snpname_stream.str() && !(GENO_bool && COUNTS_bool && MEAN_bool && SD_bool))
				{
				Rprintf("Attention. Can not find complete information for SNP \"%s\". Line %i in file. SNP skiped\n", snp->snpname.c_str(), current_line_number);
				snp_number--;
				//Start reading new snp
				snp->reset();
				snp->snpname = snpname_stream.str();
				chromosome_stream >> snp->chromosome;
				GENO_bool=false, COUNTS_bool=false, MEAN_bool=false, SD_bool=false;
				}
			
			
			if(value == geno_value_name) {g11_value_stream >> snp->GENO[0]; g12_value_stream >> snp->GENO[1]; g22_value_stream >> snp->GENO[2]; GENO_bool = true;}
			else if(value == counts_value_name) {g11_value_stream >> snp->COUNTS[0]; g12_value_stream >> snp->COUNTS[1]; g22_value_stream >> snp->COUNTS[2]; COUNTS_bool = true;}
			else if(value == mean_value_name) {g11_value_stream >> snp->MEAN[0]; g12_value_stream >> snp->MEAN[1]; g22_value_stream >> snp->MEAN[2]; MEAN_bool = true;}
			else if(value == sd_value_name) {g11_value_stream >> snp->SD[0]; g12_value_stream >> snp->SD[1]; g22_value_stream >> snp->SD[2]; SD_bool = true;}

			g11_value_stream.str(""); g11_value_stream.clear();
			g12_value_stream.str(""); g12_value_stream.clear();
			g22_value_stream.str(""); g22_value_stream.clear();
			chromosome_stream.str(""); chromosome_stream.clear();
			snpname_stream.str(""); snpname_stream.clear();
			
			line_stream.str(""); line_stream.clear(); //set eof to false
			
			value="";


			

			if(GENO_bool && COUNTS_bool && MEAN_bool && SD_bool)
 				{
				static bool if_snp_good;
				if_snp_good = unify_snp(snp);
				if(!if_snp_good) 
					{
					break; //Start reading new SNP
					}
			
				if_snp_good = include_snp(&snps_data, snp); // put new snp into the storage. If this snp is there already tham metaanalyse it
				
				if(if_snp_good)
					{	
					static double two_df_test;	
					two_df_test = perform_bartlett_test_for_snp(snp, true);
					if(!is_na(two_df_test)) {bartlet_tests_vec.push_back(two_df_test);}
					}
				GENO_bool=false, COUNTS_bool=false, MEAN_bool=false, SD_bool=false;
				break; //Start reading new SNP
				}
		
			}


		if(snp_number % step == 0) 
			{
			Rprintf("%i SNPs done\n", snp_number);
			if(step >= step*5) step *= 5;
			}


		if(file.eof()) {snp_number--; break;} //exit if EOF
		}
		
	file.close();
	file.clear();
	

	lambda = my_median(&bartlet_tests_vec);
//	lambda_vec.push_back(lambda);

	lambdas[file_num] = lambda;

	bartlet_tests_vec.clear();

	Rprintf("All SNPs done. Total amount of SNPs is %i\n", snp_number);
	Rprintf("Inflation factor for 2df bartlett tests for snps from \"%s\" is lambda=%2f \n", filenames[file_num], lambda);

	} // all files are read and snp pooled



//bartlett test for pooled snps

Snp_store_type::iterator iter_map;


bartlet_tests_vec.clear();

for(Snp_store_type::const_iterator i=snps_data.begin() ; i!=snps_data.end() ; ++i)
	{
	i->second->Z = perform_bartlett_test_for_snp(i->second, false);
	i->second->Z_2df = perform_bartlett_test_for_snp(i->second, true);	
	if(!is_na(i->second->Z_2df)) bartlet_tests_vec.push_back(i->second->Z_2df);
	}

lambda = my_median(&bartlet_tests_vec);
lambdas[file_amount] = lambda;
Rprintf("inflation factor for 2df bartlett tests for pooled snps is lambda=%2f \n", lambda);



std::string output_means_filename = output_filename + ".means";
std::string output_tests_filename = output_filename + ".tests";

Rprintf("\nwriting pooled means and variances into file %s\n", output_means_filename.c_str());
save_snps_data_into_file(&snps_data, output_means_filename.c_str(), delim);


Rprintf("\nwriting variance homogeneity tests results into file %s\n", output_tests_filename.c_str());
save_snps_tests_into_file(&snps_data, output_tests_filename.c_str(), delim);

Rprintf("Done\n");

}
//End of var_meta_plink_C
//___________________________________________________________




}







//___________________________________________________________
//Set genotypes in common view GA -> AG
bool unify_snp(snp_var_data* snp)
{


static std::string alleles_snp;
alleles_snp="";

for(unsigned i=0 ; i<GENO_TYPES_NUM ; i++) {alleles_snp += snp->GENO[i];}//collect all alleles 

alleles_snp = get_uniq_symbols(alleles_snp);
for(unsigned i=0 ; i<alleles_snp.size() ; i++) if(alleles_snp[i] == '/') {alleles_snp.erase(i,1);}

if(alleles_snp.size() > 2 && alleles_snp.size() <=1)
	{
	Rprintf("SNP %s has more than 2 alleles (it has %s). SNP skiped.\n", snp->snpname.c_str(), alleles_snp.c_str());
	return false;	
	}

if(alleles_snp.size() <=1)
	{
	Rprintf("SNP %s has less than 2 alleles (it has %s). SNP skiped.\n", snp->snpname.c_str(), alleles_snp.c_str());
	return false;	
	}


if(alleles_snp[0]=='0' && alleles_snp[1]=='0')
	{
	Rprintf("SNP %s has undefind allels. SNP skiped..\n", snp->snpname.c_str());
	return false;
	}










std::sort(alleles_snp.begin(), alleles_snp.end());


if(alleles_snp[0] == '0')
	{
	alleles_snp[0] = alleles_snp[1];
	alleles_snp[1] = '0';
	}

static std::string geno[3];


geno[0].push_back(alleles_snp[0]);
geno[0].push_back('/');
geno[0].push_back(alleles_snp[0]);
		
geno[1].push_back(alleles_snp[0]);
geno[1].push_back('/');
geno[1].push_back(alleles_snp[1]);

geno[2].push_back(alleles_snp[1]);
geno[2].push_back('/');
geno[2].push_back(alleles_snp[1]);


static std::string geno_hemozyg_switched;
geno_hemozyg_switched.push_back(geno[1][2]);
geno_hemozyg_switched.push_back('/');
geno_hemozyg_switched.push_back(geno[1][0]);


if(snp->GENO[0] == geno[0] && snp->GENO[1] == geno[1] && snp->GENO[2] == geno[2])
	{
	geno_hemozyg_switched = "";
	geno[0] = "";
	geno[1] = "";
	geno[2] = "";
	return true; //nothing to change
	}


static snp_var_data snp_new;

//snp_new->snpname = snp->snpname;
//snp_new->chromosome = snp->chromosome;


for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{
	if(snp->GENO[i] == geno[0])
		{snp_new.GENO[0] = geno[0]; snp_new.COUNTS[0]=snp->COUNTS[i]; snp_new.MEAN[0]=snp->MEAN[i]; snp_new.SD[0]=snp->SD[i];}
	if(snp->GENO[i] == geno[1] || snp->GENO[i]==geno_hemozyg_switched) 
		{snp_new.GENO[1] = geno[1]; snp_new.COUNTS[1]=snp->COUNTS[i]; snp_new.MEAN[1]=snp->MEAN[i]; snp_new.SD[1]=snp->SD[i];}
	if(snp->GENO[i] == geno[2]) 
		{snp_new.GENO[2] = geno[2]; snp_new.COUNTS[2]=snp->COUNTS[i]; snp_new.MEAN[2]=snp->MEAN[i]; snp_new.SD[2]=snp->SD[i];}
	}


geno_hemozyg_switched = "";
geno[0] = "";
geno[1] = "";
geno[2] = "";





//exclude those geno group which has only one id
for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{
	if(is_na(snp_new.SD[i]) || is_na(snp_new.MEAN[i]) || is_na(snp_new.COUNTS[i]) || snp_new.COUNTS[i]<=1 )
		{
		snp_new.SD[i] = NA_value;
		snp_new.MEAN[i] = NA_value;
		snp_new.COUNTS[i] = 0;
		continue;
		}

	static VARIABLE_TYPE var;
	var = snp_new.SD[i]*snp_new.SD[i];	
	if(var <= 1.E-32)
		{
		snp_new.SD[i] = NA_value;
		snp_new.MEAN[i] = NA_value;
		snp_new.COUNTS[i] = 0;
		
		Rprintf("warning: genotypic group %s in snp %s has too small variance (variance=%f). This genotypic group is excluded from analysis.\n", 
		snp_new.GENO[i].c_str(), snp->snpname.c_str(), var);
		
		}
	}




for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{
	snp->GENO[i] = snp_new.GENO[i];
	snp->COUNTS[i] = snp_new.COUNTS[i];
	snp->MEAN[i] = snp_new.MEAN[i];
	snp->SD[i] = snp_new.SD[i];
	}

snp_new.reset();




return true;
}
//___________________________________________________________




//___________________________________________________________
//put new snp into the storage. If this snp is there already tham metaanalyse it
bool include_snp(Snp_store_type * snps_storage, snp_var_data* snp)
{


Snp_store_type::iterator iter_map = snps_storage->find(snp->snpname);

char delim=' ';



if(iter_map == snps_storage->end()) 
	{
	snps_storage->insert(std::pair<std::string, snp_var_data*>(snp->snpname, snp));
	}
else
	{
	snp_var_data* snp_current = iter_map->second;
	static bool is_snp_ok;
	is_snp_ok = check(snp_current, snp);
	
	if(!is_snp_ok) {return false;}
	
	(*snps_storage)[snp->snpname.c_str()] = snp_var_meta(iter_map->second, snp);
	
	//bless god souls of these objects...
	delete snp_current;
	delete snp;
	}

return true;
}
//___________________________________________________________



//metaanalysis of data from two snps
//___________________________________________________________
snp_var_data* snp_var_meta(snp_var_data* snp1, snp_var_data* snp2)
{
snp_var_data* snp_meta = new snp_var_data;



snp_meta->chromosome = snp1->chromosome;
snp_meta->snpname = snp1->snpname;

VARIABLE_TYPE meta_var_beta[GENO_TYPES_NUM], meta_var_se[GENO_TYPES_NUM];
VARIABLE_TYPE meta_mean_beta[GENO_TYPES_NUM], meta_mean_se[GENO_TYPES_NUM];
	
VARIABLE_TYPE var1, var2, sd_var1, sd_var2;

static unsigned one = 1;

for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{
	//metanalysis for variance

	static std::string meta_codding;
			
	snp_meta->GENO[i] = snp1->GENO[i];
	
	if(snp1->COUNTS[i] != 0	&& snp2->COUNTS[i] != 0 )
		{
//		var1 = snp1->SD[i]*snp1->SD[i];
//		var2 = snp2->SD[i]*snp2->SD[i];
//		sd_var1 = var1*sqrt(2/snp1->COUNTS[i]);
//		sd_var2 = var1*sqrt(2/snp2->COUNTS[i]);
//
//		dometa_c(&var1, &var2,
//						 &sd_var1, &sd_var2,
//					 		NULL, NULL,
//							&one,
//							meta_var_beta,
//							meta_var_se);
//	
//		snp_meta->COUNTS[i] = snp1->COUNTS[i] + snp2->COUNTS[i];
//		snp_meta->SD[i] = sqrt(meta_var_beta[0]);
	

		static VARIABLE_TYPE SD_of_the_mean_snp1, SD_of_the_mean_snp2;
		
		SD_of_the_mean_snp1 = snp1->SD[i]/sqrt(double(snp1->COUNTS[i]));
		SD_of_the_mean_snp2 = snp2->SD[i]/sqrt(double(snp2->COUNTS[i]));


	
		//metanalysis for mean
		dometa_c(&snp1->MEAN[i], &snp2->MEAN[i],
						 &SD_of_the_mean_snp1, &SD_of_the_mean_snp2,
					 		NULL, NULL,
							&one,
							meta_mean_beta,
							meta_mean_se);
	
		snp_meta->MEAN[i] = meta_mean_beta[0];
		snp_meta->COUNTS[i] = snp1->COUNTS[i] + snp2->COUNTS[i];
		snp_meta->SD[i] = meta_mean_se[0]*sqrt(double(snp_meta->COUNTS[i]));
		}
	else
		{
		if(snp1->COUNTS[i] == 0)
			{
			snp_meta->SD[i] = snp2->SD[i];
			snp_meta->COUNTS[i] = snp2->COUNTS[i];
			snp_meta->MEAN[i] = snp2->MEAN[i];
			}
		else if(snp2->COUNTS[i] == 0)
			{
			snp_meta->SD[i] = snp1->SD[i];
			snp_meta->COUNTS[i] = snp1->COUNTS[i];
			snp_meta->MEAN[i] = snp1->MEAN[i];
			}
		else
			{
			error("Upss... something strange occured... wrong file format probably. Create small example of your data and send to developer.\n");
			}
				
		}
	
	
	}
		
return snp_meta;		
}
//___________________________________________________________







//___________________________________________________________
bool is_na(const VARIABLE_TYPE val, const VARIABLE_TYPE na_reference)
{
static VARIABLE_TYPE delta = fabs(na_reference/1E6);
if(val > na_reference-delta && val < na_reference+delta) return true;
else return false;
}
//___________________________________________________________








//save all data into plink like format
//___________________________________________________________
void save_snps_data_into_file(Snp_store_type *snps_data, const char *output_filename, char delim)
{
//print header
	
const unsigned pp_maxsnp=10;
		
std::ofstream file;
	

file.open(output_filename);


if(!file.is_open()){error("Can not open file %s\n", output_filename);}



//file<<setiosflags(std::ios::right);

file.precision(precision_output);

file << std::setw(4) << chromosome_column_name.c_str() << delim
		 << std::setw(pp_maxsnp) << snpname_column_name.c_str() << delim
		 << std::setw(6) << value_column_name.c_str() << delim
		 << std::setw(8) << g11_column_name.c_str() << delim
		 << std::setw(8) << g12_column_name.c_str() << delim
		 << std::setw(8) << g22_column_name.c_str() << "\n";


		
static Snp_store_type::iterator iter_map;

for(Snp_store_type::const_iterator i=snps_data->begin() ; i!=snps_data->end() ; ++i)
	{
	//geno:
	file << std::setw(4) << i->second->chromosome << delim
		   << std::setw(pp_maxsnp) << i->second->snpname << delim 
			 << std::setw(6) << geno_value_name << delim
			 << std::setw(8) << i->second->GENO[0] << delim
			 << std::setw(8) << i->second->GENO[1] << delim
			 << std::setw(8) << i->second->GENO[2] << "\n";



	//counts:
	file << std::setw(4) << i->second->chromosome << delim
		   << std::setw(pp_maxsnp) << i->second->snpname << delim 
			 << std::setw(6) << counts_value_name << delim
			 << std::setw(8) << (is_na(i->second->COUNTS[0])? "NA":double_2_str(i->second->COUNTS[0])) << delim
			 << std::setw(8) << (is_na(i->second->COUNTS[1])? "NA":double_2_str(i->second->COUNTS[1])) << delim
			 << std::setw(8) << (is_na(i->second->COUNTS[2])? "NA":double_2_str(i->second->COUNTS[2])) << "\n";

	
	//freq:
	static VARIABLE_TYPE total_id_num;
 	total_id_num = i->second->COUNTS[0] + i->second->COUNTS[1] + i->second->COUNTS[2];	
	
	file << std::setw(4) << i->second->chromosome << delim
		   << std::setw(pp_maxsnp) << i->second->snpname << delim 
			 << std::setw(6) << freq_value_name << delim
			 << std::setw(8) << (is_na(i->second->COUNTS[0])? "NA":double_2_str(i->second->COUNTS[0]/total_id_num)) << delim
			 << std::setw(8) << (is_na(i->second->COUNTS[1])? "NA":double_2_str(i->second->COUNTS[1]/total_id_num)) << delim
			 << std::setw(8) << (is_na(i->second->COUNTS[2])? "NA":double_2_str(i->second->COUNTS[2]/total_id_num)) << "\n";


	//mean:
	file << std::setw(4) << i->second->chromosome << delim
		   << std::setw(pp_maxsnp) << i->second->snpname << delim 
			 << std::setw(6) << mean_value_name << delim
			 << std::setw(8) << (is_na(i->second->MEAN[0])? "NA":double_2_str(i->second->MEAN[0])) << delim
			 << std::setw(8) << (is_na(i->second->MEAN[1])? "NA":double_2_str(i->second->MEAN[1])) << delim
			 << std::setw(8) << (is_na(i->second->MEAN[2])? "NA":double_2_str(i->second->MEAN[2])) << "\n";

	//sd:
	file << std::setw(4) << i->second->chromosome << delim
		   << std::setw(pp_maxsnp) << i->second->snpname << delim 
			 << std::setw(6) << sd_value_name << delim
			 << std::setw(8) << (is_na(i->second->SD[0])? "NA":double_2_str(i->second->SD[0])) << delim
			 << std::setw(8) << (is_na(i->second->SD[1])? "NA":double_2_str(i->second->SD[1])) << delim
			 << std::setw(8) << (is_na(i->second->SD[2])? "NA":double_2_str(i->second->SD[2])) << "\n";


	}


file.close();
}
//___________________________________________________________




//___________________________________________________________
void save_snps_tests_into_file(Snp_store_type *snps_data, const char *output_filename, char delim)
{
//print header
	
const unsigned pp_maxsnp=10;
		
std::ofstream file;
	

file.open(output_filename);


if(!file.is_open()){error("Can not open file %s\n", output_filename);}



//file<<setiosflags(std::ios::right);

file.precision(precision_output);

file << snpname_column_name << delim << "Z" << delim << "Z_2df\n";


		
static Snp_store_type::iterator iter_map;

for(Snp_store_type::const_iterator i=snps_data->begin() ; i!=snps_data->end() ; ++i)
	{
	//geno:
	file << i->second->snpname << delim 
			 << (is_na(i->second->Z)? "NA":double_2_str(i->second->Z)) << delim
			 << (is_na(i->second->Z_2df)? "NA":double_2_str(i->second->Z_2df)) << "\n"; 
	}


file.close();
}
//___________________________________________________________



//convert from double to string
//___________________________________________________________
std::string double_2_str(double val, const unsigned precision)
{
static std::stringstream stream;
stream.str(""); stream.clear();
stream.precision(precision);
stream << val;
return stream.str();
}
//___________________________________________________________



//___________________________________________________________
//Check whether snps have same genotypes and in same columns. If one of snp has allele 0 but another have real than replace 0 by real one.
bool check(snp_var_data* snp2, snp_var_data* snp1)
{
//check snp name
if(snp1->snpname != snp2->snpname) 
	{
	error("snp_var_meta: unexpected error; atempt to pool two different snps");
	}

//check chromosome name
if(snp1->chromosome != snp2->chromosome)
	{
	Rprintf("warning: SNP %s has different chromosome number in different files. Previos value is %i, current one is %i. Value %i is used.\n",
				 	snp2->snpname.c_str(), snp2->chromosome, snp1->chromosome, snp2->chromosome);
	snp1->chromosome = snp2->chromosome;
	}



//check coddings
if(snp1->GENO[1] == snp2->GENO[1])
	{
	return true;
	}



//snp1 - snp from next cohort. If genotypes from snp1 and snp2 don't match than skip snp1

if(snp1->GENO[1][0] == snp2->GENO[1][0]) //A1==A2?
	{
	if(snp1->GENO[1][2] == '0') //B1==0?
		{
		if(snp2->GENO[1][2] != '0') //B2!=0?
			{
			//situation when snp1 has A/A, A/0, 0/0 and snp2 has A/A, A/B, B/B. Replace 0 by B in snp1
			snp1->GENO[1][2] = snp2->GENO[1][2];
			snp1->GENO[2][0] = snp2->GENO[1][2];
			snp1->GENO[2][2] = snp2->GENO[1][2];
			return true;
			}
		}
	else //B1!=0!!!
		{
		if(snp2->GENO[1][2] == '0') // B2==0?
			{
			//situation when snp1 has A/A, A/B, B/B and snp2 has A/A, A/0, 0/0. Replace 0 by B in snp1
			snp2->GENO[1][2] = snp1->GENO[1][2];
			snp2->GENO[2][0] = snp1->GENO[1][2];
			snp2->GENO[2][2] = snp1->GENO[1][2];
			return true;
			}
		else
			{
			//B1!=0 and B2!=0 therefore codings defer. Skip snp1
			Rprintf("warning: snp %s has different genotypes in current and previous cohort. Current one is %s, previos - %s. The current one is skiped.\n", 
							snp1->snpname.c_str(), snp1->GENO[1].c_str(), snp2->GENO[1].c_str());
			return false;
			}
		}
	
	}
else if(snp1->GENO[1][0] == snp2->GENO[1][2])
		{
		//situation when snp1 has B/B, B/C, C/C and snp2 has A/A, A/B, B/B. C can be 0. Let's check it
		if(snp1->GENO[1][2] == '0')
			{
			//situation when snp1 has B/B, B/0, 0/0 and snp2 has A/A, A/B, B/B.
			snp1->GENO[1][2] = snp2->GENO[1][0];
			snp1->GENO[2][0] = snp2->GENO[1][0];
			snp1->GENO[2][2] = snp2->GENO[1][0];
			unify_snp(snp1);
			return true;
			}
		else //B1!=0!!!
			{
			//it could be snp1 has B/B, B/A, A/A and snp2 has A/A, A/0, 0/0. But snp1 has been sorted already => second allele of snp2 is C!=0
			Rprintf("warning: snp %s has different genotypes in current and previous cohort. current is %s, previos is %s. The current one is skiped.\n", 
							snp1->snpname.c_str(), snp1->GENO[1].c_str(), snp2->GENO[1].c_str());
			return false;
			}


		}
else if(snp1->GENO[1][2] == snp2->GENO[1][0])
	{
	//situation when snp1 has A/A, A/B, B/B and snp2 has B/B, B/0, 0/0. Replace 0 by B in snp1
	snp2->GENO[1][2] = snp1->GENO[1][0];
	snp2->GENO[2][0] = snp1->GENO[1][0];
	snp2->GENO[2][2] = snp1->GENO[1][0];
	unify_snp(snp2);
	
	}
else
	{
	//A1!=A2 and A1!=B2
	Rprintf("warning: snp %s has different genotypes in current and previos cohort. current is %s, previos is %s. The current one is skiped.\n", 
					snp1->snpname.c_str(), snp1->GENO[1].c_str(), snp2->GENO[1].c_str());
	return false;
	}

return true;
}
//___________________________________________________________



//___________________________________________________________
// input parameter is "A/AA/GG/G", output AG/
std::string get_uniq_symbols(std::string alleles_snp)
{
std::string uniqe_symbols="";
int size = alleles_snp.size();

for(int i=0 ; i<size ; i++)
	{
	static char symbol;
 	symbol = alleles_snp[i];
	
	static int size_uniqe;
	size_uniqe = uniqe_symbols.size();

	static bool flag;
	flag=false;
	for(int j=0 ; j<size_uniqe ; j++)
		{
		if(uniqe_symbols[j] == symbol) {flag = true; break;}
		}
	
	if(!flag) {uniqe_symbols += symbol;}

	}


return uniqe_symbols;
}
//___________________________________________________________




//___________________________________________________________
double perform_bartlett_test_for_snp(snp_var_data * snp, bool all_genogroup_only)
{
//http://en.wikipedia.org/wiki/Bartlett%27s_test


unsigned N=0;
VARIABLE_TYPE Sp_2=0;
	
VARIABLE_TYPE sum_ni_1lnSi2=0;

VARIABLE_TYPE sum_1__n_1=0;

VARIABLE_TYPE var;

unsigned geno_group_num=0;


for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{

	if(snp->COUNTS[i] == 0) continue;

	var = snp->SD[i]*snp->SD[i];
	
	
	geno_group_num++;
	N += snp->COUNTS[i];

	
	
	sum_ni_1lnSi2 += (snp->COUNTS[i]-1.)*log(var);
	
	sum_1__n_1 += 1./(snp->COUNTS[i]-1.);
	
	Sp_2 += (snp->COUNTS[i]-1.)*var;

	}

if(all_genogroup_only && geno_group_num != GENO_TYPES_NUM) return NA_value;

Sp_2 /= N - geno_group_num; 


return ((N - geno_group_num)*log(Sp_2) - sum_ni_1lnSi2)/(1 + (sum_1__n_1-1/(N-geno_group_num))/(3*(geno_group_num-1)) );
}

//___________________________________________________________




//___________________________________________________________
double my_median(std::vector<double> * vec)
{
unsigned size = vec->size();


if(size == 0) return -1;
if(size == 1) return (*vec)[0];

std::sort(vec->begin(), vec->end());

for(int i=0 ; i<size ; i++)
	{
	}


static double median;

if(size % 2 < 1E-12) {median =  ((*vec)[size/2] + (*vec)[size/2 + 1])/2;} //odd number
else {median = (*vec)[(size-1)/2];} //even number


return median;

}



//___________________________________________________________
