#ifndef __export_plink_H__
#define __export_plink_H__

#include <vector.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif

void get_snps_many(char *a, int *Nsnps, int *Nrows, int *b);

std::string* getGenotype(std::string coding, std::string sep);

SEXP export_plink(SEXP idnames, SEXP snpdata, SEXP Nsnps, SEXP NidsTotal, SEXP Coding, SEXP from, SEXP to,
		SEXP male, SEXP traits, SEXP pedfilename, SEXP plink, SEXP append);

#ifdef __cplusplus
}
#endif


#endif

