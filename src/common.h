#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "kseq.h"
#include "uthash.h"

/*!!!!!!!! this must be defined in just one .h file **!!!!!*/
KSEQ_INIT(gzFile, gzread);

struct MPM
{
	char* READ_NAME;
	char* EXON_NAME;	
	int READ_POS;
	int EXON_POS;
	int LENGTH;
};

size_t set_str_arr(char **, char **, size_t);
char* rev_com(char *s);
char* small_dna_str(char*, char*);
int score_DNA_str(char *s);
char* concat(char*, char*);
char* strToUpper(char*);
int max_of_int_array(const int *, size_t);

#endif