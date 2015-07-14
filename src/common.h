#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "kseq.h" 
#include "uthash.h"

struct MPM
{
	char* READ_NAME;
	char* EXON_NAME;	
	int READ_POS;
	int EXON_POS;
	int LENGTH;
};

void MPM_display(struct MPM *);
char* rev_com(char *s);
char* small_dna_str(char*, char*);
int score_DNA_str(char *s);
char* concat(char*, char*);
char* strToUpper(char*);
int strsplit (const char *str, size_t size, char *parts[], const char *delimiter);
int strsplit_size (const char *str, const char *delimiter);
char *pos_parser(char *str, int *i);

#endif