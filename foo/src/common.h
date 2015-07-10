#ifndef COMMON_H
#define COMMON_H


#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include "uthash.h"
#include "utlist.h"
#include "utstring.h"
#include "utarray.h"
#include "kseq.h"

#define MAX_K 100

struct kmer_uthash {
    char kmer[MAX_K];                /* key */
	char *pos;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

void kmer_table_destroy(struct kmer_uthash **table);
char* concat(char *s1, char *s2);
char* strToUpper(char* s);

#endif