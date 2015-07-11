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
	int count;
	char **pos;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

struct fasta_uthash {
    char* name;                /* key */
	char* seq;
	char* comment;
    UT_hash_handle hh;         /* makes this structure hashable */
};

void fasta_uthash_destroy(struct fasta_uthash**);
void kmer_uthash_destroy(struct kmer_uthash**);
char* concat(char*, char*);
char* strToUpper(char*);
int write_kmer_htable(struct kmer_uthash**, char*);
int strsplit (const char *str, size_t size, char *parts[], const char *delimiter);
int strsplit_size (const char *str, const char *delimiter);
char *pos_parser(char *str, int *i);

#endif