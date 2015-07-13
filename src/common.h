#ifndef COMMON_H
#define COMMON_H


#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h> 
#include "uthash.h"
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

void kmer_uthash_display(struct kmer_uthash*);	
void fasta_uthash_display(struct fasta_uthash*);
void fasta_uthash_destroy(struct fasta_uthash*);
void kmer_uthash_destroy(struct kmer_uthash*);
char* concat(char*, char*);
char* strToUpper(char*);
int strsplit (const char *str, size_t size, char *parts[], const char *delimiter);
int strsplit_size (const char *str, const char *delimiter);
char *pos_parser(char *str, int *i);

#endif