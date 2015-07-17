#ifndef KMER_UTHASH_H
#define KMER_UTHASH_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "common.h"
#include "kseq.h" 
#include "uthash.h"

#define MAX_K 100

struct kmer_uthash {
    char kmer[MAX_K];                /* key */
	int count;
	char **pos;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

struct kmer_uthash *kmer_uthash_load(char*, int*);
struct kmer_uthash *find_kmer(char* , struct kmer_uthash*);
void kmer_uthash_display(struct kmer_uthash *);	
void kmer_uthash_destroy(struct kmer_uthash **);

#endif
