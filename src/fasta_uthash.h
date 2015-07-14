#ifndef FASTA_UTHASH_H
#define FASTA_UTHASH_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "kseq.h" 
#include "uthash.h"

struct fasta_uthash {
    char* name;                /* key */
	char* seq;
    UT_hash_handle hh;         /* makes this structure hashable */
};

struct fasta_uthash *fasta_uthash_load(char *);
struct fasta_uthash *find_fasta(char* , struct fasta_uthash*);
void fasta_uthash_display(struct fasta_uthash*);
void fasta_uthash_destroy(struct fasta_uthash *table);

#endif
