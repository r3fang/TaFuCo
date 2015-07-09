#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include "uthash.h"
#include "utlist.h"
#include "utstring.h"
#include "utarray.h"
#include "kseq.h"
#include "common.h"

/* kmer length */
#define k_index 30

KSEQ_INIT(gzFile, gzread)  
	
struct kmer_uthash {
    char kmer[k_index];                /* key */
	char *pos;    
	int count;
    UT_hash_handle hh;         /* makes this structure hashable */
};

void add_to_kmer_hash(char kmer[k_index], char* pos);
void kmer_table_destroy();
int index_main(char *fasta_file);
