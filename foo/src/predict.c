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
KSEQ_INIT(gzFile, gzread)  

#define k_predict 30
	
struct kmer_uthash_predict {
    char *kmer;                /* key */
	char *name;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

/* Global variables */
struct kmer_uthash_predict *table_predict = NULL;

void kmer_table_destroy_pred() {
  struct kmer_uthash_predict *cur, *tmp;
  HASH_ITER(hh, table_predict, cur, tmp) {
      HASH_DEL(table_predict, cur);  /* delete it (users advances to next) */
      free(cur);            /* free it */
    }
}

int load_index_table(char *fname){
	gzFile fp;  
	kseq_t *seq;  
	int l;
	fp = gzopen(fname, "r");
	if (fp == NULL) {
	  fprintf(stderr, "Can't open input file %s!\n", fname);
	  exit(1);
	}
	seq = kseq_init(fp); // STEP 3: initialize seq  
	if (seq == NULL){
		return NULL;
	}
	while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence 
		char *kmer = seq->name.s;
		char *name = seq->seq.s;
		struct kmer_uthash_predict *s;
		HASH_FIND_STR(table_predict, kmer, s);
		if (s==NULL){
			s = (struct kmer_uthash_predict*)malloc(sizeof(struct kmer_uthash_predict));
			s->kmer = kmer;
			s->name = name;
			HASH_ADD_STR(table_predict, kmer, s);
		}
	}
	gzclose(fp); // STEP 6: close the file handler  
	return 0;
}

int predict_main(char *fasta_file){
	/* main function*/
	char *index_file = concat(fasta_file, ".index");
	load_index_table(index_file);

	struct kmer_uthash_predict *s, *tmp;
	HASH_ITER(hh, table_predict, s, tmp) {
		printf("%s\t%s\n", s->kmer, s->name);						
	}
	

	kmer_table_destroy_pred();
	return 0;
}

