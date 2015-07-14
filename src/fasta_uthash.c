#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include <assert.h>
#include "uthash.h"
#include "kseq.h"
#include "common.h"
#include "fasta_uthash.h"

KSEQ_INIT(gzFile, gzread);

void
fasta_uthash_display(struct fasta_uthash *_fasta_ht) {	
   	struct fasta_uthash *cur, *tmp;
	HASH_ITER(hh, _fasta_ht, cur, tmp) {
		if(cur == NULL)
			exit(-1);
		printf(">%s\n%s\n", cur->name, cur->seq);
	}	
}

void 
fasta_uthash_destroy(struct fasta_uthash *table) {
	/*free the kmer_hash table*/
  struct fasta_uthash *cur, *tmp;
  HASH_ITER(hh, table, cur, tmp) {
      HASH_DEL(table, cur);  /* delete it (users advances to next) */
      free(cur);            /* free it */
    }
}

struct fasta_uthash*
find_fasta(char* quary_name, struct fasta_uthash *tb) {
    struct fasta_uthash *s;
    HASH_FIND_STR(tb, quary_name, s);  /* s: output pointer */
    return s;
}

struct fasta_uthash* 
fasta_uthash_load(char *fname){
	gzFile fp;
	kseq_t *seq;
	int l;
	struct fasta_uthash *res = NULL;
	fp = gzopen(fname, "r");
	if(fp == NULL){
		return NULL;
	}
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0){
		struct fasta_uthash *s = malloc(sizeof(struct fasta_uthash));
		if(s == NULL)
			continue;
		char *name = seq->name.s;
		char *_seq = seq->seq.s;
		char *_seq_upper = strToUpper(_seq);
		if(name==NULL || _seq==NULL || _seq_upper==NULL)
			continue;
		s->name = strdup(name);
		s->seq = strdup(_seq_upper);
		HASH_ADD_STR(res, name, s);
	}	
	kseq_destroy(seq);
	gzclose(fp);
	return res;
}
