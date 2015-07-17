#ifndef FASTA_UTHASH_H
#define FASTA_UTHASH_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include "uthash.h"
#include "kseq.h"
#include "common.h"
#include "utils.h"

/* error code */
#define FA_ERR_NONE		     0 // no error
#define FA_ERR_PARAM		-1 // bad paraemter
#define FA_ERR_HASHITER		-2 // fail to iterate uthash
#define FA_ERR_MALLOC		-3 // fail to malloc memory uthash
#define FA_ERR_FILE 		-4 // fail to malloc memory uthash

struct fasta_uthash {
    char* name;                /* key */
	char* seq;
    UT_hash_handle hh;         /* makes this structure hashable */
};

static inline int 
fasta_uthash_load(char *fname, struct fasta_uthash **tb){
	gzFile fp;
	kseq_t *seq;
	int l;
	int error;
	if(*tb != NULL || fname == NULL) goto FAIL_PARAM; 
	fp = gzopen(fname, "r");
	if(fp == NULL) goto FAIL_FILE;		
	struct fasta_uthash *s;	
	if((seq = kseq_init(fp))==NULL) goto FAIL_MALLOC;

	while ((l = kseq_read(seq)) >= 0){
		if((s = malloc(sizeof(struct fasta_uthash))) == NULL) goto FAIL_MALLOC;
		if(seq->name.s == NULL || seq->seq.s==NULL)
			continue;
		s->name = strdup(seq->name.s);
		s->seq = strToUpper(seq->seq.s);
		HASH_ADD_STR(*tb, name, s);
	}	
	goto SUCCESS;
	
	FAIL_PARAM:
		error = FA_ERR_PARAM;
		return error;
	FAIL_FILE:
		error = FA_ERR_FILE;
		return error;
	FAIL_MALLOC:
		error = FA_ERR_MALLOC;
		goto EXIT;
	SUCCESS:
		error = FA_ERR_NONE;
		goto EXIT;
	EXIT:
		if(seq) kseq_destroy(seq);
		gzclose(fp);
		return error;
}

static inline void
fasta_uthash_display(struct fasta_uthash *_fasta_ht) {	
   	struct fasta_uthash *cur, *tmp;
	HASH_ITER(hh, _fasta_ht, cur, tmp) {
		if(cur == NULL)
			exit(-1);
		printf(">%s\n%s\n", cur->name, cur->seq);
	}	
}

static inline void 
fasta_uthash_destroy(struct fasta_uthash **table) {
	/*free the kmer_hash table*/
  struct fasta_uthash *cur, *tmp;
  HASH_ITER(hh, *table, cur, tmp) {
      HASH_DEL(*table, cur);  /* delete it (users advances to next) */
      free(cur);            /* free it */
    }
}

static inline struct fasta_uthash*
find_fasta(char* quary_name, struct fasta_uthash *tb) {
    struct fasta_uthash *s;
    HASH_FIND_STR(tb, quary_name, s);  /* s: output pointer */
    return s;
}

#endif
