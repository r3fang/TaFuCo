#ifndef FASTA_UTHASH_H
#define FASTA_UTHASH_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include "uthash.h"
#include "kseq.h"
#include "utils.h"

/* error code */
#define FA_ERR_NONE		     0 // no error

struct fasta_uthash {
    char *name;                /* key */
	char *chrom;
	int start;
	int end;
	int l;
	char *seq;
    UT_hash_handle hh;         /* makes this structure hashable */
};

static inline int 
fasta_uthash_destroy(struct fasta_uthash **tb) {
	if(*tb == NULL) die("fasta_uthash_destroy: parameter error\n");
	struct fasta_uthash *cur, *tmp;
	HASH_ITER(hh, *tb, cur, tmp) {
		if(cur == NULL) die("fasta_uthash_destroy: HASH_ITER fails\n");
		HASH_DEL(*tb, cur);  /* delete it (users advances to next) */
		free(cur);            /* free it */
	}
	return FA_ERR_NONE;
}

static inline int
fasta_uthash_display(struct fasta_uthash *tb) {
   	struct fasta_uthash *cur, *tmp;
	if(tb == NULL) die("fasta_uthash_display: parameter error\n");	
	HASH_ITER(hh, tb, cur, tmp) {
		if(cur == NULL) die("fasta_uthash_display: fail to iterate uthash table\n");
		printf(">%s\n%s\n", cur->name, cur->seq);
	}	
	return FA_ERR_NONE;
}

static inline int
fasta_uthash_write(struct fasta_uthash *tb, char* fname) {
	if(tb==NULL || fname==NULL) return -1;
	FILE *fp = fopen(fname, "w");
	if(fp==NULL) die("[%s] can't open %s", __func__, fname);
   	struct fasta_uthash *cur, *tmp;
	HASH_ITER(hh, tb, cur, tmp) {
		fprintf(fp, ">%s\n%s\n", cur->name, cur->seq);
	}	
	fclose(fp);
	return 0;
}

static inline struct fasta_uthash
*fasta_uthash_load(char *fname){
	if(fname == NULL) die("fasta_uthash_load: parameter error\n"); 
	struct fasta_uthash *tb = NULL;
	gzFile fp;
	kseq_t *seq;
	int l;
	int error;
	fp = gzopen(fname, "r");
	if(fp == NULL) die("fasta_uthash_load: fail to open %s\n", fname);		

	struct fasta_uthash *s;	
	if((seq = kseq_init(fp))==NULL) die("fasta_uthash_load: kseq_init fails\n");

	while ((l = kseq_read(seq)) >= 0){
		if((s = malloc(sizeof(struct fasta_uthash))) == NULL) die("fasta_uthash_load: fail to malloc");
		if(seq->name.s == NULL || seq->seq.s==NULL)
			continue;
		s->name = strdup(seq->name.s);
		s->seq = strToUpper(seq->seq.s);
		HASH_ADD_STR(tb, name, s);
	}
	if(seq) kseq_destroy(seq);
	gzclose(fp);
	return tb;
}


static inline struct fasta_uthash 
*find_fasta(struct fasta_uthash *tb, char* quary_name) {
	if(quary_name == NULL) die("[%s] input error", __func__);
	struct fasta_uthash* s = NULL;	
    HASH_FIND_STR(tb, quary_name, s);  /* s: output pointer */
	return s;
}
#endif
