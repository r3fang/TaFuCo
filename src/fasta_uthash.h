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

typedef struct{
    char *name;                /* key */
	char *chrom;
	int idx;
	int start;
	int end;
	int l;
	char *seq;
	char *strand; 
	int transcript_num;
	int tss_num;
	char **transcript_id;
	char **tss_id;
	char *gene_id;
	char *gene_name;
    UT_hash_handle hh;         /* makes this structure hashable */
}fasta_t;

/* initilize fasta_t object */
static inline fasta_t
*fasta_init(){
	fasta_t *res = mycalloc(1, fasta_t);
	res->name  = res->chrom = res->seq = res->strand = res->gene_name = res->gene_id = NULL;
	res->transcript_id = res->tss_id = NULL;
	res->tss_num = res->transcript_num = 0;
	return res;
}

static inline int 
fasta_destroy(fasta_t **tb) {
	if(*tb == NULL) die("fasta_uthash_destroy: parameter error\n");
	fasta_t *cur, *tmp;
	HASH_ITER(hh, *tb, cur, tmp) {
		if(cur == NULL) die("fasta_uthash_destroy: HASH_ITER fails\n");
		HASH_DEL(*tb, cur);  /* delete it (users advances to next) */
		if(cur->seq) free(cur->seq);            /* free it */
		if(cur->chrom) free(cur->chrom);            /* free it */
		if(cur->strand) free(cur->strand);            /* free it */
		if(cur->gene_id) free(cur->gene_id);            /* free it */
		if(cur->gene_name) free(cur->gene_name);            /* free it */
		if(cur)      free(cur);                 /* free it */
	}
	return FA_ERR_NONE;
}

static inline int
fasta_display(fasta_t *tb) {
   	fasta_t *cur, *tmp;
	if(tb == NULL) die("fasta_uthash_display: parameter error\n");	
	HASH_ITER(hh, tb, cur, tmp) {
		if(cur == NULL) die("fasta_uthash_display: fail to iterate uthash table\n");
		printf(">%s\n%s\n", cur->name, cur->seq);
	}	
	return FA_ERR_NONE;
}

static inline int
fasta_write(fasta_t *tb, char* fname) {
	if(tb==NULL || fname==NULL) return -1;
	FILE *fp = fopen(fname, "w");
	if(fp==NULL) die("[%s] can't open %s", __func__, fname);
   	fasta_t *cur, *tmp;
	HASH_ITER(hh, tb, cur, tmp) {
		fprintf(fp, ">%s\n%s\n", cur->name, cur->seq);
	}	
	fclose(fp);
	return 0;
}

static inline fasta_t
*fasta_read(char *fname){
	if(fname == NULL) die("[%s] input file name can't be NULL", __func__);
	fasta_t *tb = NULL;
	gzFile fp;
	kseq_t *seq;
	int l;
	int error;
	fp = gzopen(fname, "r");
	if(fp == NULL) die("[%s] fail to open %s\n", __func__, fname);		

	fasta_t *s;	
	if((seq = kseq_init(fp))==NULL) die("[%s]: kseq_init fails\n", __func__);

	while ((l = kseq_read(seq)) >= 0){
		if(seq->name.s == NULL || seq->seq.s==NULL)
			continue;
		s = mycalloc(1, fasta_t);
		s->name = strdup(seq->name.s);
		s->seq = strToUpper(seq->seq.s);
		HASH_ADD_STR(tb, name, s);
	}
	if(seq) kseq_destroy(seq);
	gzclose(fp);
	return tb;
}


static inline fasta_t 
*find_fasta(fasta_t *tb, char* quary_name) {
	if(quary_name == NULL) die("[%s] input error", __func__);
	fasta_t* s = NULL;	
    HASH_FIND_STR(tb, quary_name, s);  /* s: output pointer */
	return s;
}

static inline int 
name_sort(fasta_t *a, fasta_t *b) {
    return strcmp(a->name, b->name);
}

#endif
