
#ifndef KMER_UTHASH_H
#define KMER_UTHASH_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "kseq.h"
#include "kstring.h"
#include "uthash.h"
#include "utils.h"

#define KM_ERR_NONE					0

typedef struct{
    char* kmer;                /* key */
	int count;
	char **seq_names;    
    UT_hash_handle hh;         /* makes this structure hashable */
} kmer_t;

/*
 * delete duplicate seq_names 
 */
static inline void 
kmer_uniq(kmer_t **tb){
	if(*tb == NULL) die("kmer_uthash_uniq: input error");
	kmer_t *cur, *tmp;
	int i, j;
	int before, del;
	HASH_ITER(hh, *tb, cur, tmp) {
		if(cur == NULL) die("kmer_uthash_uniq: HASH_ITER fails\n");
		qsort(cur->seq_names, cur->count, sizeof(char*), mystrcmp);
		before = cur->count;
		del = 0;
		if(cur->count > 1){
			for(i=1; i<cur->count; i++){
				if(mystrcmp(cur->seq_names[i], cur->seq_names[i-1]) == 0){
					cur->seq_names[i-1] = NULL;
					del++;
				}
			}
		}		
		if(del > 0){ // if any element has been deleted
			char** tmp = mycalloc(before - del , char*);
			j=0; for(i=0; i<cur->count; i++){
				if(cur->seq_names[i] != NULL){
					tmp[j++] = strdup(cur->seq_names[i]);
				}
			}
			free(cur->seq_names);
			cur->seq_names = tmp;		
			cur->count = before - del;	
		}
    }	
}

/* add one kmer and its exon name to kmer_uthash table */
static void kmer_insert(kmer_t **table, char* kmer, char* name) {
	// check input param
	if(kmer==NULL || name==NULL) die("kmer_uthash_insert: input error");
	kmer_t *s;
	/* check if kmer exists in table*/
	HASH_FIND_STR(*table, kmer, s);  
	if (s==NULL){
		s = mycalloc(1, kmer_t);
		s->kmer = strdup(kmer);
		s->count = 1;                /* first pos in the list */
		s->seq_names = mycalloc(s->count, char*);
		s->seq_names[0] = strdup(name); /* first and only 1 element*/
		HASH_ADD_STR(*table, kmer, s); // add to hash table
	}else{
		char **tmp;
		s->count += 1;
		/* copy s->pos */
		tmp = mycalloc(s->count, char*);
		int i; for (i = 0; i < s->count-1; i++){
			tmp[i] = strdup(s->seq_names[i]);
		}
		free(s->seq_names);
		/* append pos */
		tmp[i] = strdup(name);
		/* assign tmp to s->pos*/
		s->seq_names = tmp;
	}
}

///* Write down kmer_uthash */
static inline void 
kmer_write(kmer_t *htable, char *fname){
	if(htable == NULL || fname == NULL) die("kmer_uthash_write: input error");
	/* write htable to disk*/
	FILE *ofp = fopen(fname, "w");
	if (ofp == NULL) die("Can't open output file %s!\n", fname);
	kmer_t *s, *tmp;
	HASH_ITER(hh, htable, s, tmp) {
		if(s == NULL) die("Fail to write down %s!\n", fname);
		fprintf(ofp, ">%s\t%d\n", s->kmer, s->count);		
		int i;
		for(i=0; i < s->count; i++){
			if(i==0){
				fprintf(ofp, "%s", s->seq_names[i]);																
			}else{
				fprintf(ofp, "|%s", s->seq_names[i]);
			}
		}
		fprintf(ofp, "\n");
	}
	fclose(ofp);
}

// destory 
static inline int 
kmer_destroy(kmer_t **tb) {
	if(*tb == NULL) die("kmer_uthash_destroy: parameter error\n");	
	/*free the kmer_hash table*/
	kmer_t *cur, *tmp;
	HASH_ITER(hh, *tb, cur, tmp) {
		if(cur == NULL) die("kmer_uthash_destroy: HASH_ITER fails\n");
		HASH_DEL(*tb, cur);  /* delete it (users advances to next) */
      	free(cur);            /* free it */
    }
	return KM_ERR_NONE;
}

static inline kmer_t
*find_kmer(kmer_t *tb, char* quary_kmer) {
	kmer_t *s = NULL;
    if(tb == NULL || quary_kmer == NULL) die("[%s] parameter error", __func__);
	HASH_FIND_STR(tb, quary_kmer, s);  /* s: output pointer */
    return s;
}

static inline int
kmer_display(kmer_t *_kmer_ht) {	
	if(_kmer_ht == NULL) die("kmer_uthash_display: input error\n");
   	kmer_t *cur, *tmp;
	HASH_ITER(hh, _kmer_ht, cur, tmp) {
		if(cur == NULL) die("kmer_uthash_display: HASH_ITER fails\n");
		printf("kmer=%s\tcount=%d\n", cur->kmer, cur->count);
		int i;
		for(i=0; i < cur->count; i++){
			printf("%s\t", cur->seq_names[i]);
		}
		printf("\n");
	}
	return KM_ERR_NONE;
}

#endif
