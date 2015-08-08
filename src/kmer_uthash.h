
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

struct kmer_uthash {
    char* kmer;                /* key */
	size_t count;
	char **seq_names;    
    UT_hash_handle hh;         /* makes this structure hashable */
};



/*
 * delete duplicate seq_names 
 */
static inline void 
kmer_uthash_uniq(struct kmer_uthash **tb){
	if(*tb == NULL) die("kmer_uthash_uniq: input error");
	struct kmer_uthash *cur, *tmp;
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
static void kmer_uthash_insert(struct kmer_uthash **table, char* kmer, char* name) {
	// check input param
	if(kmer==NULL || name==NULL) die("kmer_uthash_insert: input error");
	struct kmer_uthash *s;
	/* check if kmer exists in table*/
	HASH_FIND_STR(*table, kmer, s);  
	if (s==NULL){
		s = mycalloc(1, struct kmer_uthash);
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
kmer_uthash_write(struct kmer_uthash *htable, char *fname){
	if(htable == NULL || fname == NULL) die("kmer_uthash_write: input error");
	/* write htable to disk*/
	FILE *ofp = fopen(fname, "w");
	if (ofp == NULL) die("Can't open output file %s!\n", fname);
	struct kmer_uthash *s, *tmp;
	HASH_ITER(hh, htable, s, tmp) {
		if(s == NULL) die("Fail to write down %s!\n", fname);
		fprintf(ofp, ">%s\t%zu\n", s->kmer, s->count);		
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
kmer_uthash_destroy(struct kmer_uthash **tb) {
	if(*tb == NULL) die("kmer_uthash_destroy: parameter error\n");	
	/*free the kmer_hash table*/
	struct kmer_uthash *cur, *tmp;
	HASH_ITER(hh, *tb, cur, tmp) {
		if(cur == NULL) die("kmer_uthash_destroy: HASH_ITER fails\n");
		HASH_DEL(*tb, cur);  /* delete it (users advances to next) */
      	free(cur);            /* free it */
    }
	return KM_ERR_NONE;
}

static inline struct kmer_uthash
*find_kmer(struct kmer_uthash *tb, char* quary_kmer) {
	struct kmer_uthash *s = NULL;
    if(tb == NULL || quary_kmer == NULL) die("[%s] parameter error", __func__);
	HASH_FIND_STR(tb, quary_kmer, s);  /* s: output pointer */
    return s;
}

static inline int
kmer_uthash_display(struct kmer_uthash *_kmer_ht) {	
	if(_kmer_ht == NULL) die("kmer_uthash_display: input error\n");
   	struct kmer_uthash *cur, *tmp;
	HASH_ITER(hh, _kmer_ht, cur, tmp) {
		if(cur == NULL) die("kmer_uthash_display: HASH_ITER fails\n");
		printf("kmer=%s\tcount=%zu\n", cur->kmer, cur->count);
		int i;
		for(i=0; i < cur->count; i++){
			printf("%s\t", cur->seq_names[i]);
		}
		printf("\n");
	}
	return KM_ERR_NONE;
}

#endif
