
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
#define MAX_KMER_LEN                40

typedef struct{
    char* kmer;                /* key */
	int count;
	char **seq_names;    
    UT_hash_handle hh;         /* makes this structure hashable */
} kmer_t;


static inline kmer_t *find_kmer(kmer_t*, char*);
static inline void kmer_uniq(kmer_t**);
static void kmer_add(kmer_t**, char*, char*);
static inline void kmer_write(kmer_t*, char*);
static inline int kmer_destroy(kmer_t**);
static inline void kmer_display(kmer_t*);

static inline kmer_t 
*kmer_init(){
	kmer_t *res = mycalloc(1, kmer_t);
	res->kmer = NULL;
	res->count = 0;
	res->seq_names = NULL;
	return res;
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
    if(quary_kmer == NULL) die("[%s] parameter error", __func__);
	kmer_t *s = NULL;
	HASH_FIND_STR(tb, quary_kmer, s);  /* s: output pointer */
    return s;
}

static inline void
kmer_display(kmer_t *kmer_ht) {	
	if(kmer_ht == NULL) die("kmer_uthash_display: input error\n");
   	register kmer_t *kmer_cur;
	register int i;
	
	for(kmer_cur=kmer_ht; kmer_cur!=NULL; kmer_cur=kmer_cur->hh.next){
		printf("kmer=%s\tcount=%d\n", kmer_cur->kmer, kmer_cur->count);
		for(i=0; i < kmer_cur->count; i++){
			printf("%s\t", kmer_cur->seq_names[i]);
		}		
		printf("\n");
	}
}

/* add one kmer and its exon name to kmer_uthash table */
static void kmer_add(kmer_t **table, char* kmer, char* name) {
	// check input param
	if(kmer==NULL || name==NULL) die("[%s]: input error", __func__);

	register kmer_t *s;
	register int i;
	/* check if kmer exists in table*/
	if((s = find_kmer(*table, kmer))==NULL){
		s = kmer_init();
		s->kmer = strdup(kmer);
		s->count = 1;                /* first pos in the list */
		s->seq_names = mycalloc(s->count, char*);
		s->seq_names[0] = strdup(name); /* first and only 1 element*/
		HASH_ADD_STR(*table, kmer, s); // add to hash table
	}else{
		char **tmp;
		s->count++;
		/* copy s->pos */
		tmp = mycalloc(s->count, char*);
		for(i = 0; i < s->count-1; i++){
			tmp[i] = strdup(s->seq_names[i]);
		}
		tmp[i] = strdup(name);
		free(s->seq_names);
		s->seq_names = tmp;
	}
}

/*
 * delete duplicate seq_names in kmer->seq_names
 */
static inline void 
kmer_uniq(kmer_t **kmer_ht){
	if(*kmer_ht == NULL) die("[%s] input can't be NULL", __func__);
	register kmer_t *kmer_cur;
	char **names = NULL;
	bool repeat;
	register int i, j, count;
	for(kmer_cur=*kmer_ht; kmer_cur!=NULL; kmer_cur=kmer_cur->hh.next){
		names = mycalloc(kmer_cur->count, char*);
		count = 0;
		for(i=0; i<kmer_cur->count; i++){ /* iterate every evidence */
			repeat = false;
			for(j=0; j<count; j++){if(strcmp(kmer_cur->seq_names[i], names[j])==0){repeat = true;}}
			if(repeat==false){ // no duplicates
				names[count++] = strdup(kmer_cur->seq_names[i]);
			}
		}
		if(kmer_cur->seq_names){for(i=0; i<kmer_cur->count; i++) free(kmer_cur->seq_names[i]); free(kmer_cur->seq_names);};
		kmer_cur->seq_names = names;
		kmer_cur->count = count;
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

#endif
