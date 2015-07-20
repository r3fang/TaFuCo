#ifndef KMER_UTHASH_H
#define KMER_UTHASH_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "common.h"
#include "htslib/kseq.h"
#include "uthash.h"
#include "utils.h"

#define MAX_K 						100
#define KM_ERR_NONE					0

struct kmer_uthash {
    char kmer[MAX_K];                /* key */
	int count;
	char **pos;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

static inline int
find_kmer(struct kmer_uthash *tb, char* quary_kmer, struct kmer_uthash **s) {
	*s = NULL;
    if(tb == NULL || quary_kmer == NULL || *s != NULL) die("find_kmer: parameter error\n");
	HASH_FIND_STR(tb, quary_kmer, *s);  /* s: output pointer */
    return KM_ERR_NONE;
}

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

static inline int
kmer_uthash_display(struct kmer_uthash *_kmer_ht) {	
	if(_kmer_ht == NULL) die("kmer_uthash_display: parameter error\n");
   	struct kmer_uthash *cur, *tmp;
	HASH_ITER(hh, _kmer_ht, cur, tmp) {
		if(cur == NULL) die("kmer_uthash_display: HASH_ITER fails\n");
		printf("kmer=%s\tcount=%d\n", cur->kmer, cur->count);
		int i;
		for(i=0; i < cur->count; i++){
			printf("%s\t", cur->pos[i]);
		}
		printf("\n");
	}
	return KM_ERR_NONE;
}

static inline int
kmer_uthash_load(char *fname, int *k, struct kmer_uthash** htable){	
	if(fname == NULL || *htable != NULL) die("kmer_uthash_load: parameter error\n");
	gzFile fp;  kseq_t *seq; int l;
	
	if ((fp = gzopen(fname, "r")) == NULL) die("Can't open input file %s!\n", fname);
	if ((seq = kseq_init(fp)) == NULL) die("kmer_uthash_load: kseq_init fails\n"); // STEP 3: initialize seq  
	
	while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence 
		/* add kmer */
		char *kmer = seq->name.s;
		*k = strlen(kmer);
		struct kmer_uthash *s;
		if((s = malloc(1 * sizeof(struct kmer_uthash))) == NULL) die("kmer_uthash_load: malloc fails\n");
		strncpy(s->kmer, kmer, *k);
		s->kmer[*k] = '\0'; /* just in case*/
		if(strlen(s->kmer) != *k) die("kmer_uthash_load: strncpy fails\n");
		/* add count */
		int count = atoi(seq->comment.s);
		s->count = count;		
		/* add pos */
		char *pos = seq->seq.s;				
		s->pos = malloc((s->count) * sizeof(char*));
		/* split a string by delim */
		char *token;	    
		/* get the first token */
	    token = strtok(pos, "|");				
		/* walk through other tokens */
		int i = 0;
	    while(token != NULL) 
	    {
			s->pos[i] = malloc((strlen(token)+1) * sizeof(char));
			/*duplicate a string*/
			s->pos[i] = strdup(token);
			token = strtok(NULL, "|");
			i ++;
	    }		
		HASH_ADD_STR(*htable, kmer, s);
	}	
	gzclose(fp);
	kseq_destroy(seq);
	return KM_ERR_NONE;
}

#endif
