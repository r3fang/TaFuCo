#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include <assert.h>
#include "uthash.h"
#include "kseq.h"
#include "common.h"
#include "kmer_uthash.h"

KSEQ_INIT(gzFile, gzread);

struct kmer_uthash*
find_kmer(char* quary_kmer, struct kmer_uthash *tb) {
    struct kmer_uthash *s;
    HASH_FIND_STR(tb, quary_kmer, s);  /* s: output pointer */
    return s;
}

void 
kmer_uthash_destroy(struct kmer_uthash **table) {
	/*free the kmer_hash table*/
  struct kmer_uthash *cur, *tmp;
  HASH_ITER(hh, *table, cur, tmp) {
      HASH_DEL(*table, cur);  /* delete it (users advances to next) */
      free(cur);            /* free it */
    }
}

void
kmer_uthash_display(struct kmer_uthash *_kmer_ht) {	
   	struct kmer_uthash *cur, *tmp;
	HASH_ITER(hh, _kmer_ht, cur, tmp) {
		if(cur == NULL)
			exit(-1);
		printf("kmer=%s\tcount=%d\n", cur->kmer, cur->count);
		int i;
		for(i=0; i < cur->count; i++){
			printf("%s\t", cur->pos[i]);
		}
		printf("\n");
	}	
}

struct kmer_uthash*
kmer_uthash_load(char *fname, int *k){	
	/* load kmer_htable*/
	struct kmer_uthash *htable = NULL;
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
		/* add kmer */
		char *kmer = seq->name.s;
		*k = strlen(kmer);
		struct kmer_uthash *s;
		s = (struct kmer_uthash*)malloc(sizeof(struct kmer_uthash));
		strncpy(s->kmer, kmer, *k);
		s->kmer[*k] = '\0'; /* just in case*/
		
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
		HASH_ADD_STR(htable, kmer, s);
	}	
	gzclose(fp); // STEP 6: close the file handler  
	kseq_destroy(seq);
	return(htable);
}

	
	
	