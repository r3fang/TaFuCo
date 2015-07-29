#ifndef KMER_UTHASH_H
#define KMER_UTHASH_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "common.h"
#include "kseq.h"
#include "uthash.h"
#include "utils.h"

#define KM_ERR_NONE					0

struct kmer_uthash {
    char* kmer;                /* key */
	size_t count;
	char **seq_names;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

static inline int mystrcmp(const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

static inline int sortstring( const char *str1, const char *str2 )
{
    int result = 0;
    int val = strcmp(str1, str2);     
    if ( val < 0 ) result = -1;
    if ( val > 0 ) result = 1;
    return result;
}

static inline void kmer_uthash_uniq(struct kmer_uthash **tb){
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

/* split string*/
static inline char** strsplit(char* s, const char delim){
	if(s==NULL) die("strsplit: input error");
	kstring_t *ks = mycalloc(1, kstring_t);
	ks->s = strdup(s);
	ks->l = strlen(s);
	int *fields, n, i;
	fields = ksplit(ks, delim, &n);
	if(n==0) return NULL;
	char** ret = mycalloc(n, char*);
	for(i=0; i<n; i++) ret[i] = strdup(ks->s + fields[i]);
	return ret;
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

static inline struct kmer_uthash 
*kmer_uthash_construct(char *fasta_file, int k){
	gzFile fp;  
	kseq_t *seqs;  
	int l;
	char *kmer = mycalloc(k+1, char);	
	struct kmer_uthash *table = NULL;
	char *seq,  *name;
	seq = name = NULL;
	fp = gzopen(fasta_file, "r");
	if (fp == NULL) die("Can't open %s\n", fasta_file);
	seqs = kseq_init(fp);	
	if (seqs == NULL) die("kseq_init fails\n");
	while ((l = kseq_read(seqs)) >= 0) {
		seq = strToUpper(seqs->seq.s);
		name = strdup(seqs->name.s);		
		if(seq == NULL || name == NULL || strlen(seq) <= k){
			continue;
		}
		int i; for(i=0; i < strlen(seq)-k+1; i++){
			memset(kmer, '\0', k+1);
			strncpy(kmer, seq+i, k);
			kmer_uthash_insert(&table, kmer, strsplit(name, '.')[0]); 
		}
	}
	if(kmer) free(kmer);  	
	if(seq) free(seq);
	if(name) free(name);
	kseq_destroy(seqs);
	gzclose(fp);
	kmer_uthash_uniq(&table);
	return table;
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

static inline int
find_kmer(struct kmer_uthash *tb, char* quary_kmer, struct kmer_uthash **s) {
	*s = NULL;
    if(tb == NULL || quary_kmer == NULL || *s != NULL) die("find_kmer: parameter error\n");
	HASH_FIND_STR(tb, quary_kmer, *s);  /* s: output pointer */
    return KM_ERR_NONE;
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

//static inline int
//kmer_uthash_load(char *fname, int *k, struct kmer_uthash** htable){	
//	if(fname == NULL || *htable != NULL) die("kmer_uthash_load: parameter error\n");
//	gzFile fp;  kseq_t *seq; int l;
//	
//	if ((fp = gzopen(fname, "r")) == NULL) die("Can't open input file %s!\n", fname);
//	if ((seq = kseq_init(fp)) == NULL) die("kmer_uthash_load: kseq_init fails\n"); // STEP 3: initialize seq  
//	
//	while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence 
//		/* add kmer */
//		char *kmer = seq->name.s;
//		*k = strlen(kmer);
//		struct kmer_uthash *s;
//		if((s = malloc(1 * sizeof(struct kmer_uthash))) == NULL) die("kmer_uthash_load: malloc fails\n");
//		strncpy(s->kmer, kmer, *k);
//		s->kmer[*k] = '\0'; /* just in case*/
//		if(strlen(s->kmer) != *k) die("kmer_uthash_load: strncpy fails\n");
//		/* add count */
//		int count = atoi(seq->comment.s);
//		s->count = count;		
//		/* add pos */
//		char *pos = seq->seq.s;				
//		s->seq_names = malloc((s->count) * sizeof(char*));
//		/* split a string by delim */
//		char *token;	    
//		/* get the first token */
//	    token = strtok(seq_names, "|");				
//		/* walk through other tokens */
//		int i = 0;
//	    while(token != NULL) 
//	    {
//			s->pos[i] = malloc((strlen(token)+1) * sizeof(char));
//			/*duplicate a string*/
//			s->pos[i] = strdup(token);
//			token = strtok(NULL, "|");
//			i ++;
//	    }		
//		HASH_ADD_STR(*htable, kmer, s);
//	}	
//	gzclose(fp);
//	kseq_destroy(seq);
//	return KM_ERR_NONE;
//}

#endif
