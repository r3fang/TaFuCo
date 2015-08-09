/*--------------------------------------------------------------------*/
/* Created Date: 15JULY2015                                           */
/* Author: Rongxin Fang                                               */
/* Contact: r3fang@ucsd.edu                                           */
/* Library for Breakend Associated Graph (BAG).                       */
/* Functions it contain:                                              */
/*--------------------------------------------------------------------*/

#ifndef _BAG_UTHASH_H
#define _BAG_UTHASH_H

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include "utils.h"
#include "kseq.h"
#include "kmer_uthash.h"
/* error code */
#define BA_ERR_NONE		     0 // no error
/*
 * the BAG_uthash structure
 */
struct BAG_uthash {
	char *edge;
	size_t weight;
	char **read_names;
	char **evidence;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

// junction of gene fusion
typedef struct {
	char* idx; // determined by exon1.exon2.jump_start.jump_end
	char* exon1;
	char* exon2;
	char* s;         // string flanking junction site 
	char *transcript;        // concated exon string 
	int junc_pos;  // junction position on transcript
	int *S1;
	int *S2;
	int S1_num;
	int S2_num;
	size_t hits;         
	double likehood;       // alignment probability
    UT_hash_handle hh;
} junction_t;

static inline junction_t 
*junction_init(int seed_len){
	junction_t *junc = mycalloc(1, junction_t);
	junc->idx = NULL;
	junc->exon1 = NULL;
	junc->exon2 = NULL;
	junc->transcript = NULL;
	junc->junc_pos = 0;
	junc->S1 = NULL;
	junc->S2 = NULL;
	junc->s = mycalloc(seed_len+1, char);	
	junc->S1_num = 0;
	junc->S2_num = 0;	
	junc->hits = 0;	
	junc->likehood = 0.0;	
	return junc;
}

static inline void junction_destory(junction_t **s){
	junction_t *cur, *tmp;
	HASH_ITER(hh, *s, cur, tmp) {
		HASH_DEL(*s, cur);
	}
}

static inline junction_t 
*find_junction(junction_t *jc, char* quary) {
	if(quary == NULL) return NULL;
	junction_t *s;
    HASH_FIND_STR(jc, quary, s);  /* s: output pointer */
	return s;
}

static inline int min_mismatch(char* str, char* pattern){
	if(str == NULL || pattern == NULL) die("[%s] input error"); 
	register int i, j, n;
	int min_mis_match = strlen(pattern)+1;
	char substring[strlen(pattern)+1];
	for(i=0;i<strlen(str)-strlen(pattern);i++){
		n = 0;
		memset(substring, '\0', sizeof(substring));
		strncpy(substring, &str[i], strlen(pattern));
		for(j=0; j<strlen(pattern); j++){if(toupper(pattern[j]) != toupper(substring[j])){n++;}}
		if (n < min_mis_match){min_mis_match = n;} // update min_mis_match
	}
	return min_mis_match;
} 

/*
 * intiate BAG_uthash
 *
 */
static inline struct BAG_uthash 
*BAG_uthash_init() {
	struct BAG_uthash *t = mycalloc(1, struct BAG_uthash);
	t->edge = NULL;
	t->weight = 0;
	t->evidence = mycalloc(1, char*);
	t->read_names = mycalloc(1, char*);
	return t;
}

/*
 * destory 
 */
static inline int 
BAG_uthash_destroy(struct BAG_uthash **table) {
	if(*table == NULL) die("BAG_uthash_destroy: parameter error\n");
	/*free the kmer_hash table*/
	register struct BAG_uthash *cur, *tmp;
	HASH_ITER(hh, *table, cur, tmp) {
		if(cur==NULL) die("BAG_uthash_destroy: HASH_ITER fails\n");
		HASH_DEL(*table, cur); 
		free(cur);   
    }
	return BA_ERR_NONE;
}
/*
 * display BAG_uthash table on the screen
 *
 * PARAMETERS:	struct BAG_uthash *
 * RETURN:	error code
 */
static inline int BAG_uthash_display(struct BAG_uthash *graph_ht) {
	if(graph_ht == NULL) die("BAG_uthash_display: parameter error\n");
	/*free the kmer_hash table*/
	register struct BAG_uthash *cur, *tmp;
	HASH_ITER(hh, graph_ht, cur, tmp) {
		if(cur == NULL) die("BAG_uthash_display: HASH_ITER fails\n");
		int i; for(i=0; i < cur->weight; i++){
			printf(">%s\t%s\n%s\n", cur->edge, cur->read_names[i], cur->evidence[i]);
		}
	}
	return BA_ERR_NONE;
}

static inline struct BAG_uthash
*find_edge(struct BAG_uthash *tb, char* quary_name) {
	if(quary_name == NULL) die("[%s] input error", __func__);
	struct BAG_uthash* s = NULL;	
    HASH_FIND_STR(tb, quary_name, s);  /* s: output pointer */
	return s;
}

/*
 * add one edge to graph
 */
static inline int 
BAG_uthash_add(struct BAG_uthash** graph_ht, char* edge_name, char* read_name, char* evidence){
	/* check parameters */
	if(edge_name == NULL || evidence == NULL) die("[%s]: parameter error\n", __func__);
	struct BAG_uthash *s;
	register int n;
	if((s = find_edge(*graph_ht, edge_name)) == NULL){
		s = BAG_uthash_init();
		s->edge = strdup(edge_name);
		s->weight = 1;
		s->read_names[0] = strdup(read_name);
		s->evidence[0] = strdup(evidence); /* first and only 1 element*/
		HASH_ADD_STR(*graph_ht, edge, s);								
	}else{
		s->weight++;
		char **tmp = mycalloc(s->weight, char*);
	 	for (n = 0; n < s->weight-1; n++){
	 		tmp[n] = strdup(s->evidence[n]);
	 	}
	 	tmp[n] = strdup(evidence);
	 	free(s->evidence);
	 	s->evidence = tmp;
		
		tmp = mycalloc(s->weight, char*);
	 	for (n = 0; n < s->weight-1; n++){
	 		tmp[n] = strdup(s->read_names[n]);
	 	}
	 	tmp[n] = strdup(read_name);
	 	free(s->read_names);
	 	s->read_names = tmp;
	}
	return BA_ERR_NONE;
}

static inline int 
cmpstr(void const *a, void const *b) { 
    char const *aa = (char const *)a;
    char const *bb = (char const *)b;
    return strcmp(aa, bb);
}

/*
 * trim edges with evidence less than min_weight
 */
static inline int 
BAG_uthash_trim(struct BAG_uthash** tb, int min_weight){
	if(*tb == NULL || min_weight < 0) die("BAG_uthash_tim: wrong parameters\n");
	register struct BAG_uthash *cur, *tmp;
	HASH_ITER(hh, *tb, cur, tmp) {
		if(cur==NULL) die("BAG_uthash_tim: HASH_ITER fails\n");
		if(cur->weight < min_weight) {HASH_DEL(*tb, cur); free(cur);}
    }
	return BA_ERR_NONE;
}

/*
 * rm duplicate evidence for graph edge
 */
static inline int 
BAG_uthash_uniq(struct BAG_uthash **tb){
	if(*tb == NULL) die("BAG_uthash_uniq: input error");
	struct BAG_uthash *cur, *tmp;
	int i, j;
	int before, del;
	HASH_ITER(hh, *tb, cur, tmp) {
		if(cur == NULL) die("kmer_uthash_uniq: HASH_ITER fails\n");
		qsort(cur->evidence, cur->weight, sizeof(char*), mystrcmp);
		before = cur->weight;
		del = 0;
		if(cur->weight > 1){
			for(i=1; i<cur->weight; i++){
				if(mystrcmp(cur->evidence[i], cur->evidence[i-1]) == 0){
					cur->evidence[i-1] = NULL;
					del++;
				}
			}
		}		
		if(del > 0){ // if any element has been deleted
			char** tmp = mycalloc(before - del , char*);
			j=0; for(i=0; i<cur->weight; i++){
				if(cur->evidence[i] != NULL){
					tmp[j++] = strdup(cur->evidence[i]);
				}
			}
			free(cur->evidence);
			cur->evidence = tmp;		
			cur->weight = before - del;	
		}
    }	
	return BA_ERR_NONE;
}

/* 
 * Find all genes uniquely matched with kmers on _read.          
 * hash     - a hash table count number of matches between _read and every gene
 * _read    - inqury read
 * _k       - kmer length
 */
static inline int
find_all_genes(str_ctr **hash, struct kmer_uthash *KMER_HT, char* _read, int _k){
/*--------------------------------------------------------------------*/
	/* check parameters */
	if(_read == NULL || _k < 0) die("find_all_MEKMs: parameter error\n");
/*--------------------------------------------------------------------*/
	/* declare vaiables */
	str_ctr *s;
	int _read_pos = 0;
	int num;
	char* gene = NULL;
	register struct kmer_uthash *s_kmer = NULL; 
	char buff[_k];
/*--------------------------------------------------------------------*/
	while(_read_pos<(strlen(_read)-_k+1)){
		/* copy a kmer of string */
		strncpy(buff, _read + _read_pos, _k); buff[_k] = '\0';	
		if(strlen(buff) != _k) die("find_next_match: buff strncpy fails\n");
		/*------------------------------------------------------------*/
		if((s_kmer=find_kmer(KMER_HT, buff)) == NULL){_read_pos++; continue;} // kmer not in table but not an error
		if(s_kmer->count == 1){ // only count the uniq match 
			//gene = strdup(s_kmer->seq_names[0]);
			gene = strsplit(s_kmer->seq_names[0], '.', &num)[0];
			if(gene == NULL) die("find_next_match: get_exon_name fails\n");
			if(str_ctr_add(hash, gene) != 0) die("find_all_MEKMs: str_ctr_add fails\n");
		}
		_read_pos++;
	}
	return 0;
}


/* 
 * Find all exons uniquely matched with kmers on _read.          
 * hash     - a hash table count number of matches between _read and every gene
 * _read    - inqury read
 * _k       - kmer length
 */
static inline int
find_all_exons(str_ctr **hash, struct kmer_uthash *KMER_HT, char* _read, int _k){
/*--------------------------------------------------------------------*/
	/* check parameters */
	if(_read == NULL || _k < 0) die("find_all_MEKMs: parameter error\n");
/*--------------------------------------------------------------------*/
	/* declare vaiables */
	str_ctr *s;
	int _read_pos = 0;
	char* exon = NULL;
	register struct kmer_uthash *s_kmer = NULL; 
	char buff[_k];
/*--------------------------------------------------------------------*/
	while(_read_pos<(strlen(_read)-_k+1)){
		/* copy a kmer of string */
		strncpy(buff, _read + _read_pos, _k); buff[_k] = '\0';	
		if(strlen(buff) != _k) die("find_next_match: buff strncpy fails\n");
		/*------------------------------------------------------------*/
		if((s_kmer=find_kmer(KMER_HT, buff)) == NULL){_read_pos++; continue;} // kmer not in table but not an error
		if(s_kmer->count == 1){ // only count the uniq match 
			exon = strdup(s_kmer->seq_names[0]);
			if(exon == NULL) die("find_next_match: get_exon_name fails\n");
			if(str_ctr_add(hash, exon) != 0) die("find_all_MEKMs: str_ctr_add fails\n");
		}
		_read_pos++;
	}
	return 0;
}

#endif
