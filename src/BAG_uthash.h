/*--------------------------------------------------------------------*/
/* Created Date: 15JULY2015                                           */
/* Author: Rongxin Fang                                               */
/* Contact: r3fang@ucsd.edu                                           */
/* Library for Breakend Associated Graph (BAG).                       */
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

/* error code */
#define BA_ERR_NONE		     0 // no error

/*
 * the BAG_uthash structure
 */
struct BAG_uthash {
	char *edge;
	size_t weight;
	char **evidence;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

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
		printf(">%s\t%zu\n", cur->edge, cur->weight);
		int i; for(i=0; i<cur->weight; i++){
			printf("%s\n", cur->evidence[i]);	
		}
	}
	return BA_ERR_NONE;
}

/*
 * add one edge to graph
 */
static inline int 
BAG_uthash_add(struct BAG_uthash** graph_ht, char* edge_name, char* evidence){
	/* check parameters */
	if(edge_name == NULL || evidence == NULL) die("BAG_uthash_add: parameter error\n");
	register struct BAG_uthash *s;
	HASH_FIND_STR(*graph_ht, edge_name, s);
	if(s==NULL){
		if((s=(struct BAG_uthash*)malloc(sizeof(struct BAG_uthash))) == NULL) die("BAG_uthash_add: malloc fails\n");
		s->edge = edge_name;
		s->weight = 1;
		if((s->evidence = malloc(sizeof(char*))) == NULL) die("BAG_uthash_add: malloc fails\n");
		s->evidence[0] = evidence; /* first and only 1 element*/
		HASH_ADD_STR(*graph_ht, edge, s);						
	}else{
		s->weight ++;
		char **tmp;
		if((tmp = malloc((s->weight+1) * sizeof(char*)))==NULL) die("BAG_uthash_add: malloc fails\n");
		int n;
		for (n = 0; n < s->weight-1; n++){
			tmp[n] = strdup(s->evidence[n]);
		}
		free(s->evidence);
		tmp[n] = evidence;
		s->evidence = tmp;
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

///*
// * trim edges with evidence less than min_weight
// */
//static inline struct BAG_uthash *BAG_uthash_load(char* fname){
//	if(fname == NULL) die("[%s] input error", __func__);
//	struct BAG_uthash *tb = NULL;
//	gzFile fp;
//	kseq_t *seq;
//	int l;
//	fp = gzopen(fname, "r");
//	if(fp == NULL) die("[%s] fail to open %s", __func__, fname);		
//
//	struct fasta_uthash *s;	
//	if((seq = kseq_init(fp))==NULL) die("[%s] kseq_init fails", __func__);
//
//	while ((l = kseq_read(seq)) >= 0){
//		if(seq->name.s == NULL || seq->seq.s==NULL)
//			continue;
//		printf("%s\n", seq->name.s);
//		printf("%s\n", seq->seq.s);
//	}	
//	if(seq) kseq_destroy(seq);
//	gzclose(fp);
//	return tb;
//}
#endif
