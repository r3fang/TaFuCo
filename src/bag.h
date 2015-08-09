/*--------------------------------------------------------------------*/
/* Created Date: 15JULY2015                                           */
/* Author: Rongxin Fang                                               */
/* Contact: r3fang@ucsd.edu                                           */
/* Library for Breakend Associated Graph (bag).                       */
/*--------------------------------------------------------------------*/

#ifndef _BAG_H
#define _BAG_H

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include "utils.h"
#include "kseq.h"

/* error code */
#define BA_ERR_NONE		     0 // no error

typedef struct {
	char* idx;          /* the key of this junction, determined by exon1.start.exon2.end, must be unique */
	char* exon1;       
	char* exon2;
	char* s;            /* short junction string that flanks the junction sites*/
	char *transcript;   /* constructed transcript */ 
	int junc_pos;       /* junction position on transcript */ 
	int *S1;            /* stores exon jump position of gene1 */
	int *S2;            /* stores exon jump position of gene2 */
	int S1_num;         /* number of exon jump positions of gene1 */
	int S2_num;         /* number of exon jump positions of gene2 */
	size_t hits;        /* number of times the junction is hitted by reads */
	double likehood;    /* likelihood of the junction */
    UT_hash_handle hh;
} junction_t;

/*
 * the BAG_uthash structure
 */
typedef struct{
	char *edge;
	size_t weight;
	char **read_names;  /* stores the name of read pair that support this edge*/
	char **evidence;    /* stores the read pair that support this edge*/
	junction_t *junc;   /* stores the junctions of the edge, NULL if no junction identified */
    UT_hash_handle hh;  /* makes this structure hashable */
} bag_t;

static inline bag_t 
*bag_init() {
	bag_t *t = mycalloc(1, bag_t);
	t->edge = NULL;
	t->weight = 0;
	t->evidence = mycalloc(1, char*);
	t->read_names = mycalloc(1, char*);
	t->junc = NULL;
	return t;
}

static inline bag_t
*find_edge(bag_t *bag, char* quary) {
	if(quary == NULL) return NULL;
	bag_t* edge = NULL;	
    HASH_FIND_STR(bag, quary, edge);  /* s: output pointer */
	return edge;
}

static inline int 
bag_distory(bag_t **bag) {
	if(*bag == NULL) return 0;
	/*free the kmer_hash table*/
	register bag_t *bag_cur, *bag_tmp;
	junction_t *junc_cur, *junc_tmp;
	HASH_ITER(hh, *bag, bag_cur, bag_tmp) {
		if(bag_cur->junc != NULL){ // if the edge has junctions
			HASH_ITER(hh, bag_cur->junc, junc_cur, junc_tmp){
				HASH_DEL(bag_cur->junc, junc_cur); 			
				free(junc_cur);	
			}			
		}
		HASH_DEL(*bag, bag_cur); 
		free(bag_cur);   
    }
	return 0;
}

/*
 * display bag object
 *
 * PARAMETERS:	bag_t *
 * RETURN:	error code
 */
static inline int bag_display(bag_t *bag) {
	if(bag == NULL) return -1;
	register bag_t *bag_cur, *bag_tmp;
	HASH_ITER(hh, bag, bag_cur, bag_tmp) {
		int i; for(i=0; i < bag_cur->weight; i++){
			printf(">%s\t%s\n%s\n", bag_cur->edge, bag_cur->read_names[i], bag_cur->evidence[i]);
		}
	}
	return 0;
}

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

static inline junction_t 
*find_junction(junction_t *junc, char* quary) {
	if(quary == NULL) return NULL;
	junction_t *s;
    HASH_FIND_STR(junc, quary, s);  /* s: output pointer */
	return s;
}

static inline int 
junction_destory(junction_t **junc){
	if(*junc==NULL) return -1;
	junction_t *junc_cur, *junc_tmp;
	HASH_ITER(hh, *junc, junc_cur, junc_tmp) {
		HASH_DEL(*junc, junc_cur);
		free(junc_cur);
	}
	return 0;
}

/*
 * add one edge to graph
 */
static inline int 
bag_add(bag_t** bag, char* edge_name, char* read_name, char* evidence){
	if(edge_name == NULL || evidence == NULL) return -1;
	bag_t *bag_cur;
	register int n;
	if((bag_cur = find_edge(*bag, edge_name)) == NULL){ /* if edge does not exist */
		bag_cur = bag_init();
		bag_cur->edge = strdup(edge_name);
		bag_cur->weight = 1;
		bag_cur->read_names[0] = strdup(read_name);
		bag_cur->evidence[0] = strdup(evidence);       /* first and only 1 element*/
		HASH_ADD_STR(*bag, edge, bag_cur);								
	}else{
		bag_cur->weight++;
		char **tmp1 = mycalloc(bag_cur->weight, char*);
		char **tmp2 = mycalloc(bag_cur->weight, char*);
		
	 	for (n = 0; n < bag_cur->weight-1; n++){
	 		tmp1[n] = strdup(bag_cur->evidence[n]);
	 	}
	 	tmp1[n] = strdup(evidence);
	 	free(bag_cur->evidence);
	 	bag_cur->evidence = tmp1;
		
	 	for (n = 0; n < bag_cur->weight-1; n++){
	 		tmp2[n] = strdup(bag_cur->read_names[n]);
	 	}
	 	tmp2[n] = strdup(read_name);
	 	free(bag_cur->read_names);
	 	bag_cur->read_names = tmp2;
	}
	return 0;
}

/*
 * delete edges in bag with evidence less than min_weight
 */
static inline bag_t
*bag_trim(bag_t* bag, int min_weight){
	if(bag == NULL) return NULL;	
	
	register bag_t *bag_cur, *bag_tmp, *bag_s, *bag_res=NULL;
	HASH_ITER(hh, bag, bag_cur, bag_tmp) {
		if(bag_cur->weight >= min_weight){
			if((bag_s = find_edge(bag_res, bag_cur->edge))==NULL){ /* this edge not exist */
				HASH_ADD_STR(bag_res, edge, bag_cur);
			}
		}
    }
	return bag_res;
}

/*
 * remove duplicate reads that support graph edge, make sure read pairs 
 * that support every egde is unique.
 */
static inline bag_t 
*bag_uniq(bag_t *bag){
	if(bag == NULL) return NULL;
	bag_t *bag_cur, *bag_tmp, *bag_res=NULL;
	int i, j, repeat;
	/* iterate every edge and remove duplicates */
	for(bag_cur=bag; bag_cur != NULL; bag_cur=bag_cur->hh.next){
		if((bag_tmp = find_edge(bag_res, bag_cur->edge)) == NULL){
			bag_tmp = bag_init();
			bag_tmp->edge = strdup(bag_cur->edge);
			bag_tmp->weight = 0;
			bag_tmp->evidence = mycalloc(bag_cur->weight, char*);
			bag_tmp->read_names = mycalloc(bag_cur->weight, char*);
			bag_tmp->junc = bag_cur->junc;
			/* only get unique reads */
			for(i=0; i<bag_cur->weight; i++){ /* iterate every evidence */
				repeat = false;
				for(j=0; j<bag_tmp->weight; j++)
					if(strcmp(bag_cur->evidence[i], bag_tmp->evidence[j])==0){repeat = true;}
				if(repeat==false){ // no duplicates
					bag_tmp->evidence[bag_tmp->weight] = strdup(bag_cur->evidence[i]);
					bag_tmp->read_names[bag_tmp->weight] = strdup(bag_cur->read_names[i]);
					bag_tmp->weight ++;
				}
			}
			HASH_ADD_STR(bag_res, edge, bag_tmp);
		}		
	}
	return bag_res;
}

static inline int min_mismatch(char* str, char* pattern){
	if(str == NULL || pattern == NULL) return INT_MAX;
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
#endif
