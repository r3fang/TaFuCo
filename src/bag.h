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
	char *gname1; // gene1 and gene2 has order
	char *gname2;
	size_t weight;
	char **read_names;  /* stores the name of read pair that support this edge*/
	char **evidence;    /* stores the read pair that support this edge*/
	bool junc_flag;
	junction_t *junc;   /* stores the junctions of the edge, NULL if no junction identified */	
    UT_hash_handle hh;  /* makes this structure hashable */
} bag_t;

static inline bag_t 
*bag_init() {
	bag_t *t = mycalloc(1, bag_t);
	t->gname1 = NULL;
	t->gname1 = NULL;	
	t->edge = NULL;
	t->junc_flag = false;
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
junction_destory(junction_t **junc){
	if(*junc==NULL) return -1;
	junction_t *junc_cur, *junc_tmp;
	HASH_ITER(hh, *junc, junc_cur, junc_tmp) {
		HASH_DEL(*junc, junc_cur);
		if(junc_cur->idx)             free(junc_cur->idx);
		if(junc_cur->exon1)           free(junc_cur->exon1);
		if(junc_cur->exon2)           free(junc_cur->exon2);
		if(junc_cur->s)               free(junc_cur->s);
		if(junc_cur->transcript)      free(junc_cur->transcript);
		if(junc_cur->S1)              free(junc_cur->S1);
		if(junc_cur->S2)              free(junc_cur->S2);
		free(junc_cur);
	}
	return 0;
}

static inline int 
bag_destory(bag_t **bag) {
	if(*bag == NULL) return 0;
	/*free the kmer_hash table*/
	register bag_t *bag_cur, *bag_tmp;
	junction_t *junc_cur, *junc_tmp;
	int i;
	HASH_ITER(hh, *bag, bag_cur, bag_tmp) {
		HASH_DEL(*bag, bag_cur); 
		if(bag_cur->junc_flag==true){
			junction_destory(&bag_cur->junc);			
		}else{
			if(bag_cur->junc->idx)        free(bag_cur->junc->idx);
			if(bag_cur->junc->exon1)      free(bag_cur->junc->exon1);
			if(bag_cur->junc->exon2)      free(bag_cur->junc->exon2);
			if(bag_cur->junc->s)          free(bag_cur->junc->s);
			if(bag_cur->junc->transcript) free(bag_cur->junc->transcript);
			if(bag_cur->junc->S1)         free(bag_cur->junc->S1);
			if(bag_cur->junc->S2)         free(bag_cur->junc->S2);			
		}
		if(bag_cur->edge)                 free(bag_cur->edge);
		if(bag_cur->gname1)               free(bag_cur->gname1);
		if(bag_cur->gname2)               free(bag_cur->gname2);
		for(i=0; i<bag_cur->weight; i++)  free(bag_cur->evidence[i]);
	    for(i=0; i<bag_cur->weight; i++)  free(bag_cur->read_names[i]);
		if(bag_cur->evidence)             free(bag_cur->evidence);
		if(bag_cur->read_names)           free(bag_cur->read_names);
		if(bag_cur)                       free(bag_cur);
	}
	return 0;
}

static inline int junction_display(junction_t *junc){
	if(junc == NULL) return -1;
	junction_t *junc_cur;
	for(junc_cur=junc; junc_cur!=NULL; junc_cur=junc_cur->hh.next){
		printf("exon1=%s\texon2=%s\thits=%zu\tlikehood=%f\n", junc_cur->exon1, junc_cur->exon2, junc_cur->hits, junc_cur->likehood);
		if(junc_cur->s != NULL)      printf("junc_str=%s\n", junc_cur->s);
		if(junc_cur->transcript != NULL)  printf("transcript=%s\n", junc_cur->transcript);
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
	register bag_t *bag_cur;
	for(bag_cur=bag; bag_cur!=NULL; bag_cur=bag_cur->hh.next){
		printf("Fusion:\n---------\n");
		printf("%s\t%s\tweight=%zu\n\n", bag_cur->gname1, bag_cur->gname2, bag_cur->weight);
		printf("Junction:\n-------\n");
		//if(bag_cur->junc_flag == true) junction_display(bag_cur->junc);
		junction_display(bag_cur->junc);
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
	junc->s = mycalloc(seed_len+3, char);	
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
static inline int
bag_trim(bag_t **bag, int min_weight){
	if(*bag==NULL) return -1;
	register bag_t *bag_cur, *bag_tmp;
	HASH_ITER(hh, *bag, bag_cur, bag_tmp) {
		if(bag_cur->weight < min_weight){
			HASH_DEL(*bag, bag_cur);
			free(bag_cur);
		}
	}
	return 0;
}


/*
 * remove duplicate reads that support graph edge, make sure read pairs 
 * that support every egde is unique.
 */
static inline int 
bag_uniq(bag_t **bag){
	if(*bag==NULL) return -1;
	bag_t *bag_cur;
	int i, j;
	bool repeat;
	int weight;
	char **names, **evidence;
	/* iterate every edge and remove duplicates */
	for(bag_cur=*bag; bag_cur != NULL; bag_cur=bag_cur->hh.next){
		evidence = mycalloc(bag_cur->weight, char*);
		names = mycalloc(bag_cur->weight, char*);
		weight = 0;
		for(i=0; i<bag_cur->weight; i++){ /* iterate every evidence */
			repeat = false;
			for(j=0; j<weight; j++){if(strcmp(bag_cur->evidence[i], evidence[j])==0){repeat = true;}}
			if(repeat==false){ // no duplicates
				evidence[weight] = strdup(bag_cur->evidence[i]);
				names[weight] = strdup(bag_cur->read_names[i]);
				weight ++;
			}
		}
		if(bag_cur->evidence)    free(bag_cur->evidence);
		if(bag_cur->read_names)  free(bag_cur->read_names);
		bag_cur->evidence = evidence;
		bag_cur->read_names = names;
		bag_cur->weight = weight;
	}
	return 0;
}
/*
 * min mismatch between long string and a short pattern
 */
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
