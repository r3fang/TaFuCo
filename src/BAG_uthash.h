/* Contact: Rongxin Fang <r3fang@ucsd.edu> */
/* Last Modified: 15JULY2015 */

#ifndef _BAG_UTHASH_H
#define _BAG_UTHASH_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "common.h"
#include "utils.h"

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
BAG_uthash_tim(struct BAG_uthash** tb, int min_weight){
	if(*tb == NULL || min_weight < 0) die("BAG_uthash_tim: wrong parameters\n");
	
	register struct BAG_uthash *cur, *tmp;
	HASH_ITER(hh, *tb, cur, tmp) {
		if(cur==NULL) die("BAG_uthash_tim: HASH_ITER fails\n");
		if(cur->weight < min_weight) {HASH_DEL(*tb, cur); free(cur);}
    }
	return BA_ERR_NONE;
}
#endif
