/* Contact: Rongxin Fang <r3fang@ucsd.edu> */
/* Last Modified: 15JULY2015 */

#ifndef BAG_UTHASH_H
#define BAG_UTHASH_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "uthash.h"

#define BAG_NONE				 0 // no error
#define BAG_BADPARAM			-1 // bad paraemter
#define BAG_BADITER				-2 // fail to iterate uthash
#define BAG_BADLLOC				-3 // fail to iterate uthash

/*
 * error handling
 * all errors are returned as an integer code, and a string
 * amplifying the error is saved in here; this can then be
 * printed
 */
char qe_errbuf[256] = "no error";	/* the error message buffer */
#define ERRBUF(str)			(void) strncpy(qe_errbuf, str, sizeof(qe_errbuf))
#define ERRBUF2(str,n)		(void) sprintf(qe_errbuf, str, n)
#define ERRBUF3(str,n,m)	(void) sprintf(qe_errbuf, str, n, m)

/*
 * the BAG_uthash structure
 */
struct BAG_uthash {
	char *edge;
	size_t weight;
	char **evidence;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

static inline int BAG_uthash_destroy(struct BAG_uthash **table) {
	if(*table == NULL){
		ERRBUF("BAG_uthash_destroy: invalid NULL parameter");
		return(BAG_BADPARAM);
	}
	/*free the kmer_hash table*/
	register struct BAG_uthash *cur, *tmp;
	HASH_ITER(hh, *table, cur, tmp) {
		HASH_DEL(*table, cur);  /* delete it (users advances to next) */
		free(cur);            /* free it */
    }
	return BAG_NONE;
}

static inline int BAG_uthash_display(struct BAG_uthash *graph_ht) {
	if(graph_ht == NULL){
		ERRBUF("BAG_uthash_display: invalid NULL parameter");
		return(BAG_BADPARAM);
	}
	/*free the kmer_hash table*/
	register struct BAG_uthash *cur, *tmp;
	HASH_ITER(hh, graph_ht, cur, tmp) {
		if(cur == NULL){
			ERRBUF("BAG_uthash_display: fail to iterate BAG_uthash");
			return(BAG_BADITER);
		}
		printf(">%s\t%zu\n", cur->edge, cur->weight);
		int i;
		for(i=0; i<cur->weight; i++){
			printf("%s\n", cur->evidence[i]);	
		}
	}
	return BAG_NONE;
}

/* Add one edge to BAG.             */
static inline int BAG_uthash_add(struct BAG_uthash** graph_ht, char* edge_name, char* evidence){
	/* check parameters */
	if(*graph_ht == NULL || edge_name == NULL || evidence == NULL){
		ERRBUF("BAG_uthash_add: invalid NULL parameter");
		return(BAG_BADPARAM);
	}	
	register struct BAG_uthash *s;
	HASH_FIND_STR(*graph_ht, edge_name, s);
	if(s==NULL){
		s = (struct BAG_uthash*)malloc(sizeof(struct BAG_uthash));
		if(s == NULL){
			ERRBUF("BAG_uthash_add: fail to malloc");
			return(BAG_BADLLOC);
		}
		s->edge = edge_name;
		s->weight = 1;
		s->evidence = malloc(sizeof(char*));
		if(s->evidence == NULL){
			ERRBUF("BAG_uthash_add: fail to malloc");
			return(BAG_BADLLOC);
		}
		s->evidence[0] = evidence; /* first and only 1 element*/
		HASH_ADD_STR(*graph_ht, edge, s);						
	}else{
		s->weight ++;
		char **tmp = malloc((s->weight+1) * sizeof(char*));
		int n;
		for (n = 0; n < s->weight-1; n++){
			tmp[n] = strdup(s->evidence[n]);
		}
		tmp[n] = evidence;
		free(s->evidence);
		s->evidence = tmp;
	}
	return BAG_NONE;
}
#endif
