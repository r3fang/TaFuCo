/* Contact: Rongxin Fang <r3fang@ucsd.edu> */
/* Last Modified: 15JULY2015 */

#ifndef _BAG_UTHASH_H
#define _BAG_UTHASH_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>

/* error code */
#define BA_ERR_NONE		     0 // no error
#define BA_ERR_PARAM		-1 // bad paraemter
#define BA_ERR_HASHITER		-2 // fail to iterate uthash
#define BA_ERR_MALLOC		-3 // fail to malloc memory uthash

/*
 * error handling
 * all errors are returned as an integer code, and a string
 * amplifying the error is saved in here; this can then be
 * printed
 */
char qe_errbuf[256] = "no error\n";	/* the error message buffer */
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

/*
 * free BAG_uthash table
 *
 * PARAMETERS:	struct BAG_uthash **
 * RETURN:	error code
 */

static inline int BAG_uthash_destroy(struct BAG_uthash **table) {
	int error;
	if(*table == NULL) goto FAIL_PARAM;
	/*free the kmer_hash table*/
	register struct BAG_uthash *cur, *tmp;
	HASH_ITER(hh, *table, cur, tmp) {
		if(cur==NULL) goto FAIL_HASH_ITER;
		HASH_DEL(*table, cur);  /* delete it (users advances to next) */
		free(cur);            /* free it */
    }
	goto SUCCESS;


	FAIL_PARAM:
		ERRBUF("BAG_uthash_add: invalid NULL parameter");
		error = BA_ERR_PARAM;
		goto EXIT;		
	FAIL_MALLOC:
		ERRBUF("BAG_uthash_display: fail to malloc memory");
		error = BA_ERR_MALLOC;
		goto EXIT;
	FAIL_HASH_ITER:
		ERRBUF("BAG_uthash_display: fail to iterate BAG_uthash");
		error = BA_ERR_HASHITER;
		goto EXIT;
	SUCCESS:
		error = BA_ERR_NONE;
		goto EXIT;
	EXIT:
		return error;

}

/*
 * display BAG_uthash table on the screen
 *
 * PARAMETERS:	struct BAG_uthash *
 * RETURN:	error code
 */
static inline int BAG_uthash_display(struct BAG_uthash *graph_ht) {
	int error;
	if(graph_ht == NULL) goto FAIL_PARAM;

	/*free the kmer_hash table*/
	register struct BAG_uthash *cur, *tmp;
	HASH_ITER(hh, graph_ht, cur, tmp) {
		if(cur == NULL) goto FAIL_HASH_ITER;
		printf(">%s\t%zu\n", cur->edge, cur->weight);
		int i;
		for(i=0; i<cur->weight; i++){
			printf("%s\n", cur->evidence[i]);	
		}
	}
	goto SUCCESS;

	FAIL_PARAM:
		ERRBUF("BAG_uthash_add: invalid NULL parameter");
		error = BA_ERR_PARAM;
		goto EXIT;		
	FAIL_MALLOC:
		ERRBUF("BAG_uthash_display: fail to malloc memory");
		error = BA_ERR_MALLOC;
		goto EXIT;
	FAIL_HASH_ITER:
		ERRBUF("BAG_uthash_display: fail to iterate BAG_uthash");
		error = BA_ERR_HASHITER;
		goto EXIT;
	SUCCESS:
		error = BA_ERR_NONE;
		goto EXIT;
	EXIT:
		return error;

}

/*
 * add one edge to graph
 */
static inline int BAG_uthash_add(struct BAG_uthash** graph_ht, char* edge_name, char* evidence){
	int error;
	/* check parameters */
	if(*graph_ht == NULL || edge_name == NULL || evidence == NULL) goto FAIL_PARAM;

	register struct BAG_uthash *s;
	HASH_FIND_STR(*graph_ht, edge_name, s);

	if(s==NULL){
		if((s=(struct BAG_uthash*)malloc(sizeof(struct BAG_uthash))) == NULL) goto FAIL_MALLOC;
		s->edge = edge_name;
		s->weight = 1;
		if((s->evidence = malloc(sizeof(char*))) == NULL) goto FAIL_MALLOC;
		s->evidence[0] = evidence; /* first and only 1 element*/
		HASH_ADD_STR(*graph_ht, edge, s);						
	}else{
		s->weight ++;
		char **tmp;
		if((tmp = malloc((s->weight+1) * sizeof(char*)))==NULL) goto FAIL_MALLOC;
		int n;
		for (n = 0; n < s->weight-1; n++){
			tmp[n] = strdup(s->evidence[n]);
		}
		free(s->evidence);
		tmp[n] = evidence;
		s->evidence = tmp;
	}
	goto SUCCESS;
	
	FAIL_PARAM:
		ERRBUF("BAG_uthash_add: invalid NULL parameter");
		error = BA_ERR_PARAM;
		goto EXIT;		
	FAIL_MALLOC:
		ERRBUF("BAG_uthash_display: fail to malloc memory");
		error = BA_ERR_MALLOC;
		goto EXIT;
	FAIL_HASH_ITER:
		ERRBUF("BAG_uthash_display: fail to iterate BAG_uthash");
		error = BA_ERR_HASHITER;
		goto EXIT;
	SUCCESS:
		error = BA_ERR_NONE;
		goto EXIT;
	EXIT:
		return error;
	
}

#endif
