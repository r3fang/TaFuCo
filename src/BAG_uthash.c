#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "uthash.h"
#include "BAG_uthash.h"

void 
BAG_uthash_destroy(struct BAG_uthash **table) {
	/*free the kmer_hash table*/
  struct BAG_uthash *cur, *tmp;
  HASH_ITER(hh, *table, cur, tmp) {
      HASH_DEL(*table, cur);  /* delete it (users advances to next) */
      free(cur);            /* free it */
    }
}

void 
BAG_uthash_display(struct BAG_uthash *table) {
	/*free the kmer_hash table*/
  struct BAG_uthash *cur, *tmp;
  HASH_ITER(hh, table, cur, tmp) {
	  if(cur == NULL)
		  exit(-1);
	  printf(">%s\t%d\n", cur->edge, cur->weight);
	  int i;
	  for(i=0; i<cur->weight; i++){
	  	  printf("%s\n", cur->evidence[i]);	
	   }
    }
}

/* Add one edge to BAG.             */
int 
BAG_uthash_add(struct BAG_uthash** graph_ht, char* edge_name, char* evidence){
	assert(graph_ht != NULL);
	assert(edge_name != NULL);
	assert(evidence != NULL);
	if(edge_name == NULL){
		return 0;
	}
	struct BAG_uthash *s;
	HASH_FIND_STR(*graph_ht, edge_name, s);
	if(s==NULL){
		s = (struct BAG_uthash*)malloc(sizeof(struct BAG_uthash));
		s->edge = strdup(edge_name);
		s->weight = 1;
		s->evidence = malloc(sizeof(char*));
		s->evidence[0] = evidence; /* first and only 1 element*/
		HASH_ADD_STR(*graph_ht, edge, s);						
	}else{
		s->weight += 1;
		char **tmp = malloc((s->weight+1) * sizeof(char*));
		int n;
		for (n = 0; n < s->weight-1; n++){
			tmp[n] = strdup(s->evidence[n]);
		}
		tmp[n] = evidence;
		free(s->evidence);
		/* assign tmp to s->pos*/
		s->evidence = tmp;
	}
	return 0;
}

