#ifndef BAG_UTHASH_H
#define BAG_UTHASH_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "uthash.h"

struct BAG_uthash {
	char *edge;
	int weight;
	//char **evidence;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

int BAG_uthash_add(struct BAG_uthash** graph_ht, char* edge_name, char* evidence);
void BAG_uthash_destroy(struct BAG_uthash **table);
void BAG_uthash_display(struct BAG_uthash *table);

#endif
