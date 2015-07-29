/*--------------------------------------------------------------------*/
/* index.c                                                            */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Indexing reference exon sequences.                                 */
/*--------------------------------------------------------------------*/
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h>
#include <zlib.h> 
#include <errno.h>
#include <assert.h>
#include "uthash.h"
#include "kseq.h"
#include "kstring.h"
#include "common.h"
#include "utils.h"
#include "kmer_uthash.h"
#include "fasta_uthash.h"


/*--------------------------------------------------------------------*/
int main(int argc, char *argv[]) { 
	if (argc != 3) {  
	        fprintf(stderr, "Usage: %s <in.fa> <k>\n", argv[0]);  
	        return 1;  
	 }  	 
	char *fasta_file = argv[1];
	
	int k;
	if (sscanf(argv[2], "%i", &k)!=1) {printf ("error - k not an integer");}	
	struct kmer_uthash *tb = kmer_uthash_construct(fasta_file, k);
	/* index file */
	char *index_file = concat(fasta_file, ".index");	
	if(index_file == NULL) die("output file error\n");
	kmer_uthash_write(tb, index_file);
	kmer_uthash_destroy(&tb);
	return 0;
}