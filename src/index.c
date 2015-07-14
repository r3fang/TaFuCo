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
#include "common.h"

/*--------------------------------------------------------------------*/
/* The name of the file. */
static const char *pcPgmName="index.c";

KSEQ_INIT(gzFile, gzread);
/*--------------------------------------------------------------------*/

static void add_to_kmer_hash(struct kmer_uthash **table, char kmer[MAX_K], char* pos, int k_index);
static void kmer_uthash_write(struct kmer_uthash *htable, char *fname);
int index_main(char *fasta_file, int k);

/* add one kmer and its position to kmer_uthash table */
static void add_to_kmer_hash(struct kmer_uthash **table, char kmer[MAX_K], char* pos, int k_index) {
	struct kmer_uthash *s;	
	assert(table != NULL);
	assert(kmer != NULL);
	/* check if kmer exists in table*/
	HASH_FIND_STR(*table, kmer, s);  

	if (s==NULL){
		s = (struct kmer_uthash*)malloc(sizeof(struct kmer_uthash));
		if (s == NULL)
			{perror(pcPgmName); exit(EXIT_FAILURE);}
			
		strncpy(s->kmer, kmer, k_index);
		//strcpy(s->kmer, kmer);		
		s->kmer[k_index] = '\0';     /*IMPORTANT*/
		s->count = 1;                /* first pos in the list */
		/* an array of char pointers */
		char **tmp;
		s->pos = malloc(s->count * sizeof(char*));
		s->pos = malloc((strlen(pos)+1) * sizeof(char));
		s->pos[0] = pos; /* first and only 1 element*/
		HASH_ADD_STR(*table, kmer, s);
	}else{
		char **tmp;
		s->count += 1;
		/* copy s->pos */
		tmp = malloc(s->count * sizeof(char*));
		int i;
		for (i = 0; i < s->count-1; i++){
		    tmp[i] = malloc((strlen(s->pos[i])+1) * sizeof(char));
			tmp[i] = s->pos[i];
		}
		/* append pos */
		tmp[i] = malloc(strlen(pos) * sizeof(char));
		tmp[i] = pos;
		/* assign tmp to s->pos*/
		s->pos = tmp;
	}
}
/*--------------------------------------------------------------------*/

/* Write down kmer_uthash */
static void kmer_uthash_write(struct kmer_uthash *htable, char *fname){
	/* write htable to disk*/
	FILE *ofp = fopen(fname, "w");

	if (ofp == NULL) {
	  fprintf(stderr, "Can't open output file %s!\n", fname);
	  exit(1);
	}

	if(htable == NULL)
		exit(1);
	
	struct kmer_uthash *s, *tmp;

	HASH_ITER(hh, htable, s, tmp) {

		if(s == NULL)
			fprintf(stderr, "Fail to write down %s!\n", fname);
		
		fprintf(ofp, ">%s\t%d\n", s->kmer, s->count);		
		int i;
		for(i=0; i < s->count; i++){
			if(i==0){
				fprintf(ofp, "%s", s->pos[i]);																
			}else{
				fprintf(ofp, "|%s", s->pos[i]);
			}
		}
		fprintf(ofp, "\n");
	}
	fclose(ofp);
}

/*--------------------------------------------------------------------*/
int index_main(char *fasta_file, int k){	
	if(k > MAX_K)
		{fprintf(stderr, "ERROR: input k exceeds 100\n"); exit(EXIT_FAILURE);}	
	/* index file */
	gzFile fp;  
	kseq_t *seqs;  
	int l;
	
	struct kmer_uthash *table = NULL;

	fp = gzopen(fasta_file, "r");
	if (fp == NULL)
		{perror(fasta_file); exit(EXIT_FAILURE);}
	
	seqs = kseq_init(fp);	
	while ((l = kseq_read(seqs)) >= 0) {
		char *seq = strToUpper(seqs->seq.s);
		char *name = seqs->name.s;		
		if(seq == NULL || name == NULL || strlen(seq) <= k){
			continue;
		}
		int i;
		for(i=0; i < strlen(seq)-k+1; i++){
			char kmer[MAX_K];
			if(kmer == NULL)
				continue;
			memset(kmer, '\0', sizeof(kmer));
			memcpy(kmer, seq+i, k);
			char i_str[100];
			sprintf(i_str, "%d", i);
			add_to_kmer_hash(&table, small_dna_str(kmer, rev_com(kmer)), concat(concat(name, "_"), i_str), k); 
		}
	}  	
	char *index_file = concat(fasta_file, ".index");	
	if(index_file == NULL)
		{exit(EXIT_FAILURE);}

	kmer_uthash_write(table, index_file);
	kmer_uthash_destroy(table);
	kseq_destroy(seqs);
	gzclose(fp);
	return 0;
}