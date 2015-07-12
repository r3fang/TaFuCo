/*--------------------------------------------------------------------*/
/* index.c                                                            */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/*--------------------------------------------------------------------*/

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include <errno.h>
#include "uthash.h"
#include "utlist.h"
#include "utstring.h"
#include "utarray.h"
#include "kseq.h"
#include "common.h"

/*--------------------------------------------------------------------*/
/* The name of the file. */
static const char *pcPgmName="index.c";

KSEQ_INIT(gzFile, gzread)  
/*--------------------------------------------------------------------*/

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
		int i;
		s->count += 1;
		/* copy s->pos */
		tmp = malloc(s->count * sizeof(char*));
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
static void write_kmer_htable(struct kmer_uthash **htable, char *fname){
	/* write htable to disk*/
	FILE *ofp = fopen(fname, "w");

	if (ofp == NULL) {
	  fprintf(stderr, "Can't open output file %s!\n", fname);
	  exit(1);
	}

	if(htable == NULL)
		exit(1);
	
	struct kmer_uthash *s, *tmp;

	HASH_ITER(hh, *htable, s, tmp) {

		if(s == NULL)
			fprintf(stderr, "Fail to write down %s!\n", fname);
		
		fprintf(ofp, ">%s\t%d\n", s->kmer, s->count);		
		for(int i=0; i < s->count; i++){
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
		if(seq == NULL)
			{perror(pcPgmName); exit(EXIT_FAILURE);}
				
		char *name = seqs->name.s;
		if (name==NULL)
			{perror(pcPgmName); exit(EXIT_FAILURE);}
		
		for(int i=0; i < strlen(seq)-k+1; i++){
			char kmer[MAX_K];
			if(kmer == NULL)
				{perror(pcPgmName); exit(EXIT_FAILURE);}
			
			memset(kmer, '\0', sizeof(kmer));
			memcpy(kmer, seq+i, k);
			char i_str[100];
			sprintf(i_str, "%d", i);
			add_to_kmer_hash(&table, kmer, concat(concat(name, "_"), i_str), k); 
		}
	}  
	char *index_file = concat(fasta_file, ".index");
	if(index_file)
		{exit(EXIT_FAILURE);}
	write_kmer_htable(&table, index_file);
	kmer_uthash_destroy(&table);
	kseq_destroy(seqs);
	gzclose(fp);
	return 0;
}