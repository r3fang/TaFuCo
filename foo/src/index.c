#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include "uthash.h"
#include "utlist.h"
#include "utstring.h"
#include "utarray.h"
#include "kseq.h"
#include "common.h"

/* kmer length */
#define k 15

KSEQ_INIT(gzFile, gzread)  
	
struct kmer_uthash {
    char kmer[k];                /* key */
	char *pos;    
	int count;
    UT_hash_handle hh;         /* makes this structure hashable */
};

/* Global variables */
struct kmer_uthash *table = NULL;

void add_to_kmer_hash(char kmer[k], char* pos) {
	/* You really need to pass a pointer to the hash pointer: **table*/
	struct kmer_uthash *s;
	HASH_FIND_STR(table, kmer, s);
	if (s==NULL){
		s = (struct kmer_uthash*)malloc(sizeof(struct kmer_uthash));
		strncpy(s->kmer, kmer, k);
		s->pos = pos;
		s->count = 1;		
		HASH_ADD_STR(table, kmer, s);
	}else{
		s->count += 1;
		s->pos = concat(concat(s->pos, "|"), pos);
	}
}

void kmer_table_destroy() {
  struct kmer_uthash *cur, *tmp;
  
  HASH_ITER(hh, table, cur, tmp) {
      HASH_DEL(table, cur);  /* delete it (users advances to next) */
      free(cur);            /* free it */
    }
}

int index_main(char *fasta_file) {	
	/* index file */
	char *index_file = concat(fasta_file, ".index");
	gzFile fp;  
	kseq_t *seqs;  
	int l;
	fp = gzopen(fasta_file, "r");
	if (fp == NULL){
		return NULL;
	}
	seqs = kseq_init(fp); // STEP 3: initialize seq  
	while ((l = kseq_read(seqs)) >= 0) { // STEP 4: read sequence 
		char *seq = strToUpper(seqs->seq.s);
		int i;
		if (seqs->name.s==NULL)
			return NULL;
		char *name = seqs->name.s;
		char kmer[k];
		for(i=0; i < strlen(seq)-k+1; i++){
			memcpy(kmer, &seq[i], k);
			kmer[k] = '\0';
			/* convert i to string */
			char i_str[100];
			sprintf(i_str, "%d", i);
			add_to_kmer_hash(kmer, concat(concat(name, "."), i_str)); 
		}
	}  
	
	FILE *ofp = fopen(index_file, "w");
	if (ofp == NULL) {
	  fprintf(stderr, "Can't open output file %s!\n", index_file);
	  exit(1);
	}
	struct kmer_uthash *s, *tmp;
	HASH_ITER(hh, table, s, tmp) {
		fprintf(ofp, ">%s\n%s\n", s->kmer, s->pos);						
	}
	/*free everything*/
	fclose(ofp);
	kmer_table_destroy();
	kseq_destroy(seqs);
	gzclose(fp);
	return 0;
}
