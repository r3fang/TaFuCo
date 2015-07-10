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
KSEQ_INIT(gzFile, gzread)  

void add_to_kmer_hash(struct kmer_uthash **table, char kmer[MAX_K], char* pos, int k_index) {
	/* add one key to kmer_hash */
	struct kmer_uthash *s;
	HASH_FIND_STR(*table, kmer, s);  /* check if kmer exists*/
	if (s==NULL){
		s = (struct kmer_uthash*)malloc(sizeof(struct kmer_uthash));
		strncpy(s->kmer, kmer, k_index);
		s->kmer[k_index] = '\0';     /*IMPORTANT*/
		s->pos = pos;
		HASH_ADD_STR(*table, kmer, s);
	}else{
		s->pos = concat(concat(s->pos, "|"), pos); /*append one more position*/
	}
}

int index_main(char *fasta_file, int k){	
	/* index file */
	char *index_file = concat(fasta_file, ".index");
	gzFile fp;  
	kseq_t *seqs;  
	int l;
	fp = gzopen(fasta_file, "r");
	if (fp == NULL){
		return NULL;
	}
	struct kmer_uthash *table = NULL;
	seqs = kseq_init(fp); // STEP 3: initialize seq  
	while ((l = kseq_read(seqs)) >= 0) { // STEP 4: read sequence 
		char *seq = strToUpper(seqs->seq.s);
		int i;
		if (seqs->name.s==NULL)
			return NULL;
		char *name = seqs->name.s;
		for(i=0; i < strlen(seq)-k+1; i++){
			char kmer[MAX_K];
			memset(kmer, '\0', sizeof(kmer));
			memcpy(kmer, seq+i, k);
			char i_str[100];
			sprintf(i_str, "%d", i);
			add_to_kmer_hash(&table, kmer, concat(concat(name, "."), i_str), k); 
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
	kmer_table_destroy(&table);
	kseq_destroy(seqs);
	gzclose(fp);
	return 0;
}