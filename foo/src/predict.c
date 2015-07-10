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

struct kmer_uthash *load_kmer_htable(char *fname){
	/* load kmer_htable*/
	struct kmer_uthash *htable = NULL;
	gzFile fp;  
	kseq_t *seq;  
	int l;
	fp = gzopen(fname, "r");
	if (fp == NULL) {
	  fprintf(stderr, "Can't open input file %s!\n", fname);
	  exit(1);
	}
	seq = kseq_init(fp); // STEP 3: initialize seq  
	if (seq == NULL){
		return NULL;
	}
	while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence 
		char *kmer = seq->name.s;
		char *pos = seq->seq.s;
		int k = strlen(kmer);
		struct kmer_uthash *s;
		s = (struct kmer_uthash*)malloc(sizeof(struct kmer_uthash));
		strncpy(s->kmer, kmer, k);
		s->kmer[k] = '\0'; /* just in case*/
		s->pos = pos;
	}
	gzclose(fp); // STEP 6: close the file handler  
	kseq_destroy(seq);
	return(htable);
}

int predict_main(char *fasta_file, char *fastq_file){
	/* main function*/
	char *index_file = concat(fasta_file, ".index");
	struct kmer_uthash *htable = load_kmer_htable(index_file);
	
	gzFile fp;  
	kseq_t *seq;  
	int l;
	fp = gzopen(fastq_file, "r");
	if (fp == NULL) {
	  fprintf(stderr, "Can't open input file %s!\n", fastq_file);
	  exit(1);
	}
	seq = kseq_init(fp); // STEP 3: initialize seq  
	if (seq == NULL){
		return NULL;
	}
	while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence 
		char *name = seq->name.s;
		char *end = seq->seq.s;
		printf("%s\t%s\n", name, end);
	}
	gzclose(fp); // STEP 6: close the file handler  
	kseq_destroy(seq);
	kmer_table_destroy(&htable);
	return 0;
}

