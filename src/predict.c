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
		/* add kmer */
		char *kmer = seq->name.s;
		int k = strlen(kmer);
		struct kmer_uthash *s;
		s = (struct kmer_uthash*)malloc(sizeof(struct kmer_uthash));
		strncpy(s->kmer, kmer, k);
		s->kmer[k] = '\0'; /* just in case*/
		
		/* add count */
		int count = atoi(seq->comment.s);
		s->count = count;
		
		/* add pos */
		char *pos = seq->seq.s;				
		s->pos = malloc((s->count) * sizeof(char*));
		
		/* split a string by delim */
		char *token;	    
		/* get the first token */
	    token = strtok(pos, "|");				
		/* walk through other tokens */
		int i = 0;
	    while(token != NULL) 
	    {
			s->pos[i] = malloc((strlen(token)+1) * sizeof(char));
			/*duplicate a string*/
			s->pos[i] = strdup(token);
			token = strtok(NULL, "|");
			i ++;
	    }		
		HASH_ADD_STR(htable, kmer, s);
	}	
	gzclose(fp); // STEP 6: close the file handler  
	kseq_destroy(seq);
	return(htable);
}

struct fasta_uthash *fasta_uthash_load(char *fname){
	gzFile fp;
	kseq_t *seq;
	int l;
	struct fasta_uthash *res = NULL;
	fp = gzopen(fname, "r");
	if(fp == NULL){
		return NULL;
	}
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		struct fasta_uthash *s;
		s = (struct fasta_uthash*)malloc(sizeof(struct fasta_uthash));
		/* if we kseq_destroy(seq);, we need to duplicate the string!!!*/
		s->name = strdup(seq->name.s);
		s->seq = strdup(strToUpper(seq->seq.s));
		if(seq->comment.l) s->comment = seq->comment.s;
		HASH_ADD_STR(res, name, s);
	}	
	kseq_destroy(seq);
	gzclose(fp);
	return res;
}

char* 
find_mpm(char *_read, int *pos_read, int k, struct kmer_uthash **kmer_ht, struct fasta_uthash **fasta_ht){
	/* copy a part of string */
	char buff[k];
	strncpy(buff, _read + *pos_read, k);
	if(buff == NULL || strlen(buff) != k){
		*pos_read += 1; // move the pointer one base forward
		return NULL;
	}
	struct kmer_uthash *s_kmer;
	HASH_FIND_STR(*kmer_ht, buff, s_kmer);
	if(s_kmer==NULL){ // kmer does not exist in the hash table
		*pos_read += 1;
		return NULL;
	}
	int pos_seq;
	char *tmp = strdup(s_kmer->pos[0]);
	assert(tmp);
	char* exon = pos_parser(tmp, &pos_seq);
	/* error report*/
	if(exon==NULL || pos_seq==NULL || pos_seq<0){
		*pos_read += 1;
		return NULL;
	}
	struct fasta_uthash *s_fasta;
	HASH_FIND_STR(*fasta_ht, exon, s_fasta);
	if(s_fasta == NULL){
		* pos_read += 1;
		return NULL;
	}
	char *_seq = strdup(s_fasta->seq);
	int m = 0;
	while(*(_seq + m + pos_seq) == *(_read+ *pos_read + m)){
		m ++;
	}				
	*pos_read += m;
	return exon;
}

int find_exons_on_read(char* _exons[], char* _read, int k, struct kmer_uthash **kmer_ht, struct fasta_uthash **fasta_ht){
	int pos_read = 0;
	int i = 0; 
	while(pos_read<(strlen(_read)-k)){
		char* exon = find_mpm(_read, &pos_read, k, kmer_ht, fasta_ht);
		if (exon!=NULL){
			printf("exon=%s\n", exon);
			_exons[i++] = strdup(exon);
			//printf("i=%d\texon=%s\n", i, _exons[0]);
		}
		free(exon);
	}
	return i;
}

int predict_main(char *fasta_file, char *fastq_file, int k){
	/* load kmer hash table in the memory */
	char *index_file = concat(fasta_file, ".index");	
	/* load kmer_uthash table */
	struct kmer_uthash *kmer_ht = load_kmer_htable(index_file);
	/* load fasta_uthash table */
	struct fasta_uthash *fasta_ht = fasta_uthash_load(fasta_file);

	/* starting parsing and processing fastq file*/
	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen(fastq_file, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		char *_read = strdup(seq->seq.s);
		char **exons = calloc(100, sizeof(char *));
		assert(exons);
		int num = find_exons_on_read(exons, _read, k, &kmer_ht, &fasta_ht);
		free(exons);
	}

	
	//struct kmer_uthash *s, *tmp;
	//HASH_ITER(hh, kmer_ht, s, tmp) {
	//	if(s->count > 1){
	//		printf("%s\n", s->kmer);
	//		for(int i=0; i<s->count; i++)
	//		printf("%s\n", s->pos[i]);						
	//	}
	//}
	
	kseq_destroy(seq);
	gzclose(fp);
	kmer_uthash_destroy(&kmer_ht);	
	fasta_uthash_destroy(&fasta_ht);	
	return 0;
}

