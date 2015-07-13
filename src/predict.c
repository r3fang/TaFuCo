#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include <assert.h>
#include "uthash.h"
#include "kseq.h"
#include "common.h"
KSEQ_INIT(gzFile, gzread)  


static const char *pcPgmName="predict.c";
static struct kmer_uthash *KMER_HT = NULL;
static struct fasta_uthash *FASTA_HT = NULL;


char* 
find_mpm(char *_read, int pos_read, int k){
	/* copy a part of string */
	char _buff[k];
	strncpy(_buff, _read + pos_read, k);
	_buff[k] = '\0';
	char *buff = small_dna_str(_buff, rev_com(_buff));
	if(buff == NULL || strlen(buff) != k){
		return NULL;
	}
		
	struct kmer_uthash *s_kmer = find_kmer(buff, KMER_HT);
	if(s_kmer==NULL){
		return NULL;
	}
	
	char *max_exon;
	int *max_len_list = malloc(s_kmer->count * sizeof(int));
	int max_len = 0;
		
	/* discard if it matches more than 1 prefix kmer */
	for(int i =0; i < s_kmer->count; i++){		
		int _pos_exon; // position on exon
		
		char *exon = pos_parser(s_kmer->pos[i], &_pos_exon);
		/* error report*/
		if(exon==NULL || _pos_exon==NULL || _pos_exon<0){
			return NULL;
		}
		
		struct fasta_uthash *s_fasta = find_fasta(exon, FASTA_HT);
		if(s_fasta == NULL){
			return NULL;
		}
		char *_seq = strdup(s_fasta->seq);
		
		int m = 0;
		/* extending kmer to find MPM */
		while(*(_seq + m + _pos_exon) == *(_read+ pos_read + m)){
			m ++;
		}
		
		max_len_list[i] = m;
		if(m > max_len){
			max_len = m;
			max_exon = strdup(exon);
		}
		
		free(_seq);
		free(exon);
	}
	
	int max_count = 0; // count how many MPM found
	for(int i=0; i < s_kmer->count; i++)
		if(max_len == max_len_list[i])
			max_count ++;
	
	if(max_count == 1 && strlen(max_exon) != NULL)
		return max_exon;							

	return NULL;
}

int 
find_MPMs_on_read(struct MPM *_qptr, char* _read, char* _read_name, int _k){
	int _read_pos = 0;
	int _i = 0; 
	char** hits = malloc(strlen(_read) * sizeof(char*));  
	while(_read_pos<(strlen(_read)-_k)){
		char* _exon = find_mpm(_read, _read_pos, _k);
		_read_pos += 1;
		if (_exon != NULL){
			hits[_i] = strdup(_exon);
			free(_exon);
			_i ++;
		}
	}
	for(int j=0; j<_i; j++)
		printf("%s\t", hits[j]);
	printf("\n");
	return _i;
}

char* 
construct_BAG(char *_fastq_file, int _k){	
	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen(_fastq_file, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		char *_read = strdup(seq->seq.s);
		char *_read_name = strdup(seq->name.s);
		if(_read == NULL || _read_name == NULL)
			return NULL; 
		struct MPM *qptr = calloc(100, sizeof(struct MPM));
		if(qptr == NULL)
			return NULL;
		int num = find_MPMs_on_read(qptr, _read, _read_name, _k);
		free(qptr);
	}
	return NULL;
}

int predict_main(char *fasta_file, char *fastq_file){
	/* load kmer hash table in the memory */
	assert(fastq_file != NULL);
	assert(fasta_file != NULL);

	/* load kmer_uthash table */
	char *index_file = concat(fasta_file, ".index");
	if(index_file == NULL)
		return -1;
	
	printf("loading kmer uthash table ...\n");
	int k;
	KMER_HT = kmer_uthash_load(index_file, &k);	
	if(KMER_HT == NULL){
		fprintf(stderr, "Fail to load kmer_uthash table\n");
		exit(-1);		
	}
	printf("k=%d\n", k);
	/* MAX_K is defined in common.h */
	if(k > MAX_K){
		fprintf(stderr, "input k(%d) greater than allowed lenght - 100\n", k);
		exit(-1);		
	}			
	
	/* load fasta_uthash table */
	FASTA_HT = fasta_uthash_load(fasta_file);
	if(FASTA_HT == NULL){
		fprintf(stderr, "Fail to load fasta_uthash table\n");
		exit(-1);		
	}
	//kmer_uthash_display(KMER_HT);
	construct_BAG(fastq_file, k);
	
	kmer_uthash_destroy(KMER_HT);	
	fasta_uthash_destroy(FASTA_HT);	
	return 0;
}

