/*--------------------------------------------------------------------*/
/* predict.c                                                          */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Predict Gene Fusion by given fastq files.                          */
/*--------------------------------------------------------------------*/

#include <stdio.h>  
#include <stdlib.h>  
#include <string.h> 
#include <zlib.h>  
#include <assert.h>
#include "uthash.h"
#include "kseq.h"
#include "common.h"
KSEQ_INIT(gzFile, gzread);

/*--------------------------------------------------------------------*/
/*Global paramters.*/

static const char *pcPgmName="predict.c";
static struct kmer_uthash *KMER_HT = NULL;
static struct fasta_uthash *FASTA_HT = NULL;


/*--------------------------------------------------------------------*/

/* Find next Maximal Extended Kmer Match (MEKM) on _read at pos_read. */
char* 
find_next_MEKM(char *_read, int pos_read, int k){
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
/*--------------------------------------------------------------------*/

/* Find all Maximal Extended Kmer Matchs (MEKMs) on _read.            */
size_t 
find_all_MEKMs(char **hits, char* _read, int _k){
	int _read_pos = 0;
	size_t _i = 0; 
	while(_read_pos<(strlen(_read)-_k)){
		char* _exon = find_next_MEKM(_read, _read_pos, _k);
		_read_pos += 1;
		if (_exon != NULL){
			hits[_i] = strdup(_exon);
			free(_exon);
			_i ++;
		}
	}
	return _i;
}

/*--------------------------------------------------------------------*/

/* Construct breakend associated graph by given fq files.             */
void
construct_BAG(char *fq_file1, char *fq_file2, int _k){
	gzFile fp1, fp2;
	kseq_t *seq1, *seq2;
	int l1, l2;
	fp1 = gzopen(fq_file1, "r");
	fp2 = gzopen(fq_file2, "r");
	seq1 = kseq_init(fp1);
	seq2 = kseq_init(fp2);
	while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2))) {
		char *_read1 = strdup(seq1->seq.s);
		char *_read2 = strdup(seq2->seq.s);
		char *_read_name1 = strdup(seq1->name.s);
		char *_read_name2 = strdup(seq2->name.s);
		if(_read1 == NULL || _read_name1 == NULL || _read2 == NULL || _read_name2 == NULL)
			return NULL;
		if(strcmp(_read_name1, _read_name2) != 0){
			fprintf(stderr, "ERROR: %s and %s read name not matching\n", fq_file1, fq_file2);
			exit(-1);					
		}				
		char** hits1 = malloc(strlen(_read1) * sizeof(char*));  
		char** hits2 = malloc(strlen(_read2) * sizeof(char*));  
		
		size_t num1 = find_all_MEKMs(hits1, _read1, _k);
		size_t num2 = find_all_MEKMs(hits2, _read2, _k);
		
		free(_read1);
		free(_read2);
		free(_read_name1);
		free(_read_name2);
		free(hits1);
		free(hits2);
	}
	kseq_destroy(seq1);
	kseq_destroy(seq2);	
	gzclose(fp1);
	gzclose(fp2);
	return NULL;
}

/*--------------------------------------------------------------------*/

/* main function. */
int 
predict_main(char *fasta_file, char *fq_file1, char *fq_file2){
	/* load kmer hash table in the memory */
	assert(fasta_file != NULL);
	assert(fq_file1 != NULL);
	assert(fq_file2 != NULL);

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
	
	construct_BAG(fq_file1, fq_file2, k);
	
	kmer_uthash_destroy(KMER_HT);	
	fasta_uthash_destroy(FASTA_HT);	
	return 0;
}

