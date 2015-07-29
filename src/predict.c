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
#include <math.h>
#include <regex.h>
#include "kseq.h"
#include "common.h"
#include "kstring.h"
#include "uthash.h"
#include "kmer_uthash.h"
#include "BAG_uthash.h"
#include "fasta_uthash.h"
#include "utils.h"

/* error code */
#define PR_ERR_NONE		     	 		0 		// no error
/*--------------------------------------------------------------------*/
/*Global paramters.*/
static struct kmer_uthash *KMER_HT  	= NULL;
static struct fasta_uthash *FASTA_HT 	= NULL;
static struct BAG_uthash *BAG_HT 		= NULL;
/*--------------------------------------------------------------------*/
#define k                               20

typedef struct
{
	char 	*KEY;
	size_t  SIZE;
	UT_hash_handle hh;
} str_ctr;

int str_ctr_add(str_ctr** tb, char* key){
	if(key == NULL) die("str_ctr_add: prameter error\n");
	str_ctr *s;
	HASH_FIND_STR(*tb, key, s);
	if(s == NULL){
		s = mycalloc(1, str_ctr);
		s->KEY = key;
		s->SIZE = 1;
		HASH_ADD_STR(*tb, KEY, s);
	}else{
		s->SIZE++;
	}
	return PR_ERR_NONE;
}

int str_ctr_sort(str_ctr *a, str_ctr *b) {
    return (a->SIZE >= b->SIZE);
}

/* 
 * Find all genes uniquely matched with kmers on _read.          
 * hash     - a hash table count number of matches between _read and every gene
 * _read    - inqury read
 * _k       - kmer length
 */
int
find_all_matched_genes(str_ctr **hash, char* _read, int _k){
/*--------------------------------------------------------------------*/
	/* check parameters */
	if(_read == NULL || _k < 0) die("find_all_MEKMs: parameter error\n");
/*--------------------------------------------------------------------*/
	/* declare vaiables */
	str_ctr *s;
	int _read_pos = 0;
	char* gene = NULL;
	struct kmer_uthash *s_kmer = NULL; 
	char buff[_k];
/*--------------------------------------------------------------------*/
	while(_read_pos<(strlen(_read)-_k-1)){
		/* copy a kmer of string */
		strncpy(buff, _read + _read_pos, _k); buff[_k] = '\0';	
		if(strlen(buff) != _k) die("find_next_match: buff strncpy fails\n");
		/*------------------------------------------------------------*/
		if(find_kmer(KMER_HT, buff, &s_kmer) != PR_ERR_NONE) die("find_next_match: find_kmer fails\n");
		if(s_kmer == NULL){_read_pos++; continue;} // kmer not in table but not an error
		if(s_kmer->count == 1){ // only count the uniq match 
			gene = strdup(s_kmer->seq_names[0]);
			if(gene == NULL) die("find_next_match: get_exon_name fails\n");
			if(str_ctr_add(hash, gene) != PR_ERR_NONE) die("find_all_MEKMs: str_ctr_add fails\n");
		}
		_read_pos++;
	}
	return PR_ERR_NONE;
}

/* 
 * Construct breakend associated graph.             
 * fq_file1 - fastq file name for R1
 * fq_file2 - fastq file name for R2
 * _k       - kmer length for hash table
 * cutoff   - min number of kmer matches between a read again a gene.
 * bag      - BAG_uthash object (breakend associated graph)
 */
int
construct_BAG(char *fq_file1, char *fq_file2, int _k, int cutoff, struct BAG_uthash **bag){
	if(fq_file1 == NULL || fq_file2 == NULL) die("construct_BAG: parameter error\n");
	int error;
	gzFile fp1, fp2;
	kseq_t *seq1, *seq2;
	int l1, l2;
	char *_read1, *_read2, *edge_name;
	_read1 = _read2 = edge_name = NULL;
	
	fp1 = gzopen(fq_file1, "r");
	fp2 = gzopen(fq_file2, "r");
	seq1 = kseq_init(fp1);
	seq2 = kseq_init(fp2);
	
	if(fp1 == NULL || fp2 == NULL || seq1 == NULL || seq2 == NULL) die("construct_BAG: fail to read fastq files\n");

	while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >= 0 ) {
		_read1 = rev_com(seq1->seq.s); // reverse complement of read1
		_read2 = seq2->seq.s;		
		if(_read1 == NULL || _read2 == NULL) die("construct_BAG: fail to get _read1 and _read2\n");
		if(strcmp(seq1->name.s, seq2->name.s) != 0) die("construct_BAG: read pair not matched\n");		
		if(strlen(_read1) < _k || strlen(_read2) < _k){continue;}
		str_ctr* gene_counter = NULL;
		str_ctr *s, *tmp;
		find_all_matched_genes(&gene_counter, _read1, _k);
		find_all_matched_genes(&gene_counter, _read2, _k);
		
		HASH_SORT(gene_counter, str_ctr_sort);
		unsigned int num = HASH_COUNT(gene_counter);
		
		char** hits = mycalloc(num, char*);
		int i=0; if(num > 1){
			HASH_ITER(hh, gene_counter, s, tmp) { 
				if(s->SIZE >= cutoff){hits[i++] = strdup(s->KEY);}
			}			
		}
		int m, n; for(m=0; m < i; m++){for(n=m+1; n < i; n++){
				int rc = strcmp(hits[m], hits[n]);
				if(rc<0)  edge_name = concat(concat(hits[m], "_"), hits[n]);
				if(rc>0)  edge_name = concat(concat(hits[n], "_"), hits[m]);
				if(rc==0) edge_name = NULL;
				if(edge_name!=NULL){
					if(BAG_uthash_add(bag, edge_name, concat(concat(_read1, "_"), _read2)) != PR_ERR_NONE) die("BAG_uthash_add fails\n");							
				}
		}}
		if(hits)		 free(hits);
		if(gene_counter) free(gene_counter);
	}
	if(edge_name)   free(edge_name);
	kseq_destroy(seq1);
	kseq_destroy(seq2);	
	gzclose(fp1);
	gzclose(fp2);
	return PR_ERR_NONE; // no error raised
}

/*--------------------------------------------------------------------*/

/* main function. */
int main(int argc, char *argv[]) {
	if (argc != 6) {  
	        fprintf(stderr, "Usage: %s <in.fa> <read_R1.fq> <read_R2.fq> <int k> <int min_match> <int min_weight>\n", argv[0]);  
	        return 1;  
	 }
 	int  min_match, min_weight;
	char *fasta_file = argv[1];
	char *fq_file1 = argv[2];
	char *fq_file2 = argv[3];
	if (sscanf (argv[4], "%d", &min_match)!=1)	die("Input error: wrong type for k\n");
	if (sscanf (argv[5], "%d", &min_weight)!=1)	die("Input error: wrong type for min_weight\n");
	/* load kmer hash table in the memory */
	///* load kmer_uthash table */
	printf("Generating kmer hash table ... \n");
	KMER_HT = kmer_uthash_construct(fasta_file, k);	
	if(KMER_HT == NULL) die("Fail to load the index\n");	
	///* load fasta_uthash table */
	printf("Loading fasta hash table ... \n");
	if((fasta_uthash_load(fasta_file, &FASTA_HT)) != PR_ERR_NONE) die("main: fasta_uthash_load fails\n");	
	if((construct_BAG(fq_file1, fq_file2, k, min_match, &BAG_HT)) != PR_ERR_NONE)	die("main: construct_BAG fails\n");	
	
	//timeUpdate();
	//if(BAG_uthash_trim(&BAG_HT, min_weight) != PR_ERR_NONE)	die("main: BAG_uthash_trim\n");		
	if(BAG_uthash_display(BAG_HT)   != PR_ERR_NONE)	die("main: kmer_uthash_destroy\n");	
	//timeUpdate();
	///*--------------------------------------------------------------------*/	
	//if(kmer_uthash_destroy(&KMER_HT)   != PR_ERR_NONE)	die("main: kmer_uthash_destroy\n");	
	//if(fasta_uthash_destroy(&FASTA_HT) != PR_ERR_NONE)	die("main: fasta_uthash_destroy fails\n");		
	//if(BAG_uthash_destroy(&BAG_HT)     != PR_ERR_NONE)	die("main: BAG_uthash_destroy\n");	
	return 0;
}