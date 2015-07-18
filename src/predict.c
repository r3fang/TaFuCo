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
#include "uthash.h"
#include "kmer_uthash.h"
#include "BAG_uthash.h"
#include "fasta_uthash.h"
#include "utils.h"

/* error code */
#define PR_ERR_NONE		     	 0 		// no error
/*--------------------------------------------------------------------*/
/*Global paramters.*/
static const char *pcPgmName			="predict.c";
static struct kmer_uthash *KMER_HT  	= NULL;
static struct fasta_uthash *FASTA_HT 	= NULL;
static struct BAG_uthash *BAG_HT 		= NULL;
/*--------------------------------------------------------------------*/

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

char* 
get_gene_name(char* s){
	if(s == NULL) die("get_exon_name: parameter error\n");
	const char delim[2] = ".";
	char *token;
	token = strdup(strtok(s, delim));
	if(token == NULL) die("get_exon_name: strtok fails\n");
	free(s);
	return token;
}
/* Find all Maximal Extended Kmer Matchs (MEKMs) on _read.            */
int
find_all_MEKMs(str_ctr **hash, char* _read, int _k){
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
		if(s_kmer->count == 1){
			gene = get_gene_name(strdup(s_kmer->pos[0]));
			if(gene == NULL) die("find_next_match: get_exon_name fails\n");
			if(str_ctr_add(hash, gene) != PR_ERR_NONE) die("find_all_MEKMs: str_ctr_add fails\n");
		}
		_read_pos++;
	}
	return PR_ERR_NONE;
}

/* Construct breakend associated graph by given fq files.             */
int
construct_BAG(char *fq_file1, char *fq_file2, int _k, int min_match, struct BAG_uthash **graph_ht){
	if(fq_file1 == NULL || fq_file2 == NULL || *graph_ht != NULL) die("construct_BAG: parameter error\n");
	int error;
	gzFile fp1, fp2;
	kseq_t *seq1, *seq2;
	int l1, l2;
	char *_read1, *_read2, *edge_name, *gene1, *gene2;	
	_read1 = _read2 = edge_name = gene1 = gene2 = NULL;
		
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
		//printf("%s\t%s\n", seq1->name.s, seq2->name.s);
		str_ctr* gene_counter = NULL;
		str_ctr *s, *tmp;
		find_all_MEKMs(&gene_counter, _read1, _k);
		find_all_MEKMs(&gene_counter, _read2, _k);
		
		HASH_SORT(gene_counter, str_ctr_sort);
		unsigned int num;
		num = HASH_COUNT(gene_counter);
		if(num > 1){
			HASH_ITER(hh, gene_counter, s, tmp) {
			    printf("KEY=%s: SIZE=%zu\n", s->KEY, s->SIZE);
			}			
		}		
		if(gene_counter) free(gene_counter);
		//printf("%d\t%d\n", num1, num2);
		//size_t size1 = set_str_arr(hits1, hits_uniq1, num1);
		//size_t size2 = set_str_arr(hits2, hits_uniq2, num2);
		//int i, j; for(i=0; i<size1; i++){ for(j=0; j<size2; j++){
		//		size_t len1 = strsplit_size(hits_uniq1[i], ".");
		//		size_t len2 = strsplit_size(hits_uniq2[j], ".");
		//		strsplit(hits_uniq1[i], len1, parts1, ".");  
		//		strsplit(hits_uniq2[j], len2, parts2, ".");  
		//		printf("%s\t%s\n", hits_uniq1[i], hits_uniq1[j]);
		//		if(parts1[0] == NULL || parts2[0] == NULL) continue;
		//		gene1 = strdup(parts1[0]);
		//		gene2 = strdup(parts2[0]);
		//		//printf("%s\t%s\n", gene1, gene2);
		//		int rc = strcmp(gene1, gene2);
		//		if(rc<0)  edge_name = concat(concat(gene1, "_"), gene2);
		//		if(rc>0)  edge_name = concat(concat(gene2, "_"), gene1);
		//		if(rc==0) edge_name = NULL;
		//		//if(edge_name!=NULL){
		//		//	printf("%s\t%s\n", edge_name, concat(concat(_read1, "_"), _read2));
		//			//if((error = BAG_uthash_add(&graph_ht, edge_name, concat(concat(_read1, "_"), _read2))) != PR_ERR_NONE)
		//			//	die("BAG_uthash_add fails\n");					
		//		//}
		//}}
	}
	if(gene1)		free(gene1);
	if(gene2)		free(gene2);
	kseq_destroy(seq1);
	kseq_destroy(seq2);	
	gzclose(fp1);
	gzclose(fp2);
	return PR_ERR_NONE;
}

/*--------------------------------------------------------------------*/

/* main function. */
int main(int argc, char *argv[]) { 
	if (argc != 5) {  
	        fprintf(stderr, "Usage: %s <in.fa> <read_R1.fq> <read_R2.fq> <int k>\n", argv[0]);  
	        return 1;  
	 }
	char *fasta_file = argv[1];
	char *fq_file1 = argv[2];
	char *fq_file2 = argv[3];
	int min_mtch;		
	if (sscanf (argv[4], "%d", &min_mtch)!=1) die("Input error: wrong type for k\n");
	/* load kmer hash table in the memory */
	int error;
	///* load kmer_uthash table */
	char *index_file = concat(fasta_file, ".index");
	if(index_file == NULL) die("Fail to concate index_file\n");	
	
	
	
	int k; if((kmer_uthash_load(index_file, &k, &KMER_HT)) != PR_ERR_NONE) die("main: kmer_uthash_load fails\n");	
	timeUpdate();
    
	if(KMER_HT == NULL) die("Fail to load the index\n");
	if(k > MAX_K) die("input k(%d) greater than allowed lenght - 100\n", k);	
	/* load fasta_uthash table */
	if((error=fasta_uthash_load(fasta_file, &FASTA_HT)) != PR_ERR_NONE) 			   die("main: fasta_uthash_load fails\n");	
	timeUpdate();
    
	if((error=construct_BAG(fq_file1, fq_file2, k, min_mtch, &BAG_HT)) != PR_ERR_NONE) die("main: construct_BAG fails\n");	
	timeUpdate();
	
	if((error=kmer_uthash_destroy(&KMER_HT))   != PR_ERR_NONE) 						   die("main: kmer_uthash_destroy\n");	
	if((error=fasta_uthash_destroy(&FASTA_HT)) != PR_ERR_NONE) 						   die("main: fasta_uthash_destroy fails\n");		
	
	//if((error=BAG_uthash_destroy(&BAG_HT))     != PR_ERR_NONE) 						   die("main: BAG_uthash_destroy\n");	
	if((error=fasta_uthash_display(FASTA_HT)) != PR_ERR_NONE) 			die("main: fasta_uthash_display fails\n");	
	//if((error=BAG_uthash_display(BAG_HT)) != PR_ERR_NONE) 				die("main: BAG_uthash_display fails\n");		
	return 0;
}

