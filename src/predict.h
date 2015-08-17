/*--------------------------------------------------------------------*/
/* predict.c                                                          */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Created Date: 08-07-2015                                           */
/* Predict Gene Fusion by RNA-seq data.                               */
/*--------------------------------------------------------------------*/
#ifndef _PREDICT_H
#define _PREDICT_H

#include <stdio.h>  
#include <stdlib.h>  
#include <string.h> 
#include <zlib.h>  
#include <assert.h>
#include <math.h>
#include <regex.h>
#include "kseq.h"
#include "alignment.h"
#include "kmer_uthash.h"
#include "bag.h"
#include "fasta_uthash.h"
#include "utils.h"
#include "uthash.h"

#define MAX_KMER_LEN                40
#define MIN_KMER_LEN                10
#define MIN_KMER_MATCH              1
#define MIN_EDGE_WEIGHT             1
#define MIN_HITS                    1
#define MIN_ALIGN_SCORE             0
#define MAX_ALIGN_SCORE             1
#define EPSILON                     0.1
//gene_t
typedef struct {
	char* name; // gap open
	int len;
	int exon_num;
	int hits;
    UT_hash_handle hh;
} gene_t;

//opt
typedef struct {
	char* fq1; // gap open
	char* fq2; // gap open
	char* fa; // gap open
	int k; // gap extension
	int min_kmer_match; // match
	int min_edge_weight; // unmatch
	int match;
	int mismatch;
	int gap;
	int extension;
	int jump_gene;
	int jump_exon;
	int min_hits;
	int seed_len;
	int max_mismatch;
	double min_align_score;
} opt_t;

static          fasta_t   *EXON_HT     = NULL;  // stores sequences in in.fa
static           kmer_t   *KMER_HT     = NULL;  // kmer hash table by indexing in.fa
static            bag_t   *BAGR_HT     = NULL;  // Breakend Associated Graph (BAG)
static           gene_t   *GENE_HT     = NULL; 
static  solution_pair_t   *SOLU_HT     = NULL;  // alignment solition of reads against JUN0_HT
static  solution_pair_t   *SOLU_UNIQ_HT     = NULL;  // alignment solition of reads against JUN0_HT


static inline opt_t *opt_init(){
	opt_t *opt = mycalloc(1, opt_t);
	opt->fq1 = NULL;
	opt->fq2 = NULL;
	opt->fa = NULL;
	opt->k = 15;
	opt->min_kmer_match = 10;
	opt->min_edge_weight = 3;	
	opt->match = 2;
	opt->mismatch = -2.0;
	opt->gap = -5.0;
	opt->extension = -1.0;
	opt->jump_gene = -10.0;
	opt->jump_exon = -8.0;
	opt->min_hits = 3;
	opt->min_align_score = 0.8;
	opt->seed_len = 20;
	opt->max_mismatch = 2;
	return opt;
}
static inline void destory_opt(opt_t *opt){
	if(opt->fq1) free(opt->fq1);
	if(opt->fq2) free(opt->fq2);
	if(opt->fa)  free(opt->fa);
	free(opt);
}

static inline gene_t *gene_init(){
	gene_t *instance = mycalloc(1, gene_t);
	instance->name = NULL;
	instance->exon_num = 0;
	instance->hits = 0;
	instance->len = 0;
	return instance;
}

static inline gene_t *find_gene(gene_t *gene, char *quary){
	if(quary==NULL) return NULL;
	gene_t *s = NULL;
	HASH_FIND_STR(gene, quary, s);
	return s;
}

static inline int gene_display(gene_t *instance){
	if(instance==NULL) return -1;
	gene_t *gene_cur;
	for(gene_cur=instance; gene_cur!=NULL; gene_cur=gene_cur->hh.next){
		printf("name=%s\texon_num=%d\thits=%d\tl=%d\n", gene_cur->name, gene_cur->exon_num, gene_cur->hits, gene_cur->len);
	}
	return 0;
}
/*
 * destory the gene_t
 */
static inline int gene_destory(gene_t **instance){
	if(*instance==NULL) return -1;
	gene_t *gene_cur, *gene_tmp;
	HASH_ITER(hh, *instance, gene_cur, gene_tmp) {
	    HASH_DEL(*instance, gene_cur);  /* delete; users advances to next */
		free(gene_cur);            /* optional- if you want to free  */
	}
	return 0;
}
/*
 * Description:
 *------------
 * index input sequences by kmer hash table

 * Input: 
 *-------
 * fa        - fasta_t hash table contains sequences to be indexed
 * k         - length of kmer

 * Output: 
 *-------
 * kmer_t hash table that contains kmer and its occurnace positions on input seq.
 */
static kmer_t *kmer_index(fasta_t *fa, int k);

/*
 * Description:
 *------------
 * construct breakend associated graph (BAG) by kmer_hash table and RNA-seq reads

 * Input: 
 *-------
 * kmer_uthash        - kmer hash table returned by kmer_uthash_construct
 * fq1                - 5' to 3' end of read
 * fq2                - the other end of read
 * min_kmer_matches   - min number kmer matches between a gene and read needed 
 * min_edge_weight    - edges in the graph with weight smaller than min_edge_weight will be deleted
 * k                  - length of kmer
 * Output: 
 *-------
 * BAG_uthash object that contains the graph.
 */
static bag_t *bag_construct(kmer_t *kmer_uthash, fasta_t *fasta_ht, gene_t **gene_ht, char* fq1, char* fq2, int min_kmer_match, int min_edge_weight, int k);
/*
 * Description:
 *------------
 * find junction sites by aligning supportive reads to concatnated string of gene1 and gene2

 * Input: 
 *-------
 * bag        - BAG_utash object: breakend associated graph returned by BAG_uthash_construct
 * fa         - fasta_uthash object: input sequence returned by fasta_uthash_load
 * opt        - opt_t object: contains all input parameters

 * Output: 
 *-------
 * junction_t object that contains identified junctions.
 */
static int bag_junction_gen(bag_t **bag, fasta_t *fa, kmer_t *kmer, opt_t *opt);
/*
 * Description:
 *------------
 * construct fused transcript by identified fusion

 * Input: 
 *-------
 * junc_ht        - junction_t object: return by junction_construct
 * exon_ht        - fasta_uthash object: input sequence returned by fasta_uthash_load

 * Output: 
 *-------
 * junction_t object that contains identified junctions with one more property -> transcript.
 */
static junction_t *transcript_construct(junction_t *junc_ht, fasta_t *exon_ht);
/*
 * Description:
 *------------
 * 1) find subset of pairs that contain 20bp junction string by at most 2 mismatches
 * 2) align those reads to constructed transcript returned by transcript_construct

 * Input: 
 *-------
 * junc_ht        - junction_t object: return by transcript_construct
 * opt            - opt_t object: contains all input parameters

 * Output: 
 *-------
 * solution_pair_t object that contains alignment results of all reads.
 */
static int align_reads_to_transcript(solution_pair_t **res, bag_t *bag, opt_t *opt);

/*
 * Description:
 *------------
 * revisit junction sites and score them based on alignment results
 
 * Input: 
 *-------
 * sol              - alignment results returned by align_to_transcript
 * junc             - previously identified junction
 * min_align_score  - min accepted alignment identity, alignment with identity < min_align_score will be filtered
 * junc_str_len     - length of junction string 

 * Output: 
 *-------
 * junction_t object that contains identified junctions with one more property -> transcript.
 */

static junction_t *junction_score(solution_pair_t *sol, junction_t *junc, double min_align_score, int junc_str_len);
/*
 * usage info
 */
static int pred_usage(opt_t *opt);
/*
 * main function, called by main.c
 */
int predict(int argc, char *argv[]);

#endif