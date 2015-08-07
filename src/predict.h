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
#include "kmer_uthash.h"
#include "BAG_uthash.h"
#include "fasta_uthash.h"
#include "utils.h"
#include "alignment.h"
#include "junction.h"

#define MAX_ALLOWED_K                     50 // max allowed kmer length for kmer hash table 

static struct fasta_uthash *EXON_HT     = NULL;  // stores sequences in in.fa
static struct kmer_uthash  *KMER_HT     = NULL;  // kmer hash table by indexing in.fa
static struct BAG_uthash   *BAGR_HT     = NULL;  // Breakend Associated Graph (BAG)
static        junction_t   *JUN0_HT     = NULL;  // rough junctions identified from BAG
static   solution_pair_t   *SOLU_HT     = NULL;  // alignment solition of reads against JUN0_HT
static        junction_t   *JUN1_HT     = NULL;  // refined junctions with scores

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

/*
 * Description:
 *------------
 * index input sequences by kmer hash table

 * Input: 
 *-------
 * tb        - fasta_uthash object contains sequences to be indexed
 * k         - length of kmer

 * Output: 
 *-------
 * kmer_uthash object that contains kmer and its occurnace positions on input seq.
 */
static struct kmer_uthash *kmer_uthash_construct(struct fasta_uthash *tb, int k);

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
static struct BAG_uthash *BAG_uthash_construct(struct kmer_uthash *kmer_uthash, char* fq1, char* fq2, int min_kmer_match, int min_edge_weight, int k);
static int find_junction_one_edge(struct BAG_uthash *eg, struct fasta_uthash *fasta_u, opt_t *opt, junction_t **ret);
static junction_t *junction_construct(struct BAG_uthash *, struct fasta_uthash *, opt_t *);
static char *concat_exons(char* _read, struct fasta_uthash *fa_ht, struct kmer_uthash *kmer_ht, int _k, char *gname1, char* gname2, char** ename1, char** ename2, int *junction);
static solution_pair_t *align_to_transcript(junction_t *junc, opt_t *opt);
static int junction_rediscover_unit(junction_t *junc, opt_t *opt, solution_pair_t **sol_pair);
static junction_t *transcript_construct(junction_t *junc_ht, struct fasta_uthash *exon_ht);
static junction_t *junction_score(solution_pair_t *sol, junction_t *junc, opt_t *opt);
/*
 * usage info
 */
static int pred_usage(opt_t *opt);
/*
 * main function
 */
int predict(int argc, char *argv[]);

#endif