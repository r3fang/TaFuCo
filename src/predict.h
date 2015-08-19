/*--------------------------------------------------------------------*/
/* predict.c                                                          */
/* Author: Rongxin Fang                                               */
/* Contact: r3fang@ucsd.edu                                           */
/* Created Date: 08-07-2015                                           */
/* Targeted Gene Fusion Calling (tfc) from RNA-seq data.              */
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

//define input parameter valid range 
#define MAX_KMER_LEN                40
#define MIN_KMER_LEN                10
#define MIN_MIN_KMER_MATCH          1
#define MIN_MIN_EDGE_WEIGHT         1
#define MIN_MIN_HITS                1
#define MIN_MIN_ALIGN_SCORE         0
#define MAX_MIN_ALIGN_SCORE         1
#define EPSILON                     0.1

//gene_t
typedef struct {
	char* name;
	int len;
	int exon_num;
	int hits;
    UT_hash_handle hh;
} gene_t;

//opt
typedef struct {
	char* fq1; 
	char* fq2; 
	char* fa; 
	int k;
	int min_kmer_match; 
	int min_edge_weight;
	int match;
	int mismatch;
	int gap;
	int extension;
	int jump_gene;
	int jump_exon;
	int min_hits;
	int seed_len;
	int max_mismatch;
	int alpha;
	double min_align_score;
	double pvalue;
} opt_t;

/* global variables */
static          fasta_t   *EXON_HT     = NULL;  // stores sequences in in.fa
static           kmer_t   *KMER_HT     = NULL;  // kmer hash table by indexing in.fa
static            bag_t   *BAGR_HT     = NULL;  // Breakend Associated Graph (BAG)
static           gene_t   *GENE_HT     = NULL; 
static  solution_pair_t   *SOLU_HT     = NULL;  // alignment solition of reads against JUN0_HT
static  solution_pair_t   *SOLU_UNIQ_HT     = NULL;  // alignment solition of reads against JUN0_HT

/* intitlize opt_t object */
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
	opt->pvalue=0.05;
	opt->alpha=3;
	return opt;
}

/* destory opt_t object */
static inline void destory_opt(opt_t *opt){
	if(opt->fq1) free(opt->fq1);
	if(opt->fq2) free(opt->fq2);
	if(opt->fa)  free(opt->fa);
	free(opt);
}

/* intitlize gene_t object */
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

/* destory the gene_t */
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
 * usage info
 */
static int pred_usage(opt_t *opt);
/*
 * main function, called by main.c
 */
int predict(int argc, char *argv[]);

#endif