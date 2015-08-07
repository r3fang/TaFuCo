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

static int find_junction_one_edge(struct BAG_uthash *eg, struct fasta_uthash *fasta_u, opt_t *opt, junction_t **ret);
static junction_t *junction_construct(struct BAG_uthash *, struct fasta_uthash *, opt_t *);
static char *concat_exons(char* _read, struct fasta_uthash *fa_ht, struct kmer_uthash *kmer_ht, int _k, char *gname1, char* gname2, char** ename1, char** ename2, int *junction);
static struct kmer_uthash *kmer_uthash_construct(struct fasta_uthash *tb, int k);
static struct BAG_uthash *BAG_uthash_construct(struct kmer_uthash *kmer_uthash, char* fq1, char* fq2, int _min_match, int _k);
static solution_pair_t *align_to_transcript(junction_t *junc, opt_t *opt);
static int junction_rediscover_unit(junction_t *junc, opt_t *opt, solution_pair_t **sol_pair);
static junction_t *transcript_construct(junction_t *junc_ht, struct fasta_uthash *exon_ht);
static junction_t *junction_score(solution_pair_t *sol, junction_t *junc, opt_t *opt);
static int pred_usage(opt_t *opt);
int predict(int argc, char *argv[]);

#endif