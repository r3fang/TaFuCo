/*--------------------------------------------------------------------*/
/* main.c                                                             */
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
#include "kmer_uthash.h"
#include "BAG_uthash.h"
#include "fasta_uthash.h"
#include "utils.h"
#include "alignment.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.30-r15"
#endif

/* error code */
#define PR_ERR_NONE		     	 		0 		// error code
/*--------------------------------------------------------------------*/
/*Global paramters.*/
static struct kmer_uthash  *KMER_HT       = NULL;
static struct fasta_uthash *FASTA_HT      = NULL;
static struct BAG_uthash   *BAG_HT        = NULL;
static        junction_t   *JUNCTION_HT   = NULL;

static inline int tfc_usage(opt_t *opt){
	fprintf(stderr, "\n");
			fprintf(stderr, "Usage:   tfc [options] <in.fa> <R1.fq> <R2.fq>\n\n");
			fprintf(stderr, "Options: --------------------   Graph Options  -----------------------\n");
			fprintf(stderr, "         -k INT   kmer length [%d]\n", opt->k);
			fprintf(stderr, "         -n INT   min number kmer matches [%d]\n", opt->min_match);
			fprintf(stderr, "         -w INT   min weight for an edge [%d]\n", opt->min_weight);
			fprintf(stderr, "\n");
			fprintf(stderr, "         --------------------  Alignment Options  --------------------\n");
			fprintf(stderr, "         -m INT   match score [%d]\n", opt->match);
			fprintf(stderr, "         -u INT   penality for mismatch [%d]\n", opt->mismatch);
			fprintf(stderr, "         -o INT   penality for gap open [%d]\n", opt->gap);
			fprintf(stderr, "         -e INT   penality for gap extension [%d]\n", opt->extension);
			fprintf(stderr, "         -j INT   penality for jump between genes [%d]\n", opt->jump_gene);
			fprintf(stderr, "         -s INT   penality for jump between exons [%d]\n", opt->jump_exon);
			fprintf(stderr, "\n");
			fprintf(stderr, "         --------------------  Junction Options  --------------------\n");										
			fprintf(stderr, "         -h INT   min number of hits for a junction to be called [%d]\n", opt->min_hits);					
			fprintf(stderr, "         -a FLOAT min alignment score [%.2f]\n", opt->min_align_score);
			fprintf(stderr, "\n");
			return 1;
}
/*--------------------------------------------------------------------*/
/* main function. */
int main(int argc, char *argv[]) {
	opt_t *opt = init_opt(); // initlize options with default settings
	int c, i;
	srand48(11);
	while ((c = getopt(argc, argv, "m:w:k:n:u:o:e:g:s")) >= 0) {
				switch (c) {
				case 'n': opt->min_match = atoi(optarg); break;
				case 'w': opt->min_weight = atoi(optarg); break;
				case 'k': opt->k = atoi(optarg); break;
				case 'm': opt->match = atoi(optarg); break;
				case 'u': opt->mismatch = atoi(optarg); break;
				case 'o': opt->gap = atoi(optarg); break;
				case 'e': opt->extension = atoi(optarg); break;
				case 'j': opt->jump_gene = atoi(optarg); break;
				case 's': opt->jump_exon = atoi(optarg); break;
				case 'h': opt->min_hits = atoi(optarg); break;
				case 'a': opt->min_align_score = atof(optarg); break;
				default: return 1;
		}
	}
	if (optind + 3 > argc) return tfc_usage(opt);
	opt->fa = argv[optind];
	opt->fq1 = argv[optind+1];
	opt->fq2 = argv[optind+2];
	/* load kmer hash table in the memory */
	/* load kmer_uthash table */
	fprintf(stderr, "[%s] generating kmer hash table (K=%d) ... \n",__func__, opt->k);
	KMER_HT = kmer_uthash_construct(opt->fa, opt->k);	
	if(KMER_HT == NULL) die("Fail to load the index\n");	
	///* load fasta_uthash table */
	fprintf(stderr, "[%s] loading fasta hash table ... \n", __func__);
	// load fasta sequences
	if((FASTA_HT=fasta_uthash_load(opt->fa)) == NULL) die("main: fasta_uthash_load fails\n");	
	// construct break-end associated graph
	fprintf(stderr, "[%s] constructing break-end associated graph ... \n", __func__);
	if((construct_BAG(&BAG_HT, KMER_HT, opt->fq1, opt->fq2, opt->min_match, opt->k)) != PR_ERR_NONE)	die("main: construct_BAG fails\n");		
	// delete edges with weight < opt->min_weight
	if(BAG_uthash_trim(&BAG_HT, opt->min_weight) != PR_ERR_NONE)	die("main: BAG_uthash_trim\n");		
	
	fprintf(stderr, "[%s] identifying junction sites ... \n", __func__);
	JUNCTION_HT = junction_gen(BAG_HT, FASTA_HT, opt);
	junction_t *cur_junction, *tmp_junction;
	HASH_ITER(hh, JUNCTION_HT, cur_junction, tmp_junction) {
		printf("name=%s: start=%d\tend=%d\thits=%zu\tlikelihood=%f\nstr=%s\n", cur_junction->name, cur_junction->start, cur_junction->end, cur_junction->hits,cur_junction->likehood,cur_junction->s);
	}
	//*--------------------------------------------------------------------*/	
	// clear up the masses
	if(kmer_uthash_destroy(&KMER_HT)   != PR_ERR_NONE)	die("main: kmer_uthash_destroy\n");	
	if(fasta_uthash_destroy(&FASTA_HT) != PR_ERR_NONE)	die("main: fasta_uthash_destroy fails\n");		
	if(BAG_uthash_destroy(&BAG_HT)     != PR_ERR_NONE)	die("main: BAG_uthash_destroy\n");	
	/*--------------------------------------------------------------------*/	
	fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
	fprintf(stderr, "[%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n");
	return 0;
}
