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

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.30-r15"
#endif

/* error code */
#define PR_ERR_NONE		     	 		0 		// no error

/*--------------------------------------------------------------------*/
/*Global paramters.*/
static struct kmer_uthash *KMER_HT      = NULL;
static struct fasta_uthash *FASTA_HT    = NULL;
static struct BAG_uthash *BAG_HT        = NULL;

//opt
typedef struct {
	char* fq1; // gap open
	char* fq2; // gap open
	char* fa; // gap open
	int k; // gap extension
	int min_match; // match
	int min_weight; // unmatch
	int match;
	int mismatch;
	int gap;
	int extension;
	int jump_gene;
	int jump_exon;
	int min_hits;
	double min_align_score;
} opt_t;

opt_t *init_opt(){
	opt_t *opt = mycalloc(1, opt_t);
	opt->fq1 = NULL;
	opt->fq2 = NULL;
	opt->fa = NULL;
	opt->k = 15;
	opt->min_match = 10;
	opt->min_weight = 3;	
	opt->match = 2;
	opt->mismatch = -2;
	opt->gap = -5;
	opt->extension = -1;
	opt->jump_gene = -10;
	opt->jump_exon = -8;
	opt->min_hits = 3;
	opt->min_align_score = 0.8;
	return opt;
}

/* 
 * Find all genes uniquely matched with kmers on _read.          
 * hash     - a hash table count number of matches between _read and every gene
 * _read    - inqury read
 * _k       - kmer length
 */
int
find_all_matches(str_ctr **hash, char* _read, int _k){
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
construct_BAG(char *fq_file1, char *fq_file2, int _k, int _min_match, struct BAG_uthash **bag){
	if(fq_file1 == NULL || fq_file2 == NULL) die("construct_BAG: parameter error\n");
	int error;
	gzFile fp1, fp2;
	register int l1, l2;
	register kseq_t *seq1, *seq2;
	register char *_read1, *_read2;
	register char *edge_name;
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
		find_all_matches(&gene_counter, _read1, _k);
		find_all_matches(&gene_counter, _read2, _k);
		
		HASH_SORT(gene_counter, str_ctr_sort);
		unsigned int num = HASH_COUNT(gene_counter);
		
		char** hits = mycalloc(num, char*);
		int i=0; if(num > 1){
			HASH_ITER(hh, gene_counter, s, tmp) { 
				if(s->SIZE >= _min_match){hits[i++] = strdup(s->KEY);}
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
	if (optind + 3 > argc) {
			fprintf(stderr, "\n");
					fprintf(stderr, "Usage:   tfc [options] <in.fa> <R1.fq> <R2.fq>\n\n");
					fprintf(stderr, "Options: --------------------  BAG Graph Options  -------------------\n");
					fprintf(stderr, "         -n INT   min number kmer match for a hit between read and gene [%d]\n", opt->min_match);
					fprintf(stderr, "         -w INT   min weight of edge on BAG [%d]\n", opt->min_weight);
					fprintf(stderr, "         -k INT   kmer length for indexing reference [%d]\n", opt->k);
					fprintf(stderr, "\n");
					fprintf(stderr, "         --------------------  Alignment Options  --------------------\n");
					fprintf(stderr, "         -m INT   match score [%d]\n", opt->match);
					fprintf(stderr, "         -u INT   mismatch score [%d]\n", opt->mismatch);
					fprintf(stderr, "         -o INT   gap open penality [%d]\n", opt->gap);
					fprintf(stderr, "         -e INT   gap extension penality [%d]\n", opt->extension);
					fprintf(stderr, "         -j INT   jump penality between genes [%d]\n", opt->jump_gene);
					fprintf(stderr, "         -s INT   jump penality between exons [%d]\n", opt->jump_exon);
					fprintf(stderr, "\n");
					fprintf(stderr, "         --------------------  Junction Options  --------------------\n");										
					fprintf(stderr, "         -h INT   min number of hits for a junction [%d]\n", opt->min_hits);					
					fprintf(stderr, "         -a FLOAT min alignment score [%f]\n", opt->min_align_score);
					fprintf(stderr, "\n");
					return 1;
	}
	opt->fa = argv[argc-3];
	opt->fq1 = argv[argc-2];
	opt->fq2 = argv[argc-1];
	/* load kmer hash table in the memory */
	/* load kmer_uthash table */
	fprintf(stderr, "[%s] Generating kmer hash table (K=%d) ... \n",__func__, opt->k);
	KMER_HT = kmer_uthash_construct(opt->fa, opt->k);	
	if(KMER_HT == NULL) die("Fail to load the index\n");	
	///* load fasta_uthash table */
	fprintf(stderr, "[%s] Loading fasta hash table ... \n", __func__);
	// load fasta sequences
	if((fasta_uthash_load(opt->fa, &FASTA_HT)) != PR_ERR_NONE) die("main: fasta_uthash_load fails\n");	
	// construct break-end associated graph
	fprintf(stderr, "[%s] constructing break-end associated graph ... \n", __func__);
	if((construct_BAG(opt->fq1, opt->fq2, opt->k, opt->min_match, &BAG_HT)) != PR_ERR_NONE)	die("main: construct_BAG fails\n");		
	// rm duplicate evidence for edges on BAG
	if((BAG_uthash_uniq(&BAG_HT)) != PR_ERR_NONE) die("main: BAG_uthash_uniq fails\n");
	// delete edges with weight < opt->min_weight
	if(BAG_uthash_trim(&BAG_HT, opt->min_weight) != PR_ERR_NONE)	die("main: BAG_uthash_trim\n");		
	if(BAG_uthash_display(BAG_HT) != PR_ERR_NONE)	die("main: BAG_uthash_trim\n");		
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
