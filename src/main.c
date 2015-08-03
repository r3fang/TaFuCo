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

static inline solution_pair_t* align_edge(struct BAG_uthash*, struct fasta_uthash *, opt_t *);
static inline int junction_edge(solution_pair_t*, char*, junction_t**, double, double);
static inline junction_t *junction_gen(struct BAG_uthash *, struct fasta_uthash *, opt_t *);
static inline ref_t *ref_generate(struct fasta_uthash *, char*, char*);

static inline ref_t 
*ref_generate(struct fasta_uthash *tb, char* gname1, char* gname2){
	if(tb==NULL || gname1==NULL || gname2==NULL) die("[%s] input error", __func__);
	ref_t *ref = ref_init();
	struct fasta_uthash *s, *tmp;
	register int i = 0;	
	register char* str_tmp;
	// count the number of exons
	HASH_ITER(hh, tb, s, tmp) {
		if(strcmp(strsplit(s->name, '.')[0], gname1) == 0) ref->S1_l += 2;
		if(strcmp(strsplit(s->name, '.')[0], gname2) == 0) ref->S2_l += 2;
	}
	ref->S1 = mycalloc(ref->S1_l, int);
	ref->S2 = mycalloc(ref->S2_l, int);
	
	HASH_ITER(hh, tb, s, tmp) {
		if(s!=NULL){
			if(strcmp(strsplit(s->name, '.')[0], gname1) == 0){				
				if(ref->s==NULL){
					ref->s = strdup(s->seq);
					ref->S1[i++] = 0;
					ref->S1[i++] = strlen(ref->s);
				}else{
					ref->S1[i++] = strlen(ref->s);
					str_tmp = concat(ref->s, s->seq);
					free(ref->s); ref->s=strdup(str_tmp);
					free(str_tmp);
					ref->S1[i++] = strlen(ref->s);
				}			
			}
		}
	}
	i = 0;
	ref->J = strlen(ref->s);
	HASH_ITER(hh, tb, s, tmp) {
		if(s!=NULL){
			if(strcmp(strsplit(s->name, '.')[0], gname2) == 0){				
				if(ref->s==NULL){
					ref->s = strdup(s->seq);
					ref->S2[i++] = 0;
					ref->S2[i++] = strlen(ref->s);
				}else{
					ref->S2[i++] = strlen(ref->s);
					str_tmp = concat(ref->s, s->seq);
					free(ref->s); ref->s=strdup(str_tmp);
					free(str_tmp);
					ref->S2[i++] = strlen(ref->s);
				}			
			}
		}
	}
	ref->l = strlen(ref->s);
	return ref;
}

/*
 * generate junction sites from breakend associated graph 
 */
static inline junction_t 
*junction_gen(struct BAG_uthash *tb, struct fasta_uthash *fa, opt_t *opt){
	if(tb == NULL || fa == NULL || opt==NULL) die("[%s] input error", __func__);	
	struct BAG_uthash *edge, *tmp_bag;
	register int i;
	solution_pair_t *p;
	junction_t *cur_junction, *tmp_junction, *ret = NULL;
	HASH_ITER(hh, tb, edge, tmp_bag){ // iterate every edge
		fprintf(stderr, "edge=%s \n", edge->edge);
		if(((p = align_edge(edge, fa, opt)))==NULL) return NULL;
		if(junction_edge(p, edge->edge, &ret, opt->min_align_score, opt->min_hits) != 0) return NULL;		
	}
	return ret;
}

/*
 *
 * Align reads that support e(Vi, Vj) to the string concated by exons of 
 * Vi and Vj.
 *
 */
static inline solution_pair_t* 
align_edge(struct BAG_uthash *eg, struct fasta_uthash *fasta_u, opt_t *opt){
	int _k = 15;
	char* gname1 = strsplit(eg->edge, '_')[0];
	char* gname2 = strsplit(eg->edge, '_')[1];
	
	ref_t *ref1 = ref_generate(fasta_u, gname1, gname2);
	ref_t *ref2 = ref_generate(fasta_u, gname2, gname1);
	// because we don't know the order of gene fusion
	// therefore, we need align E1, E2 to both order 
	// and decide the order of gene fusion based on alignment score
	solution_pair_t *sol_pairs_r1 = NULL;
	solution_pair_t *sol_pairs_r2 = NULL;
	solution_pair_t *s, *tmp; 
	register int i, j, m, n;
	char* idx_r1, *idx_r2;
	register struct kmer_uthash *s_kmer;
	register char buff[_k];
	for(i=0; i<eg->weight; i++){
		//printf("%d\t%d\n", i, eg->weight);
		//char* aaa = strsplit(eg->evidence[i], '_')[0];
		//for(m=0; m<strlen(aaa)-_k+1; m++){
		//	strncpy(buff, aaa + m, _k); buff[_k] = '\0';	
		//	if((s_kmer=find_kmer(KMER_HT, buff)) == NULL) {continue;}	
		//	for(n=0;n<s_kmer->count;n++) printf("%s\t", s_kmer->seq_names[n]);
		//	printf("\n");
		//}
		//
		//printf("%s\n", aaa);
		solution_t *a = align(strsplit(eg->evidence[i], '_')[0], ref1, opt);
		solution_t *b = align(strsplit(eg->evidence[i], '_')[1], ref1, opt);
		solution_t *c = align(strsplit(eg->evidence[i], '_')[0], ref2, opt);
		solution_t *d = align(strsplit(eg->evidence[i], '_')[1], ref2, opt);
		idx_r1 = idx_md5(eg->edge, a->pos, b->pos);
		idx_r2 = idx_md5(eg->edge, c->pos, d->pos);
		// sol_pairs_r1
		HASH_FIND_STR(sol_pairs_r1, idx_r1, s);
		if(s==NULL){
			s = solution_pair_init();
			s->idx = strdup(idx_r1);
			s->r1 = a; s->r2 = b;
			s->prob = a->prob * b->prob;
			HASH_ADD_STR(sol_pairs_r1, idx, s);  /* idx: name of key field */
		}else{ // update with higher score
			if(s->prob < a->prob * b->prob){
				s->prob =  a->prob * b->prob;
				s->r1 = a; s->r2 = b;
			}
		}
		// sol_pairs_r2
		HASH_FIND_STR(sol_pairs_r2, idx_r2, s);
		if(s==NULL){
			s = solution_pair_init();
			s->idx = strdup(idx_r2);
			s->r1 = c; s->r2 = d;
			s->prob = c->prob * d->prob;
			HASH_ADD_STR(sol_pairs_r2, idx, s);  /* idx: name of key field */
		}else{ // update with higher score
			if(s->prob < c->prob * d->prob){
				s->prob = c->prob * d->prob;
				s->r1 = c; s->r2 = d;
			}
		}
	}
	if(idx_r1) free(idx_r1);
	if(idx_r2) free(idx_r2);
	/*------------------------------------------------------------------------------*/	
	// make decision of gene order, chose the one with larger likelihood
	register double likehood1, likehood2;
	likehood1 = likehood2 = 0;
    HASH_ITER(hh, sol_pairs_r1, s, tmp) {
		likehood1 += 10*log(s->r1->prob);
		likehood1 += 10*log(s->r2->prob);		
    }
    HASH_ITER(hh, sol_pairs_r2, s, tmp) {
		likehood2 += 10*log(s->r1->prob);
		likehood2 += 10*log(s->r2->prob);		
    }
	if(ref1) ref_destory(ref1);  if(ref2) ref_destory(ref2);
	solution_pair_t * ret;
	if(likehood1 >= likehood2){
		ret = sol_pairs_r1;
	}else{
		ret = sol_pairs_r2;
	}
	if(gname1)         free(gname1);     
	if(gname2)         free(gname2);
	return ret;
}

/*
 * generate junction sites from solution_pair_t
 */
static inline int 
junction_edge(solution_pair_t *p, char* name, junction_t **ret, double MIN_ALIGN_SCORE, double MIN_HITS){
	if(p == NULL) die("[%s] input error", __func__);
	solution_pair_t *s, *tmp;
	junction_t *m, *n;
	char* idx;
    HASH_ITER(hh, p, s, tmp) {
		// one read
		if(s->r1->jump == true && s->r1->prob >= MIN_ALIGN_SCORE){
			idx = idx_md5(name, s->r1->jump_start, s->r1->jump_end);
			HASH_FIND_STR(*ret, idx, m);
			if(m==NULL){ // this junction not in ret
				m = junction_init();
				m->idx = strdup(idx);
				m->name = strdup(name);
				m->start = s->r1->jump_start;
				m->end = s->r1->jump_end;
				m->hits = 1;
				m->likehood = 10*log(s->r1->prob); 				
				memcpy( m->s, &s->r1->s2[m->start-HALF_JUNCTION_LEN-1], HALF_JUNCTION_LEN);
				memcpy( &m->s[HALF_JUNCTION_LEN], &s->r1->s2[m->end], HALF_JUNCTION_LEN);
				HASH_ADD_STR(*ret, idx, m);
			}else{
				m->hits ++;
				m->likehood += 10*log(s->r1->prob); 
			}
		}
		// the other read
		if(s->r2->jump == true && s->r2->prob >= MIN_ALIGN_SCORE){
			idx = idx_md5(name, s->r2->jump_start, s->r2->jump_end);			
			HASH_FIND_STR(*ret, idx, m);
			if(m==NULL){ // this junction not in ret
				m = junction_init();
				m->idx = strdup(idx);
				m->name = strdup(name);
				m->start = s->r2->jump_start;
				m->end = s->r2->jump_end;
				m->hits = 1;
				m->likehood = 10*log(s->r2->prob); 
				memcpy( m->s, &s->r2->s2[m->start-HALF_JUNCTION_LEN-1], HALF_JUNCTION_LEN);
				memcpy( &m->s[HALF_JUNCTION_LEN], &s->r2->s2[m->end], HALF_JUNCTION_LEN);
				HASH_ADD_STR(*ret, idx, m);
			}else{
				m->hits ++;
				m->likehood += 10*log(s->r2->prob); 
			}
		}	
    }
	if(idx) free(idx);
	// delete those junctions with hits < MIN_HITS
	HASH_ITER(hh, *ret, m, n){
		if(m != NULL){
			m->likehood = m->likehood/m->hits;
			if(m->hits < MIN_HITS){
				HASH_DEL(*ret,m);
				free(m);
			}
		}
	}
	return 0;
}



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
