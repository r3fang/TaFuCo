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

#define PR_ERR_NONE                       0
#define MAX_EXON_NUM                      5000
#define EXON_FLANK_LEN                    0

#define pair(k1, k2)  ((k1 + k2)*(k1 + k2 + 1)/2 + k2)


static struct kmer_uthash *KMER_HT      = NULL;
static struct fasta_uthash *FASTA_HT    = NULL;
static struct BAG_uthash *BAG_HT        = NULL;

typedef struct {
	unsigned long idx;
	char *name;
	int start;
	int end;
	char* str;
	double prob;
    UT_hash_handle hh;
} junction_t;

// alingment soulution
typedef struct
{
	unsigned long idx;
	solution *r1;
	solution *r2;
	double prob;
	UT_hash_handle hh;
} solution_pair_t;

void solution_pair_destory(solution_pair_t *s){
	if(s!=NULL){
		if(s->r1) solution_destory(s->r1);
		if(s->r2) solution_destory(s->r2);
	}
	free(s);
}
ref_t* ref_generate(struct fasta_uthash *tb, char* gname1, char* gname2){
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
					ref->S1[i++] = EXON_FLANK_LEN;
					ref->S1[i++] = strlen(ref->s)-EXON_FLANK_LEN;
				}else{
					ref->S1[i++] = strlen(ref->s)+EXON_FLANK_LEN;
					str_tmp = concat(ref->s, s->seq);
					free(ref->s); ref->s=strdup(str_tmp);
					free(str_tmp);
					ref->S1[i++] = strlen(ref->s)-EXON_FLANK_LEN;
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
					ref->S2[i++] = EXON_FLANK_LEN;
					ref->S2[i++] = strlen(ref->s)-EXON_FLANK_LEN;
				}else{
					ref->S2[i++] = strlen(ref->s)+EXON_FLANK_LEN;
					str_tmp = concat(ref->s, s->seq);
					free(ref->s); ref->s=strdup(str_tmp);
					free(str_tmp);
					ref->S2[i++] = strlen(ref->s)-EXON_FLANK_LEN;
				}			
			}
		}
	}
	ref->l = strlen(ref->s);
	return ref;
}
/*
 * Align reads that support e(Vi, Vj) to the string concated by exons of 
 * Vi and Vj.
 */
void edge_align(struct BAG_uthash *eg, struct fasta_uthash *fasta_u){
	char* gname1 = strsplit(eg->edge, '_')[0];
	char* gname2 = strsplit(eg->edge, '_')[1];
	ref_t *ref1 = ref_generate(FASTA_HT, gname1, gname2);
	ref_t *ref2 = ref_generate(FASTA_HT, gname2, gname1);
	// because we don't know the order of gene fusion
	// therefore, we need align E1, E2 to both order 
	// and decide the order of gene fusion based on alignment score
	solution_pair_t *sol_pairs_r1 = NULL;
	solution_pair_t *sol_pairs_r2 = NULL;
	solution_pair_t *s, *tmp; 
	register int i, j;
	int idx_r1, idx_r2;
	
	for(i=0; i<eg->weight; i++){
		solution *a = align(strsplit(eg->evidence[i], '_')[0], ref1);
		solution *b = align(strsplit(eg->evidence[i], '_')[1], ref1);
		solution *c = align(strsplit(eg->evidence[i], '_')[0], ref2);
		solution *d = align(strsplit(eg->evidence[i], '_')[1], ref2);
		idx_r1 = pair(a->pos, b->pos); idx_r2 = pair(c->pos, d->pos);
		// sol_pairs_r1
		HASH_FIND_INT(sol_pairs_r1, &idx_r1, s);
		if(s==NULL){
			s = mycalloc(1, solution_pair_t);
			s->idx = idx_r1;
			s->r1 = a; s->r2 = b;
			s->prob = a->prob * b->prob;
			HASH_ADD_INT(sol_pairs_r1, idx, s );  /* idx: name of key field */
		}else{ // update with higher score
			if(s->prob < a->prob * b->prob){
				s->prob =  a->prob * b->prob;
				s->r1 = a;
				s->r2 = b;
			}
		}
		// sol_pairs_r2
		HASH_FIND_INT(sol_pairs_r2, &idx_r2, s);
		if(s==NULL){
			s = mycalloc(1, solution_pair_t);
			s->idx = idx_r2;
			s->r1 = c; s->r2 = d;
			s->prob = c->prob * d->prob;
			HASH_ADD_INT(sol_pairs_r2, idx, s);  /* idx: name of key field */
		}else{ // update with higher score
			if(s->prob < c->prob * d->prob){
				s->prob = c->prob * d->prob;
				s->r1 = c; s->r2 = d;
			}
		}
		
	}
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
	printf("%f\t%f\n", likehood1, likehood2);
	/*------------------------------------------------------------------------------*/	
	// find uniquely aligned reads

	/*------------------------------------------------------------------------------*/		
	if(sol_pairs_r1) solution_pair_destory(sol_pairs_r1);
	if(sol_pairs_r2) solution_pair_destory(sol_pairs_r2);
	if(gname1) free(gname1);     if(gname2) free(gname2);
	if(ref1) ref_destory(ref1);  if(ref2) ref_destory(ref2);
}

/* main function. */
int main(int argc, char *argv[]) {
	if((fasta_uthash_load("sample_data/exons.flank_0.fa.gz", &FASTA_HT)) != PR_ERR_NONE) die("main: fasta_uthash_load fails\n");	
	struct BAG_uthash *tb = BAG_uthash_load("graph_flank0.fa");
	struct BAG_uthash *s, *tmp;
	register int i;
	HASH_ITER(hh, tb, s, tmp) {
		edge_align(s, FASTA_HT);
		break;
	}
	BAG_uthash_destroy(&tb);
	fasta_uthash_destroy(&FASTA_HT);
	return 0;
}
