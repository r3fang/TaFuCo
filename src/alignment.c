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
#define EXON_FLANK_LEN                    50

static struct kmer_uthash *KMER_HT      = NULL;
static struct fasta_uthash *FASTA_HT    = NULL;
static struct BAG_uthash *BAG_HT        = NULL;

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
	int gene1_S[MAX_EXON_NUM]; // junction sites for gene1
	int gene2_S[MAX_EXON_NUM]; // junction sites for gene2
	ref_t *ref1 = ref_generate(FASTA_HT, gname1, gname2);
	ref_t *ref2 = ref_generate(FASTA_HT, gname2, gname1);
	solution **sol1_E1 = mycalloc(eg->weight, solution*);
	solution **sol1_E2 = mycalloc(eg->weight, solution*);
	solution **sol2_E1 = mycalloc(eg->weight, solution*);
	solution **sol2_E2 = mycalloc(eg->weight, solution*);
	
	double score1 = 0;
	double score2 = 0;
	
	register int i; for(i=0; i<eg->weight; i++){
		sol1_E1[i] = align(strsplit(eg->evidence[i], '_')[0], ref1);
		sol1_E2[i] = align(strsplit(eg->evidence[i], '_')[1], ref1);
		sol2_E1[i] = align(strsplit(eg->evidence[i], '_')[0], ref2);
		sol2_E2[i] = align(strsplit(eg->evidence[i], '_')[1], ref2);
		score1 += log(sol1_E1[i]->score);
		score1 += log(sol1_E2[i]->score);
		score2 += log(sol2_E1[i]->score);
		score2 += log(sol2_E2[i]->score);
		
		//if(a->jump || b->jump){
		//	printf("ref1=%f\tref2=%f\n", a->score, b->score);			
			//}
		//if(a->score > b->score && a->jump == true){
		//	printf("score=%f\tpos=%d\tjunction=%d\n", a->score, a->pos, a->jump_pos);
		//	printf("%s\n", a->s1);
		//	printf("%s\n", a->s2);
		//}
		//if(b->score >= a->score && b->jump == true){
		//if(b->jump == true || a->jump == true)
			//printf("score_a=%f\tpos_a=%d\tjump_start_a=%d\tjump_end_a=%d\tscore_b=%f\tpos_b=%d\tjump_start_b=%d\tjump_end_b=%d\n", a->score, a->pos, a->jump_start, a->jump_end, b->score, b->pos, b->jump_start, b->jump_end);
		//	printf("score_a=%f\tpos_a=%d\tinsertion=%d\tdeletion=%d\tjump_start_a=%d\tjump_end_a=%d\n", a->score, a->pos, a->insertion, a->deletion, a->jump_start, a->jump_end);
			//printf("score_a=%f\tpos_a=%d\tinsertion=%d\tdeletion=%d\tjump_start_a=%d\tjump_end_a=%d\n", b->score, b->pos, b->insertion, b->deletion, b->jump_start, b->jump_end);

			//printf("%s\n", b->s1);
			//printf("%s\n", b->s2);
		//}
	}
	printf("score1=%f\tscore1=%f\n", score1, score2);
	for(i=0; i<eg->weight; i++) free(sol1_E1[i]); free(sol1_E1);
	for(i=0; i<eg->weight; i++) free(sol1_E2[i]); free(sol1_E2);
	for(i=0; i<eg->weight; i++) free(sol2_E1[i]); free(sol2_E1);
	for(i=0; i<eg->weight; i++) free(sol2_E2[i]); free(sol2_E2);
		
	if(gname1) free(gname1);
	if(gname2) free(gname2);
	if(ref1) ref_destory(ref1);
	if(ref2) ref_destory(ref2);
}

/* main function. */
int main(int argc, char *argv[]) {
	if((fasta_uthash_load("sample_data/exons.flank_0.fa.gz", &FASTA_HT)) != PR_ERR_NONE) die("main: fasta_uthash_load fails\n");	
	struct BAG_uthash *tb = BAG_uthash_load("graph_flank0.fa");
	struct BAG_uthash *s, *tmp;
	register int i;
	HASH_ITER(hh, tb, s, tmp) {
		edge_align(s, FASTA_HT);
	}
	BAG_uthash_destroy(&tb);
	fasta_uthash_destroy(&FASTA_HT);
	return 0;
}
