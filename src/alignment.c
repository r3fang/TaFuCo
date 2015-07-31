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


static struct kmer_uthash *KMER_HT      = NULL;
static struct fasta_uthash *FASTA_HT    = NULL;
static struct BAG_uthash *BAG_HT        = NULL;

typedef struct {
	unsigned long idx; // determined by pair(start, end)
	char *name;        // name of edge
	int start;         
	int end;
	char* str;         // reference string
	size_t hits;         // reference string
	double likehood;       // alignment probability
    UT_hash_handle hh;
} junction_t;

///*
// * generate junction sites from solution_pair_t
// */
//junction_t * junction_gen(solution_pair_t *p, char* name){
//	if(p == NULL) return NULL;
//	solution_pair_t *s, *tmp;
//	junction_t *m, *n, *ret = NULL;
//	unsigned int idx;
//    HASH_ITER(hh, p, s, tmp) {
//		// one read
//		if(s->r1->jump == true && s->r1->prob >= MIN_ALIGN_SCORE){
//			idx = pair(s->r1->jump_start, s->r1->jump_end);
//			HASH_FIND_INT(ret, &idx, m);
//			if(m==NULL){ // this junction not in ret
//				m = mycalloc(1, junction_t);
//				m->idx = idx;
//				m->name = strdup(name);
//				m->start = s->r1->jump_start;
//				m->end = s->r1->jump_end;
//				m->hits = 1;
//				m->likehood = 10*log(s->r1->pos); 
//				HASH_ADD_INT(ret, idx, m);
//			}else{
//				m->hits ++;
//				m->likehood += 10*log(s->r1->pos); 
//			}
//		}
//		// the other read
//		if(s->r2->jump == true && s->r2->prob >= MIN_ALIGN_SCORE){
//			idx = pair(s->r2->jump_start, s->r2->jump_end);
//			HASH_FIND_INT(ret, &idx, m);
//			if(m==NULL){ // this junction not in ret
//				m = mycalloc(1, junction_t);
//				m->idx = idx;
//				m->name = strdup(name);
//				m->start = s->r2->jump_start;
//				m->end = s->r2->jump_end;
//				m->hits = 1;
//				m->likehood = 10*log(s->r2->pos); 
//				HASH_ADD_INT(ret, idx, m);
//			}else{
//				m->hits ++;
//				m->likehood += 10*log(s->r2->pos); 
//			}
//		}	
//    }
//	// delete those junctions with hits < MIN_HITS
//	HASH_ITER(hh, ret, m, n){
//		if(m != NULL){
//			m->likehood = m->likehood/m->hits;
//			if(m->hits < MIN_HITS){
//				HASH_DEL(ret,m);
//				free(m);
//			}
//		}
//	}
//	return ret;	
//}
//

/* main function. */
int main(int argc, char *argv[]) {
	if((fasta_uthash_load("sample_data/exons.flank_0.fa.gz", &FASTA_HT)) != PR_ERR_NONE) die("main: fasta_uthash_load fails\n");	
	struct BAG_uthash *tb = BAG_uthash_load("graph_flank0.fa");
	struct BAG_uthash *s, *tmp;
	register int i;

	HASH_ITER(hh, tb, s, tmp) {
			solution_pair_t *t = edge_align(s, FASTA_HT);
			if(t != NULL){
				printf("%f\n", t->prob);
			}
		break;
	}
	BAG_uthash_destroy(&tb);
	fasta_uthash_destroy(&FASTA_HT);
	return 0;
}
