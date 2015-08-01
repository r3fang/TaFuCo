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

static struct kmer_uthash *KMER_HT      = NULL;
static struct fasta_uthash *FASTA_HT    = NULL;
static struct BAG_uthash *BAG_HT        = NULL;

/* main function. */
int main(int argc, char *argv[]) {
	if((fasta_uthash_load("sample_data/exons.flank_0.fa.gz", &FASTA_HT)) != PR_ERR_NONE) die("main: fasta_uthash_load fails\n");	
	struct BAG_uthash *tb = BAG_uthash_load("graph_flank0.fa");
	struct BAG_uthash *s_bag, *tmp_bag;
	register int i;
	solution_pair_t *p;
	junction_t *s_junction, *cur_junction, *tmp_junction;
	HASH_ITER(hh, tb, s_bag, tmp_bag){
		if(((p = edge_align(s_bag, FASTA_HT)))==NULL) return 0;
		if((s_junction = junction_gen(p, s_bag->edge))==NULL) return 0;
		ASH_ITER(hh, s_junction, cur_junction, tmp_junction) {
			printf("name=%s: start=%d\tend=%d\tstr=%s\n", cur_junction->name, cur_junction->start, cur_junction->end, cur_junction->str);
			junction_destory(s_junction);
		}				
	solution_pair_destory(p);
	break;
	}
	BAG_uthash_destroy(&tb);
	fasta_uthash_destroy(&FASTA_HT);
	return 0;
}
