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

char* exon_concat(struct fasta_uthash *tb, char* gname, int int_arr[MAX_EXON_NUM], size_t *n){
	if(tb==NULL || gname==NULL) die("[%s] input error", __func__);
	struct fasta_uthash *s, *tmp;
	char *ref = NULL;
	int i = 0;
	HASH_ITER(hh, FASTA_HT, s, tmp) {
		if(s!=NULL){
			if(strcmp(strsplit(s->name, '.')[0], gname) == 0){				
				if(ref==NULL){
					ref = strdup(s->seq);
					int_arr[i++] = EXON_FLANK_LEN-2;
					int_arr[i++] = EXON_FLANK_LEN-1;
					int_arr[i++] = EXON_FLANK_LEN;
					int_arr[i++] = EXON_FLANK_LEN+1;
					int_arr[i++] = EXON_FLANK_LEN+2;

					int_arr[i++] = strlen(ref)-EXON_FLANK_LEN-2;
					int_arr[i++] = strlen(ref)-EXON_FLANK_LEN-1;
					int_arr[i++] = strlen(ref)-EXON_FLANK_LEN;
					int_arr[i++] = strlen(ref)-EXON_FLANK_LEN+1;
					int_arr[i++] = strlen(ref)-EXON_FLANK_LEN+2;
				}else{
					int_arr[i++] = strlen(ref)+EXON_FLANK_LEN-2;
					int_arr[i++] = strlen(ref)+EXON_FLANK_LEN-1;
					int_arr[i++] = strlen(ref)+EXON_FLANK_LEN;
					int_arr[i++] = strlen(ref)+EXON_FLANK_LEN+1;
					int_arr[i++] = strlen(ref)+EXON_FLANK_LEN+2;

					char *tmp = concat(ref, s->seq);
					free(ref); ref=strdup(tmp);
					free(tmp);
					int_arr[i++] = strlen(ref)-EXON_FLANK_LEN-2;
					int_arr[i++] = strlen(ref)-EXON_FLANK_LEN-1;
					int_arr[i++] = strlen(ref)-EXON_FLANK_LEN;
					int_arr[i++] = strlen(ref)-EXON_FLANK_LEN+1;
					int_arr[i++] = strlen(ref)-EXON_FLANK_LEN+2;
				}			
			}
		}
	}
	*n = i;
	return ref;
}

void edge_align(struct BAG_uthash *eg, struct fasta_uthash *fasta_u){
	char* gname1 = strsplit(eg->edge, '_')[0];
	char* gname2 = strsplit(eg->edge, '_')[1];
	int gene1_S[MAX_EXON_NUM]; // junction sites for gene1
	int gene2_S[MAX_EXON_NUM]; // junction sites for gene2
	ref_t *ref = ref_init();
	ref->str1 = exon_concat(FASTA_HT, gname1, gene1_S, &ref->n1);
	ref->str2 = exon_concat(FASTA_HT, gname2, gene2_S, &ref->n2);
	ref->S1 = gene1_S;
	ref->S2 = gene2_S;
	
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	
	ks2->s = concat(ref->str1, ref->str2);
	int i; for(i=0; i>ref->n2; i++) ref->S2[i] += strlen(ref->str1); 
	ks2->l = strlen(ks2->s);
	
	i; for(i=0; i<eg->weight; i++){
		ks1->s = strsplit(eg->evidence[i], '_')[0]; ks1->l = strlen(ks1->s);
		if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
		if(ks1->l > ks2->l) die("first sequence must be shorter than the second\n");
		kstring_t *r1 = mycalloc(1, kstring_t);
		kstring_t *r2 = mycalloc(1, kstring_t);
		r1->s = mycalloc(ks1->l + ks2->l, char);
		r2->s = mycalloc(ks1->l + ks2->l, char);
		printf("score=%f\n", align(ks1, ks2, r1, r2, ref));
		printf("%s\n%s\n", r1->s, r2->s);	
		//kstring_destory(ks1);
	}
		
	kstring_destory(ks2);
	//kstring_destory(r1);
	//kstring_destory(r2);

	if(gname1) free(gname1);
	if(gname2) free(gname2);
	if(ref) free(ref);
}
/* main function. */
int main(int argc, char *argv[]) {
	struct BAG_uthash *tb = mycalloc(1, struct BAG_uthash);
	tb->edge = "EGFR_PLAG1";
	tb->weight = 3;
	tb->evidence = mycalloc(3, char*);
	tb->evidence[0] = "GGAAAAAAATTTTAGCCTTATCTTAATCTGTCCCAACAGCAATGTGACGGATTTTTGCAGATTCAAAATCTGCAATGGTTATTTACAAGTCAATCT";
	tb->evidence[1] = "GATCATTCTACAAGATGTCAGTGCACTGAAACATGCAGGGGCGTGTTGAGTGTGGAAGGATCTTGACAAGTTGTTT_CCTGCCCGGCAGGAGTCATGGGAGAAAACAACACCCTGGTCTGGAAGTACGCAGACGCCGGCCATGTGTGCCACCT";
	tb->evidence[2] = "ATTCTACAAGATGTCAGTGCACTGAAACATGCAGGGGCGTGTTGAGTGTGGAAGGATCTTGACAAGTTGTTTTGAA_ATCCAAACTGCACCTACGGATGCACTGGGCCAGGTCTTGAAGGCTGTCCAACGAATGGAAGCTACATAGTGTCTC";
	if((fasta_uthash_load("sample_data/exons.fa", &FASTA_HT)) != PR_ERR_NONE) die("main: fasta_uthash_load fails\n");	
	edge_align(tb, FASTA_HT);
	
	return 0;
}
