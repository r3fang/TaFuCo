/*--------------------------------------------------------------------*/
/* alingment.h                                                        */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Created Date: 07-23-2015                                           */
/* Pair wise fit alignment with affine gap and jump state.            */
/* This could be used to align RNA-seq reads with intron splicing     */
/* and fusion between two candidate genes.                            */
/*                                                                    */
/* initilize L(i,j), M(i,j), U(i,j), J(i,j), G1(i,j), G2(i,j):        */
/*--------------------------------------------------------------------*/
/* M(0,0)  = -INF                                                     */
/* M(i,0)  = -INF               (i>0 )                                */
/* M(0,j)  = 0                  (j>0 )                                */
/* L(0,j)  = -INF               (j>=0)                                */
/* L(i,0)  = -INF               (i>=0)                                */
/* U(i,0)  = -INF               (i>0 )                                */
/* U(0,j)  = 0                  (j>=0)                                */
/* G1(0,j) = -INF               (j>=0)                                */
/* G1(i,0) = -INF               (i>=0)                                */
/* G2(0,j) = -INF               (j>=0)                                */
/* G2(i,0) = -INF               (i>=0)                                */
/* J(0,j)  = -INF               (j>=0)                                */
/* J(i,0)  = -INF               (i>=0)                                */
/*                                                                    */
/* reccurrance relations:                                             */
/*--------------------------------------------------------------------*/
/*  M(i,j) = max{M(i-1, j-1) + s(x,y), U(i-1, j-1) + s(x,y),          */
/*              L(i-1, j-1) + s(x,y),                                 */
/*              J(i-1, j-1) + s(x,y) (j>x),                           */
/*              G1(i-1, j-1) + s(x,y) if j in S1,                     */
/*              G2(i-1, j-1) + s(x,y) if j in S2}                     */                                          
/*  U(i,j) = max{M(i-1, j) + GAP, U(i-1, j) + EXTENSION}              */
/*  L(i,j) = max{M(i, j-1) + GAP, L(i, j-1) + EXTENSION}              */
/*--------------------------------------------------------------------*/
/* G1(i,j) = max{M(i-1, j) + JUMP_EXON (j in S1), G1(i-1, j)}         */
/* G2(i,j) = max{M(i-1, j) + JUMP_EXON (j in S2), G2(i-1, j)}         */
/*  J(i,j) = max{M(i, j-1) + JUMP_GENE (j < x), J(i-1, j)}            */
/* Traceback:                                                         */
/*--------------------------------------------------------------------*/
/* Start at max(M(m, j_max), L(n, j_max)),                            */
/* Stop at topmost row of M or L                                      */
/* The rational behind this is alignment should not start with gap    */
/* read.                                                              */
/*--------------------------------------------------------------------*/
/* NOTE:                                                              */
/* 1. WE ONLY ALLOW STATE CHANGE FROM M(MATCH) TO J(JUMP) AT SPECIFIC */
/* POSITIONS ON S2.                                                   */ 
/* 2. WE ALLOW STATE CHANGE FROM JUMP TO MATCH AT ANYWHERE            */
/* 3. THE READ SHOULD NOT START WITH JUMP AND ALSO SHOULD NOT END     */
/* WITH JUMP STATE.                                                   */ 
/*--------------------------------------------------------------------*/

#ifndef _ALIGNMENT_
#define _ALIGNMENT_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <math.h>
#include "utils.h"

// constant define
typedef enum { true, false } bool;
// input for alignment
#define GAP 					-5.0
#define MATCH 					 2.0
#define MISMATCH 				-1.0
#define EXTENSION               -2.0
#define JUMP_EXON               -10.0
#define JUMP_GENE               -15.0

// input for junction identification
#define MIN_ALIGN_SCORE          0.7
#define MIN_HITS                 3
#define EXON_FLANK_LEN           0

// alignment state
#define LOW                     500
#define MID                     600
#define UPP                     700
#define JUMP                    800
#define GENE1                   900
#define GENE2                   1000

#define pair(k1, k2)  ((k1 + k2)*(k1 + k2 + 1)/2 + k2)

// dynamic programming matrices
typedef struct {
  unsigned int m;
  unsigned int n;
  double **L;
  double **M;
  double **U;
  double **J;
  double **G1;
  double **G2;
  int  **pointerL;
  int  **pointerM;
  int  **pointerU;
  int  **pointerJ;
  int  **pointerG1;
  int  **pointerG2;
} matrix_t;

/*
 * create matrix, allocate memor
 */
static inline matrix_t 
*create_matrix(size_t m, size_t n){
	size_t i, j; 
	matrix_t *S = mycalloc(1, matrix_t);
	S->m = m;
	S->n = n;
	S->L = mycalloc(m, double*);
	S->M = mycalloc(m, double*);
	S->U = mycalloc(m, double*);
	S->J = mycalloc(m, double*);
	S->G1 = mycalloc(m, double*);
	S->G2 = mycalloc(m, double*);
	
	for (i = 0; i < m; i++) {
      S->M[i] = mycalloc(n, double);
      S->L[i] = mycalloc(n, double);
      S->U[i] = mycalloc(n, double);
      S->J[i] = mycalloc(n, double);
      S->G1[i] = mycalloc(n, double);
      S->G2[i] = mycalloc(n, double);
	  
    }	
	
	S->pointerM = mycalloc(m, int*);
	S->pointerU = mycalloc(m, int*);
	S->pointerL = mycalloc(m, int*);
	S->pointerJ = mycalloc(m, int*);
	S->pointerG1 = mycalloc(m, int*);
	S->pointerG2 = mycalloc(m, int*);
	
	for (i = 0; i < m; i++) {
       	S->pointerU[i] = mycalloc(n, int);
        S->pointerM[i] = mycalloc(n, int);
        S->pointerL[i] = mycalloc(n, int);
		S->pointerJ[i] = mycalloc(n, int);
		S->pointerG1[i] = mycalloc(n, int);
		S->pointerG2[i] = mycalloc(n, int);		
    }
	return S;
}

/*
 * destory matrix
 */
static inline void 
destory_matrix(matrix_t *S){
	if(S == NULL) die("destory_matrix: parameter error\n");
	int i;
	for(i = 0; i < S->m; i++){
		if(S->L[i]) free(S->L[i]);
		if(S->M[i]) free(S->M[i]);
		if(S->U[i]) free(S->U[i]);
		if(S->J[i]) free(S->J[i]);
		if(S->G1[i]) free(S->G1[i]);
		if(S->G2[i]) free(S->G2[i]);
	}
	for(i = 0; i < S->m; i++){
		if(S->pointerL[i]) free(S->pointerL[i]);
		if(S->pointerM[i]) free(S->pointerM[i]);
		if(S->pointerU[i]) free(S->pointerU[i]);
		if(S->pointerJ[i]) free(S->pointerJ[i]);
		if(S->pointerG1[i]) free(S->pointerG1[i]);
		if(S->pointerG2[i]) free(S->pointerG2[i]);	
	}
	free(S);
}

// junction of gene fusion
typedef struct {
	unsigned long idx; // determined by pair(start, end)
	char *name;        // name of edge
	int start;         
	int end;
	char* str;         // string flanking junction site 
	size_t hits;         // reference string
	double likehood;       // alignment probability
    UT_hash_handle hh;
} junction_t;

// initlize junction_t
static inline junction_t 
*junction_init(){
	junction_t *j = mycalloc(1, junction_t);
	j->name = NULL;
	j->str = NULL;
	j->likehood = 0;
	j->hits = 0;
	return j;
}
// destory junction
static inline void 
junction_destory(junction_t *s){
	if(s==NULL) die("[%s] input error", __func__);
	if(s->name) free(s->name);
	if(s->str) free(s->str);
	free(s);
}
/*
 * concatenated string of exons.
 */ 
typedef struct
{
	char *s;      // string
	size_t l;     // length of the string
	int *S1;      // exon junction sites
	int *S2;      // exon junction sites
	size_t S1_l;  // number of exon junction sites
	size_t S2_l;  // number of exon junction sites
	size_t J;     // junction site between 2 genes
} ref_t;

// initilize ref_t
static inline ref_t 
*ref_init(){
	ref_t *r = mycalloc(1, ref_t);
	r->S1_l  = 0;
	r->S2_l  = 0;
	r->s = NULL;
	r->l = 0;
	r->J = 0;
	return r;
}

static inline void 
ref_destory(ref_t* t){
	if(t==NULL) die("[%s] input error", __func__);
	if(t->S1) free(t->S1);
	if(t->S2) free(t->S2);
	if(t->s)  free(t->s);
	free(t);
}

// alingment soulution of a single read
typedef struct
{
	char* s1; 
	char* s2;
	double score;
	bool jump;
	int jump_start;
	int jump_end;
	int pos;
	int match;
	int insertion;
	int deletion;
	double prob;
} solution_t;

// initilize solution_t
static inline solution_t 
*solution_init(){
	solution_t *t = mycalloc(1, solution_t);
	t->s1 = NULL;
	t->s2 = NULL;
	t->score = 0;
	t->jump = false;
	t->jump_start = 0;
	t->jump_end = 0;
	t->pos = 0;
	t->match = 0;
	t->insertion = 0;
	t->deletion = 0;
	t->prob = 0.0;
	return t;
}
// destory solution_t
static inline void solution_destory(solution_t *s){
	if(s->s1) free(s->s1);
	if(s->s2) free(s->s2);
	free(s);
}

// alingment soulution for a read pair
typedef struct
{
	unsigned long idx;
	solution_t *r1;
	solution_t *r2;
	double prob;
	UT_hash_handle hh;
} solution_pair_t;

static inline solution_pair_t 
*solution_pair_init(){
	solution_pair_t* s = mycalloc(1, solution_pair_t);
	s->r1 = solution_init();
	s->r2 = solution_init();
	return s;
}

static inline void 
solution_pair_destory(solution_pair_t *s){
	if(s==NULL) die("[%s] input error", __func__);	
	if(s->r1) solution_destory(s->r1);
	if(s->r2) solution_destory(s->r2);
	free(s);
}

/* max of fix values */
static inline int 
max6(double *res, double a1, double a2, double a3, double a4, double a5, double a6){
	*res = -INFINITY;
	int state;
	if(a1 > *res){*res = a1; state = 0;}
	if(a2 > *res){*res = a2; state = 1;}
	if(a3 > *res){*res = a3; state = 2;}	
	if(a4 > *res){*res = a4; state = 3;}	
	if(a5 > *res){*res = a5; state = 4;}	
	if(a6 > *res){*res = a6; state = 5;}	
	return state;
}

static inline char 
*strrev(char *s){
	if(s == NULL) return NULL;
	int l = strlen(s);
	char *ss = strdup(s);
	free(s);
	s = mycalloc(l, char);
	int i; for(i=0; i<l; i++){
		s[i] = ss[l-i-1];
	}
	s[l] = '\0';
	return s;
}


static inline bool 
isvalueinarray(int val, int *arr, int size){
    int i;
    for (i=0; i < size; i++) {
        if (arr[i] == val)
            return TRUE;
    }
    return FALSE;
}

static inline solution_t* 
trace_back(matrix_t *S, char *s1, char *s2, int state, int i, int j){
	if(S == NULL || s1 == NULL || s2 == NULL) die("trace_back: paramter error");
	solution_t *s = solution_init();
	s->jump = false;
	char *res_ks1 = mycalloc(strlen(s1)+strlen(s2), char); 
	char *res_ks2 = mycalloc(strlen(s1)+strlen(s2), char); 
	int cur = 0; 
	s->jump_start = s->jump_end = 0;
	while(i>0){
		switch(state){
			case LOW:
				s->deletion ++;
				state = S->pointerL[i][j]; // change to next state
				res_ks1[cur] = s1[--i];
				res_ks2[cur++] = '-';
				break;
			case MID:
				s->match ++;
				state = S->pointerM[i][j]; // change to next state
                res_ks1[cur] = s1[--i];
                res_ks2[cur++] = s2[--j];
				break;
			case UPP:
				s->insertion ++;
				state = S->pointerU[i][j];
				res_ks1[cur] = '-';
            	res_ks2[cur++] = s2[--j];
				break;
			case JUMP:
				s->jump = true;
				if(j > s->jump_end) s->jump_end = j;
				s->jump_start=j;
				state = S->pointerJ[i][j];
				res_ks1[cur] = '-';
	           	res_ks2[cur++] = s2[--j];
				break;
			case GENE1:
				state = S->pointerG1[i][j];
				res_ks1[cur] = '-';
		        res_ks2[cur++] = s2[--j];
				break;
			case GENE2:
				state = S->pointerG2[i][j];
				res_ks1[cur] = '-';
			    res_ks2[cur++] = s2[--j];
				break;
			default:
				break;
			}
	}
	s->s1 = s1;
	s->s2 = s2;
	s->pos = j;
	return s;
}

static inline solution_t 
*align(char *s1, ref_t *ref){
	if(s1 == NULL || ref == NULL) die("align: parameter error\n");
	char *s2 = ref->s;
	int  *S1 = ref->S1;
	int  *S2 = ref->S2;
	int  JUNCTION = ref->J;
	if(strlen(s1) > strlen(s2)) die("first sequence must be shorter than the second to do fitting alignment"); 
	size_t m   = strlen(s1) + 1; 
	size_t n   = strlen(s2) + 1;
	matrix_t *S = create_matrix(m, n);
	// initlize leftmost column
	int i, j;
	for(i=0; i<S->m; i++){
		S->M[i][0] = -INFINITY;
		S->U[i][0] = -INFINITY;
		S->L[i][0] = -INFINITY;
		S->J[i][0] = -INFINITY;
		S->G1[i][0] = -INFINITY;
		S->G2[i][0] = -INFINITY;
	}
	// initlize first row
	for(j=0; j<S->n; j++){
		S->M[0][j] = 0.0;
		S->U[0][j] = 0.0;
		S->L[0][j] = -INFINITY;
		S->J[0][j] = -INFINITY;
		S->G1[0][j] = -INFINITY;
		S->G2[0][j] = -INFINITY;
	}
	double delta, tmp_J, tmp_G1, tmp_G2, tmp_M;
	int idx;
	
	// recurrance relation
	for(i=1; i<=strlen(s1); i++){
		for(j=1; j<=strlen(s2); j++){
			// MID any state can goto MID
			delta = ((toupper(s1[i-1]) - toupper(s2[j-1])) == 0) ? MATCH : MISMATCH;
			tmp_J = (j > JUNCTION) ?  S->J[i-1][j-1]+delta : -INFINITY;
			tmp_G1 = (isvalueinarray(j, S1, ref->S1_l)) ?  S->G1[i-1][j-1]+delta : -INFINITY;
			tmp_G2 = (isvalueinarray(j, S2, ref->S2_l)) ?  S->G2[i-1][j-1]+delta : -INFINITY;
			idx = max6(&S->M[i][j], S->L[i-1][j-1]+delta, S->M[i-1][j-1]+delta, S->U[i-1][j-1]+delta, tmp_J, tmp_G1, tmp_G2);
			if(idx == 0) S->pointerM[i][j]=LOW;
			if(idx == 1) S->pointerM[i][j]=MID;
			if(idx == 2) S->pointerM[i][j]=UPP;
			if(idx == 3) S->pointerM[i][j]=JUMP;			 				
			if(idx == 4) S->pointerM[i][j]=GENE1;			 				
			if(idx == 5) S->pointerM[i][j]=GENE2;			 				
			// LOW
			idx = max6(&S->L[i][j], S->L[i-1][j]+EXTENSION, S->M[i-1][j]+GAP, -INFINITY, -INFINITY, -INFINITY,  -INFINITY);
			if(idx == 0) S->pointerL[i][j]=LOW;
			if(idx == 1) S->pointerL[i][j]=MID;
			// UPP
			idx = max6(&S->U[i][j], S->M[i][j-1]+GAP, S->U[i][j-1]+EXTENSION,  -INFINITY, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerU[i][j]=MID;
			if(idx == 1) S->pointerU[i][j]=UPP;
			// JUMP 
			tmp_M = (j < JUNCTION) ?  S->M[i][j-1]+JUMP_GENE : -INFINITY;
			idx = max6(&S->J[i][j], tmp_M, S->J[i][j-1], -INFINITY,  -INFINITY,  -INFINITY, -INFINITY);
			if(idx == 0) S->pointerJ[i][j] = MID;			
			if(idx == 1) S->pointerJ[i][j] = JUMP;			
			// GENE1
			tmp_M = (isvalueinarray(j, S1, ref->S1_l)) ?  S->M[i][j-1]+JUMP_EXON : -INFINITY;
			idx = max6(&S->G1[i][j], tmp_M, S->G1[i][j-1], -INFINITY, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerG1[i][j] = MID;			
			if(idx == 1) S->pointerG1[i][j] = GENE1;
			// GENE2
			tmp_M = (isvalueinarray(j, S2, ref->S2_l)) ?  S->M[i][j-1]+JUMP_EXON : -INFINITY;
			idx = max6(&S->G2[i][j], tmp_M, S->G2[i][j-1], -INFINITY, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerG2[i][j] = MID;			
			if(idx == 1) S->pointerG2[i][j] = GENE2;
			}
		}
	// find trace-back start point
	// NOTE: ALWAYS STARTS TRACING BACK FROM MID OR LOW
	int i_max, j_max;
	double max_score = -INFINITY;
	int max_state;
	i_max = strlen(s1);
	for(j=0; j<strlen(s2); j++){
		if(max_score < S->M[i_max][j]){
			max_score = S->M[i_max][j];
			j_max = j;
			max_state = MID;
		}
	}
	for(j=0; j<strlen(s2); j++){
		if(max_score < S->L[i_max][j]){
			max_score = S->L[i_max][j];
			j_max = j;
			max_state = LOW;
		}
	}
	solution_t *s = trace_back(S, s1, s2, max_state, i_max, j_max);	
	s->score = max_score;		
	destory_matrix(S);
	if(s->jump == false) {s->prob = s->score/(strlen(s1)*MATCH);}
	if(s->jump == true)  {s->prob = s->score/(strlen(s1)*MATCH+JUMP_GENE);}
	return s;
}

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
 *
 * Align reads that support e(Vi, Vj) to the string concated by exons of 
 * Vi and Vj.
 *
 */
static inline solution_pair_t* 
edge_align(struct BAG_uthash *eg, struct fasta_uthash *fasta_u){
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
	register int i, j;
	int idx_r1, idx_r2;

	for(i=0; i<eg->weight; i++){
		solution_t *a = align(strsplit(eg->evidence[i], '_')[0], ref1);
		solution_t *b = align(strsplit(eg->evidence[i], '_')[1], ref1);
		solution_t *c = align(strsplit(eg->evidence[i], '_')[0], ref2);
		solution_t *d = align(strsplit(eg->evidence[i], '_')[1], ref2);
		idx_r1 = pair(a->pos, b->pos); idx_r2 = pair(c->pos, d->pos);
		// sol_pairs_r1
		HASH_FIND_INT(sol_pairs_r1, &idx_r1, s);
		if(s==NULL){
			s = solution_pair_init();
			s->idx = idx_r1;
			s->r1 = a; s->r2 = b;
			s->prob = a->prob * b->prob;
			HASH_ADD_INT(sol_pairs_r1, idx, s );  /* idx: name of key field */
		}else{ // update with higher score
			if(s->prob < a->prob * b->prob){
				s->prob =  a->prob * b->prob;
				s->r1 = a; s->r2 = b;
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
	if(ref1) ref_destory(ref1);  if(ref2) ref_destory(ref2);
	printf("likehood1=%f\tlikehood1=%f\n", likehood1, likehood2);
	solution_pair_t * ret;
	if(likehood1 >= likehood2){
		//if(sol_pairs_r2) solution_pair_destory(sol_pairs_r2);
		ret = sol_pairs_r1;
	}else{
		//if(sol_pairs_r1) solution_pair_destory(sol_pairs_r1);
		ret = sol_pairs_r2;
	}
	
	if(gname1) free(gname1);     
	if(gname2) free(gname2);
	return ret;
}
/*
 * generate junction sites from solution_pair_t
 */
static inline junction_t 
*junction_gen(solution_pair_t *p, char* name){
	if(p == NULL) return NULL;
	solution_pair_t *s, *tmp;
	junction_t *m, *n, *ret = NULL;
	unsigned int idx;
    HASH_ITER(hh, p, s, tmp) {
		// one read
		if(s->r1->jump == true && s->r1->prob >= MIN_ALIGN_SCORE){
			idx = pair(s->r1->jump_start, s->r1->jump_end);
			HASH_FIND_INT(ret, &idx, m);
			if(m==NULL){ // this junction not in ret
				m = junction_init();
				m->idx = idx;
				m->name = strdup(name);
				m->start = s->r1->jump_start;
				m->end = s->r1->jump_end;
				m->hits = 1;
				m->likehood = 10*log(s->r1->pos); 
				m->str = strdup(s->r1->s2);
				HASH_ADD_INT(ret, idx, m);
			}else{
				m->hits ++;
				m->likehood += 10*log(s->r1->pos); 
			}
		}
		// the other read
		if(s->r2->jump == true && s->r2->prob >= MIN_ALIGN_SCORE){
			idx = pair(s->r2->jump_start, s->r2->jump_end);
			HASH_FIND_INT(ret, &idx, m);
			if(m==NULL){ // this junction not in ret
				m = mycalloc(1, junction_t);
				m->idx = idx;
				m->name = strdup(name);
				m->start = s->r2->jump_start;
				m->end = s->r2->jump_end;
				m->hits = 1;
				m->likehood = 10*log(s->r2->pos); 
				m->str = strdup(s->r2->s2);
				HASH_ADD_INT(ret, idx, m);
			}else{
				m->hits ++;
				m->likehood += 10*log(s->r2->pos); 
			}
		}	
    }
	// delete those junctions with hits < MIN_HITS
	HASH_ITER(hh, ret, m, n){
		if(m != NULL){
			m->likehood = m->likehood/m->hits;
			if(m->hits < MIN_HITS){
				HASH_DEL(ret,m);
				free(m);
			}
		}
	}
	return ret;	
}

#endif