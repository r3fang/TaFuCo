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

#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include <stdio.h>  
#include <stdlib.h>  
#include <string.h> 
#include <zlib.h>  
#include <assert.h>
#include <math.h>
#include <regex.h>
#include "utils.h"

// alignment state
#define LOW                     500
#define MID                     600
#define UPP                     700
#define JUMP                    800
#define GENE1                   900
#define GENE2                   1000

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

// alingment of a single read
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

// alingment of a read pair
typedef struct
{
	char* idx;
	solution_t *r1;
	solution_t *r2;
	double prob;
	UT_hash_handle hh;
} solution_pair_t;

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
	S->G1= mycalloc(m, double*);
	S->G2= mycalloc(m, double*);
	
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

static inline char* 
idx2str(char* name, int i, int j){
	if(name == NULL) die("[%s] input error", __func__);
	// convert alignment information to md5
	char istr[1000];
	char jstr[1000];
	sprintf(istr, "%d", i);
	sprintf(jstr, "%d", j);
	return concat(concat(concat(concat(name, "."), istr), "."), jstr); 
}

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
	if(s->idx) free(s->idx);
	if(s->r1) solution_destory(s->r1);
	if(s->r2) solution_destory(s->r2);
	free(s);
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
			default:
				break;
			}
	}
	s->s1 = strrev(res_ks1);
	s->s2 = strrev(res_ks2);
	s->pos = j;
	return s;
}

static inline solution_t 
*align(char *s1, char *s2, int junction, opt_t* opt){
	if(s1 == NULL || s2 == NULL || opt==NULL) die("[%s] parameter error", __func__);
	double MATCH = opt->match;
	double MISMATCH = opt->mismatch;
	double GAP = opt->gap;
	double EXTENSION = opt->extension;
	double JUMP_GENE = opt->jump_gene;
	
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
	}
	// initlize first row
	for(j=0; j<S->n; j++){
		S->M[0][j] = 0.0;
		S->U[0][j] = 0.0;
		S->L[0][j] = -INFINITY;
		S->J[0][j] = -INFINITY;
	}
	double delta, tmp_J, tmp_M;
	int idx;
	
	// recurrance relation
	for(i=1; i<=strlen(s1); i++){
		for(j=1; j<=strlen(s2); j++){
			// MID any state can goto MID
			delta = ((toupper(s1[i-1]) - toupper(s2[j-1])) == 0) ? MATCH : MISMATCH;
			tmp_J = (j > junction) ?  S->J[i-1][j-1]+delta : -INFINITY;
			idx = max6(&S->M[i][j], S->L[i-1][j-1]+delta, S->M[i-1][j-1]+delta, S->U[i-1][j-1]+delta, tmp_J, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerM[i][j]=LOW;
			if(idx == 1) S->pointerM[i][j]=MID;
			if(idx == 2) S->pointerM[i][j]=UPP;
			if(idx == 3) S->pointerM[i][j]=JUMP;			 				
			// LOW
			idx = max6(&S->L[i][j], S->L[i-1][j]+EXTENSION, S->M[i-1][j]+GAP, -INFINITY, -INFINITY, -INFINITY,  -INFINITY);
			if(idx == 0) S->pointerL[i][j]=LOW;
			if(idx == 1) S->pointerL[i][j]=MID;
			// UPP
			idx = max6(&S->U[i][j], S->M[i][j-1]+GAP, S->U[i][j-1]+EXTENSION,  -INFINITY, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerU[i][j]=MID;
			if(idx == 1) S->pointerU[i][j]=UPP;
			// JUMP 
			tmp_M = (j < junction) ?  S->M[i][j-1]+JUMP_GENE : -INFINITY;
			idx = max6(&S->J[i][j], tmp_M, S->J[i][j-1], -INFINITY,  -INFINITY,  -INFINITY, -INFINITY);
			if(idx == 0) S->pointerJ[i][j] = MID;			
			if(idx == 1) S->pointerJ[i][j] = JUMP;			
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

static inline solution_t
*trace_back_exon_jump(matrix_t *S, char *s1, char *s2, int state, int i, int j){
	if(S == NULL || s1 == NULL || s2 == NULL) return NULL;
	solution_t *s = solution_init();
	char *res_ks1 = mycalloc(strlen(s1)+strlen(s2), char); 
	char *res_ks2 = mycalloc(strlen(s1)+strlen(s2), char); 
	int cur = 0; 
	while(i>0){
		switch(state){
			case LOW:
				state = S->pointerL[i][j]; // change to next state
				res_ks1[cur] = s1[--i];
				res_ks2[cur++] = '-';
				break;
			case MID:
				state = S->pointerM[i][j]; // change to next state
                res_ks1[cur] = s1[--i];
                res_ks2[cur++] = s2[--j];
				break;
			case UPP:
				state = S->pointerU[i][j];
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
	s->s1 = strrev(res_ks1);
	s->s2 = strrev(res_ks2);
	s->pos = j;
	return s;
}

static inline solution_t *align_exon_jump(char *s1, char *s2, int *S1, int *S2, int S1_num, int S2_num, opt_t *opt){
	if(s1 == NULL || s2 == NULL || opt==NULL) return NULL;
	if(strlen(s1) > strlen(s2)) return NULL; 
	
	size_t m   = strlen(s1) + 1; 
	size_t n   = strlen(s2) + 1;
	matrix_t *S = create_matrix(m, n);
	// initlize leftmost column
	int i, j;
	for(i=0; i<S->m; i++){
		S->M[i][0] = -INFINITY;
		S->U[i][0] = -INFINITY;
		S->L[i][0] = -INFINITY;
		S->G1[i][0] = -INFINITY;
		S->G2[i][0] = -INFINITY;
	}
	// initlize first row
	for(j=0; j<S->n; j++){
		S->M[0][j] = 0.0;
		S->U[0][j] = 0.0;
		S->L[0][j] = -INFINITY;
		S->G1[0][j]  = -INFINITY;
		S->G2[0][j]  = -INFINITY;
	}
	double delta, tmp_M, tmp_G1, tmp_G2; 
	int idx;
	// recurrance relation
	for(i=1; i<=strlen(s1); i++){
		for(j=1; j<=strlen(s2); j++){
			// MID any state can goto MID
			delta = ((toupper(s1[i-1]) - toupper(s2[j-1])) == 0) ? opt->match : opt->mismatch;
			tmp_G1 = -INFINITY;
			tmp_G2 = -INFINITY;
			//tmp_G1 = (isvalueinarray(j, S1, S1_num)) ?  S->G1[i-1][j-1]+delta : -INFINITY;
			//tmp_G2 = (isvalueinarray(j, S2, S2_num)) ?  S->G2[i-1][j-1]+delta : -INFINITY;
			idx = max6(&S->M[i][j], S->L[i-1][j-1]+delta, S->M[i-1][j-1]+delta, S->U[i-1][j-1]+delta, tmp_G1, tmp_G2, -INFINITY);
			if(idx == 0) S->pointerM[i][j]=LOW;
			if(idx == 1) S->pointerM[i][j]=MID;
			if(idx == 2) S->pointerM[i][j]=UPP;
			if(idx == 3) S->pointerM[i][j]=GENE1;
			if(idx == 4) S->pointerM[i][j]=GENE2;
			// LOW
			idx = max6(&S->L[i][j], S->L[i-1][j]+opt->extension, S->M[i-1][j]+opt->gap, -INFINITY, -INFINITY, -INFINITY,  -INFINITY);
			if(idx == 0) S->pointerL[i][j]=LOW;
			if(idx == 1) S->pointerL[i][j]=MID;
			// UPP
			idx = max6(&S->U[i][j], S->M[i][j-1]+opt->gap, S->U[i][j-1]+opt->extension,  -INFINITY, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerU[i][j]=MID;
			if(idx == 1) S->pointerU[i][j]=UPP;
			// G1
			//tmp_M = (isvalueinarray(j, S1, S1_num)) ?  S->M[i][j-1]+opt->jump_exon : -INFINITY;
			tmp_M = -INFINITY;
			idx = max6(&S->G1[i][j], tmp_M, S->G1[i][j-1], -INFINITY, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerG1[i][j] = MID;			
			if(idx == 1) S->pointerG1[i][j] = GENE1;
			//G2
			//tmp_M = (isvalueinarray(j, S2, S2_num)) ?  S->M[i][j-1]+opt->jump_exon : -INFINITY;
			tmp_M = -INFINITY;
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
	solution_t *s = trace_back_exon_jump(S, s1, s2, max_state, i_max, j_max);	
	destory_matrix(S);
	s->score = max_score;	
	//s->prob = max_score/(MATCH*strlen(s1));	
	return s;
}


#endif