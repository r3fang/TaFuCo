/*--------------------------------------------------------------------*/
/* alingment.h	                                              */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Date: 07-23-2015                                                   */
/* Pair wise fit alignment with affine gap and jump state.            */
/* This could be used to align RNA-seq reads with intron splicing as  */
/* jump state.                                                        */
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
/* G1(i,j) = max{M(i-1, j) + JUMP_EXON (j in S1), G1(i-1, j)}        */
/* G2(i,j) = max{M(i-1, j) + JUMP_EXON (j in S2), G2(i-1, j)}        */
/*  J(i,j) = max{M(i, j-1) + JUMP_GENE (j < x), J(i-1, j)}            */

/* Traceback:                                                         */
/*--------------------------------------------------------------------*/
/* start at max(M(m,j_max), L(n, j_max)), Stop at any of i=0 on M/L;  */
/* The rational behind this is no gap allowed to flank s1             */
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
#include "kseq.h"
#include "kstring.h"

KSEQ_INIT(gzFile, gzread);

int S1[6] = {130, 131, 132, 407, 408, 409};
int S2[6] = {761, 762, 763, 844, 845, 846};
int JUNCTION = 612;

typedef enum { true, false } bool;

#define GAP 					-3.0
#define MATCH 					 2.0
#define MISMATCH 				-0.5
#define EXTENSION               -1.0
#define JUMP_GENE               -15.0
#define JUMP_EXON               -10.0

// POINTER STATE
#define LOW                     500
#define MID                     600
#define UPP                     700
#define JUMP                    800
#define GENE1                   900
#define GENE2                   1000

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

//for alignment allows jump state with junctions
typedef struct {
	size_t size;
	int *pos;
} junction_t;

//opt
typedef struct {
	int o; // gap open
	int e; // gap extension
	int m; // match
	int u; // unmatch
	int j; // jump penality
	bool s;
	junction_t G1;
	junction_t G2;
	int mid;
} opt_t;

static inline opt_t 
*init_opt(){
	opt_t *opt = mycalloc(1, opt_t);
	opt->o = -5.0;
	opt->e = -1.0;
	opt->m =  1.0;
	opt->u = -2.0;
	opt->j = -10.0;
	opt->s = false;
	opt->G1.size = opt->G2.size = 0;	
	opt->G1.pos =  opt->G2.pos = NULL;	
	opt->mid = 0;
	return opt;
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

/*
 * destory kstring
 */
static inline void 
kstring_destory(kstring_t *ks){
	free(ks->s);
	free(ks);
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

static inline void 
trace_back(matrix_t *S, kstring_t *s1, kstring_t *s2, kstring_t *res_ks1, kstring_t *res_ks2, int state, int i, int j){
	if(S == NULL || s1 == NULL || s2 == NULL || res_ks1 == NULL || res_ks2 == NULL) die("trace_back: paramter error");
	int cur = 0; 
	while(i>0){
		switch(state){
			case LOW:
				state = S->pointerL[i][j]; // change to next state
				res_ks1->s[cur] = s1->s[--i];
				res_ks2->s[cur++] = '-';
				break;
			case MID:
				state = S->pointerM[i][j]; // change to next state
                res_ks1->s[cur] = s1->s[--i];
                res_ks2->s[cur++] = s2->s[--j];
				break;
			case UPP:
				state = S->pointerU[i][j];
				res_ks1->s[cur] = '-';
            	res_ks2->s[cur++] = s2->s[--j];
				break;
			case JUMP:
				state = S->pointerJ[i][j];
				res_ks1->s[cur] = '-';
	           	res_ks2->s[cur++] = s2->s[--j];
				break;
			case GENE1:
				state = S->pointerG1[i][j];
				res_ks1->s[cur] = '-';
		        res_ks2->s[cur++] = s2->s[--j];
				break;
			case GENE2:
				state = S->pointerG2[i][j];
				res_ks1->s[cur] = '-';
			    res_ks2->s[cur++] = s2->s[--j];
				break;
			default:
				break;
			}
	}
	res_ks1->l = cur;
	res_ks2->l = cur;	
	res_ks1->s = strrev(res_ks1->s);
	res_ks2->s = strrev(res_ks2->s);
}
double align(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL) die("align: parameter error\n");
	if(s1->l > s2->l) die("first sequence must be shorter than the second to do fitting alignment"); 
	size_t m   = s1->l + 1; size_t n   = s2->l + 1;
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
	for(i=1; i<=s1->l; i++){
		for(j=1; j<=s2->l; j++){
			// MID any state can goto MID
			delta = ((s1->s[i-1] - s2->s[j-1]) == 0) ? MATCH : MISMATCH;
			tmp_J = (j > JUNCTION) ?  S->J[i-1][j-1]+delta : -INFINITY;
			tmp_G1 = (isvalueinarray(j, S1, 6)) ?  S->G1[i-1][j-1]+delta : -INFINITY;
			tmp_G2 = (isvalueinarray(j, S2, 6)) ?  S->G2[i-1][j-1]+delta : -INFINITY;
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
			tmp_M = (isvalueinarray(j, S1, 6)) ?  S->M[i][j-1]+JUMP_EXON : -INFINITY;
			idx = max6(&S->G1[i][j], tmp_M, S->G1[i][j-1], -INFINITY, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerG1[i][j] = MID;			
			if(idx == 1) S->pointerG1[i][j] = GENE1;
			// GENE2
			tmp_M = (isvalueinarray(j, S2, 6)) ?  S->M[i][j-1]+JUMP_EXON : -INFINITY;
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
	i_max = s1->l;
	for(j=0; j<s2->l; j++){
		if(max_score < S->M[i_max][j]){
			max_score = S->M[i_max][j];
			j_max = j;
			max_state = MID;
		}
	}
	for(j=0; j<s2->l; j++){
		if(max_score < S->L[i_max][j]){
			max_score = S->L[i_max][j];
			j_max = j;
			max_state = LOW;
		}
	}
	trace_back(S, s1, s2, r1, r2, max_state, i_max, j_max);	
	destory_matrix(S);
	return max_score;
}
#endif