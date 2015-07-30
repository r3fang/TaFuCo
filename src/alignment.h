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

typedef enum { true, false } bool;

#define GAP 					-5.0
#define MATCH 					 2.0
#define MISMATCH 				-1.0
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

typedef struct
{
	char *s;   // string
	size_t l;  // length of the string
	int *S1;    // exon junction sites
	int *S2;    // exon junction sites
	size_t S1_l; // number of exon junction sites
	size_t S2_l; // number of exon junction sites
	size_t J;     // junction site between 2 genes
} ref_t;

typedef struct
{
	char* s1;
	char* s2;
	double score;
} solution;

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
	if(t->S1) free(t->S1);
	if(t->S2) free(t->S2);
	if(t->s)  free(t->s);
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

static inline solution* 
trace_back(matrix_t *S, char *s1, char *s2, int state, int i, int j){
	int flag = 0;
	if(S == NULL || s1 == NULL || s2 == NULL) die("trace_back: paramter error");
	solution *s = mycalloc(1, solution);
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
			case JUMP:
				flag = 1;
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
	s->s1 = res_ks1;
	s->s2 = res_ks2;
	if(flag == 0) s = NULL;
	return s;
}

static inline solution 
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
			delta = ((s1[i-1] - s2[j-1]) == 0) ? MATCH : MISMATCH;
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
	solution *s = trace_back(S, s1, s2, max_state, i_max, j_max);	
	if(s != NULL){
		s->score = max_score;		
	}
	destory_matrix(S);
	return s;
}
#endif