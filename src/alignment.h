/*--------------------------------------------------------------------*/
/* alignment.h 		                                                  */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Date: 07-22-2015                                                   */
/* Basic functions/variables commonly used.                           */
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

typedef enum { true, false } bool;

#define GAP 					-3.0
#define MATCH 					 2.0
#define MISMATCH 				-0.5
#define EXTENSION               -1.0
#define JUMP_GENE               -15.0
#define JUMP_EXON               -10.0

// POINTER STATE
#define LEFT 					100
#define DIAGONAL 				200
#define RIGHT	 				300
#define HOME                    400
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
 * change string to upper case
 */
static inline char 
*str_toupper(char* s){
	char *r = mycalloc(strlen(s), char);
	int i = 0;
	char c;
	while(s[i])
	{
		r[i] = toupper(s[i]);
		i++;
	}
	r[strlen(s)] = '\0';
	return r;
}

/*
 * destory kstring
 */
static inline void 
kstring_destory(kstring_t *ks){
	free(ks->s);
	free(ks);
}

/*
 * read two sequences as str1 and str2 from the fasta file;
 * read junctions sites to opt->sites;
 * 
 */
//static inline void 
//kstring_read(char* fname, kstring_t *str1, kstring_t *str2, opt_t *opt){
//	// input check
//	if(fname == NULL || str1 == NULL || str2 == NULL || opt == NULL) 
//		die("kstring_read: input error");
//	
//	// variables declarision
//	int i, l; gzFile fp; kseq_t *seq;
//	char **tmp_seq = mycalloc(3, char*);
//	char **tmp_comment = mycalloc(3, char*);	
//	// parser fasta
//	fp = gzopen(fname, "r");
//	seq = kseq_init(fp);
//	if(fp == NULL || seq == NULL) die("Can't open %s\n", fname);
//	
//	i = 0; while((l=kseq_read(seq)) >= 0){
//		if(i >= 2) die("input fasta file has more than 2 sequences");
//		tmp_seq[i] = strdup(seq->seq.s);
//		if(seq->comment.s) tmp_comment[i] = strdup(seq->comment.s);
//		i++;
//	}
//	// read sequence
//	if(tmp_seq[0] == NULL || tmp_seq[1] == NULL) die("read_kstring: fail to read sequence");
//	(str1)->s = strdup(tmp_seq[0]); (str1)->l = strlen((str1)->s);
//	(str2)->s = strdup(tmp_seq[1]); (str2)->l = strlen((str2)->s);
//	// read the junctions sites if opt != NULL and opt->s==ture
//	if(opt != NULL && opt->s == true){
//		if(tmp_comment[1] == NULL) die("fail to read junction sites");
//		kstring_t *tmp = mycalloc(1, kstring_t);
//		tmp->s = strdup(tmp_comment[1]);
//		tmp->l = strlen(tmp->s);
//		int *fields, i, n;
//		fields = ksplit(tmp, '|', &n);
//		opt->sites.size = n;
//		opt->sites.pos = mycalloc(n, int);
//		for (i = 0; i < n; ++i) opt->sites.pos[i] = atoi(tmp->s + fields[i]);
//		if(tmp) kstring_destory(tmp);
//		if(fields) free(fields);
//	}
//	for(; i >=0; i--) if(tmp_seq[i]) free(tmp_seq[i]);	
//	free(tmp_seq);
//	if(tmp_comment) free(tmp_comment);
//	if(seq) kseq_destroy(seq);
//	gzclose(fp);
//}

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
			default:
				break;
			}
	}
	res_ks1->l = cur;
	res_ks2->l = cur;	
	res_ks1->s = strrev(res_ks1->s);
	res_ks2->s = strrev(res_ks2->s);
}

#endif