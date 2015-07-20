/*--------------------------------------------------------------------*/
/* GlobalJumpAligner.c                                                */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Predict Gene Fusion by given fastq files.                          */
/*--------------------------------------------------------------------*/
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

#define GL_ERR_NONE 			0
#define GAP 					-1.0
#define MATCH 					2.0
#define MISMATCH 				-0.5

typedef enum { true, false } bool;

typedef struct {
  char *a;
  unsigned int alen;
  char *b;
  unsigned int blen;
} seq_pair;
typedef seq_pair *seq_pair_t;

typedef struct {
  double score;
  unsigned int prev[2];
} entry;
typedef entry *entry_t;

typedef struct {
  unsigned int m;
  unsigned int n;
  entry_t **mat;
} matrix;
typedef matrix *matrix_t;

matrix_t create_matrix(size_t m, size_t n){
	if(m < 0 || n < 0) die("create_matrix: parameter error\n");
	size_t i, j; 
	matrix_t S = mycalloc(1, matrix);
	S->m = m;
	S->n = n;
	S->mat = mycalloc(m*n, entry_t);
    
	for (i = 0; i < m; i++) {
      S->mat[i] = mycalloc(n, entry_t);
    }

	for(i = 0; i < S->m; i++){
		for(j = 0; j < S->n; j++)
			S->mat[i][j] = mycalloc(1, entry);
	}
	return S;
}

int destory_matrix(matrix_t S){
	if(S == NULL) die("destory_matrix: parameter error\n");
	int i, j;
	for(i = 0; i < S->m; i++){
		for(j = 0; j < S->n; j++){
			free(S->mat[i][j]);
		}
	}
	free(S);
	return GL_ERR_NONE;
}

int destory_seq_pair(seq_pair_t p){
	if(p == NULL) die("destory_seq_pair: parameter error\n");
	free(p->a);
	free(p->b);
	free(p);
	return GL_ERR_NONE;
}

double max3(double a1, double a2, double a3){
	double res = DBL_MIN;
    res = (a1 >= res) ? a1 : res;
    res = (a2 >= res) ? a2 : res;	
    res = (a3 >= res) ? a3 : res;	
	return res;
}

double smith_waterman(seq_pair *problem, bool local){
	size_t m   = problem->alen + 1;
	size_t n   = problem->alen + 1;
	matrix_t S = create_matrix(m, n);
	size_t i, j, k, l;
	S->mat[0][0]->score		=	0;
	
	for(i=0; i < S->m; i++) S->mat[i][0]->score = 0.0;
	for(j=0; j < S->n; j++) S->mat[0][j]->score = 0.0;
	
	for(i = 1; i <= problem->alen; i++){
		for(j = 1; j <= problem->blen; j++){
	        int new_score = (strncmp(problem->a+(i-1), problem->b+(j-1), 1) == 0) ? MATCH : MISMATCH;
			S->mat[i][j]->score = DBL_MIN;
			S->mat[i][j]->score = max3(S->mat[i-1][j-1]->score + new_score, S->mat[i-1][j]->score + GAP, S->mat[i][j-1]->score + GAP);
		}
	}
	double res = S->mat[problem->alen][problem->blen]->score;
	if(destory_matrix(S) != GL_ERR_NONE) die("smith_waterman: fail to destory matrix");
	return res;
}
/* main function. */
int main(int argc, char *argv[]) {
	char *str1 = "AAAATCTTATCATCTATCATCTACTAAAAAAgggggggggggggggggggggCAGCTAATCGATATTATAAAAAAA";
	char *str2 = "CATCTACTAAAAAAgggggggggggggggggggggCAGCT";	
    seq_pair *problem;
    problem->a = str1; problem->alen = strlen(problem->a);
    problem->b = str2; problem->blen = strlen(problem->b);	
	printf("score=%f\n", smith_waterman(problem, false));
	if(destory_seq_pair(problem) != GL_ERR_NONE) die("destory_seq_pair fails\n");
	return 0;
}
