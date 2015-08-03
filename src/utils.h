#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <stdarg.h>
#include "zlib.h"
#include "kseq.h"
#include "uthash.h"

KSEQ_INIT(gzFile, gzread);

#ifndef BOOL_DEFINED
#define BOOL_DEFINED
typedef char BOOL ;
#define TRUE 1
#define FALSE 0
#endif

// constant define
typedef enum { true, false } bool;
#define pair(k1, k2)  ((k1 + k2)*(k1 + k2 + 1)/2 + k2)

typedef struct
{
	char 	*KEY;
	size_t  SIZE;
	UT_hash_handle hh;
} str_ctr;

static inline void str_ctr_destory(str_ctr **s){
	str_ctr *cur, *tmp;
	HASH_ITER(hh, *s, cur, tmp) {
		HASH_DEL(*s, cur);  /* delete; users advances to next */
		free(cur);            /* optional- if you want to free  */
	}
}

unsigned long
hash(char *str)
{
    unsigned long hash = 5381;
    int c;

    while (c == *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    return hash;
}

static inline int 
mystrcmp(const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

static inline char* 
concat(char *s1, char *s2)
{
	if(s1 == NULL) return s2;
	if(s2 == NULL) return s1;	
    char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the zero-terminator
    strcpy(result, s1);
	strcat(result, s2);
	return result;
}

/* Reverse complement of DNA-seq */
static inline char 
*rev_com(char *s){
	/* Reverse complement of given DNA sequence*/
	assert(s != NULL);
	int n, c, d;
	n=strlen(s);
	char* r;
	r = (char*)malloc((n+1) * sizeof(char)); // always allocate N+1
	for (c = n - 1, d = 0; c >= 0; c--, d++){
		switch(toupper(s[c])){
			case 'A':
				r[d] = 'T';
				break;
			case 'T':
				r[d] = 'A';
				break;
			case 'C':
				r[d] = 'G';
				break;
			case 'G':
				r[d] = 'C';
				break;
			default:
				r[d] = s[c];
				break;
		}
	}
	r[n] = '\0';
	return r;
}

static inline char
*strToUpper(char* s){
	if(s == NULL)
		return NULL;
	int n;
	n=strlen(s);
	char* r;
	r = (char*)malloc((n+1) * sizeof(char));
	if(r == NULL)
		return NULL;	
	int i;
	for(i=0; i < n; i++){
		r[i] = toupper(s[i]);
	}
	r[n] = '\0';
	return r;
}


static inline void die (char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "FATAL ERROR: ") ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;
  exit (-1) ;
}

static inline void *_mycalloc (long number, int size)
{
  void *p = (void*) calloc (number, size) ;
  if (!p) die ("mycalloc failure requesting %d of size %d bytes", number, size) ;
  return p ;
}
#define mycalloc(n,type) (type*)_mycalloc(n,sizeof(type))

static inline int 
str_ctr_add(str_ctr** tb, char* key){
	if(key == NULL) die("str_ctr_add: prameter error\n");
	str_ctr *s;
	HASH_FIND_STR(*tb, key, s);
	if(s == NULL){
		s = mycalloc(1, str_ctr);
		s->KEY = key;
		s->SIZE = 1;
		HASH_ADD_STR(*tb, KEY, s);
	}else{
		s->SIZE++;
	}
	return 0;
}

static inline int  
str_ctr_sort(str_ctr *a, str_ctr *b) {
    return (a->SIZE >= b->SIZE);
}

//opt
typedef struct {
	char* fq1; // gap open
	char* fq2; // gap open
	char* fa; // gap open
	int k; // gap extension
	int min_match; // match
	int min_weight; // unmatch
	int match;
	int mismatch;
	int gap;
	int extension;
	int jump_gene;
	int jump_exon;
	int min_hits;
	double min_align_score;
} opt_t;

static inline opt_t *init_opt(){
	opt_t *opt = mycalloc(1, opt_t);
	opt->fq1 = NULL;
	opt->fq2 = NULL;
	opt->fa = NULL;
	opt->k = 15;
	opt->min_match = 10;
	opt->min_weight = 3;	
	opt->match = 2;
	opt->mismatch = -2.0;
	opt->gap = -5.0;
	opt->extension = -1.0;
	opt->jump_gene = -10.0;
	opt->jump_exon = -8.0;
	opt->min_hits = 3;
	opt->min_align_score = 0.8;
	return opt;
}

static inline void destory_opt(opt_t *opt){
	if(opt->fq1) free(opt->fq1);
	if(opt->fq2) free(opt->fq2);
	if(opt->fa)  free(opt->fa);
	free(opt);
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

static inline char 
*strrev(char *s){
	if(s == NULL) return NULL;
	int l = strlen(s);
	char* r = mycalloc(l+1, char);
	int i; for(i=0; i<l; i++){
		r[i] = s[l-i-1];
	}
	r[i] = '\0';
	return r;
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
*substr(char* s, int i, int len){
	if(s==NULL) die("[%s] error input", __func__);
	if(i < 0 || i > strlen(s)) die("[%s] error input", __func__);
	if(len < 0 || len > strlen(s)) die("[%s] error input", __func__);
	int j = i + len;
	char *r = mycalloc(len+1, char);
	register int pt;
	for(pt=i; pt<j; pt++){
		r[pt-i] = s[pt];
	}
	r[pt] = '\0';
	return r;
}
#endif