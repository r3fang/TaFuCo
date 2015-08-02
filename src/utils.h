#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include "zlib.h"
#include "kseq.h"
#include "uthash.h"

#if defined(__APPLE__)
#  define COMMON_DIGEST_FOR_OPENSSL
#  include <CommonCrypto/CommonDigest.h>
#  define SHA1 CC_SHA1
#else
#  include <openssl/md5.h>
#endif

KSEQ_INIT(gzFile, gzread);

#ifndef BOOL_DEFINED
#define BOOL_DEFINED
typedef char BOOL ;
#define TRUE 1
#define FALSE 0
#endif

typedef struct
{
	char 	*KEY;
	size_t  SIZE;
	UT_hash_handle hh;
} str_ctr;

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

char *str2md5(const char *str, int length) {
    int n;
    MD5_CTX c;
    unsigned char digest[16];
    char *out = (char*)malloc(33);

    MD5_Init(&c);

    while (length > 0) {
        if (length > 512) {
            MD5_Update(&c, str, 512);
        } else {
            MD5_Update(&c, str, length);
        }
        length -= 512;
        str += 512;
    }

    MD5_Final(digest, &c);

    for (n = 0; n < 16; ++n) {
        snprintf(&(out[n*2]), 16*2, "%02x", (unsigned int)digest[n]);
    }

    return out;
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
	opt->mismatch = -2;
	opt->gap = -5;
	opt->extension = -1;
	opt->jump_gene = -10;
	opt->jump_exon = -8;
	opt->min_hits = 3;
	opt->min_align_score = 0.8;
	return opt;
}
#endif