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

#endif