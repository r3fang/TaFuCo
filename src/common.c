#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include <assert.h>
#include "uthash.h"
#include "kseq.h"
#include "common.h"


size_t set_str_arr(char **arr, char **arr_uniq, size_t size){
	if(arr == NULL) {return 0;}
	struct str_count {
	    char* STR;               
		int COUNT;
	    UT_hash_handle hh;
	};
	
	struct str_count *strs = NULL;
	
	int i;
	for(i=0; i<size; i++){
		struct str_count *s;
		HASH_FIND_STR(strs, arr[i], s);
		/* this is a new str*/
		if(s==NULL){
			s = (struct str_count*)malloc(sizeof(struct str_count));
			s->STR = strdup(arr[i]);
			s->COUNT = 1;
			HASH_ADD_STR(strs, STR, s);
		}
	}	
	size_t num = HASH_COUNT(strs);
	// copy to arr_uniq
    struct str_count *s;
	int j = 0;
    for(s=strs; s != NULL; s=s->hh.next) {
        arr_uniq[j] = strdup(s->STR);
		j++;
    }
	return num;		
}

/* Reverse complement of DNA-seq */
char *rev_com(char *s){
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

char*
small_dna_str(char* s1, char* s2){
	assert(s1 != NULL);
	assert(s2 != NULL);
	int score1 = score_DNA_str(s1);
	int score2 = score_DNA_str(s2);
	if(score1 <= score2){
		return s1;
	}else{
		return s2;		
	}
}

int 
score_DNA_str(char *s){
	assert(s != NULL);
	int score = 0;
	int flag;
	int i;
	for(i = 0; i < strlen(s); i++){
		switch(toupper(s[i])){
			case 'A':
				flag = 0;
				break;
			case 'T':
				flag = 1;
				break;
			case 'C':
				flag = 2;
				break;
			case 'G':
				flag = 3;
				break;
			default:
				flag = 5;
				break;
		}
		score = 5*score + flag;
	}
	return score;
}

char* 
concat(char *s1, char *s2)
{
    char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the zero-terminator
    strcpy(result, s1);
	strcat(result, s2);
	return result;
}

char* 
strToUpper(char* s){
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


int max_of_int_array(const int *arr, size_t length) {
    size_t i;
    int minimum = arr[0];
    for (i = 1; i < length; ++i) {
        if (minimum > arr[i]) {
            minimum = arr[i];
        }
    }
    return minimum;
}