#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include <assert.h>
#include "uthash.h"
#include "kseq.h"
#include "common.h"
KSEQ_INIT(gzFile, gzread); 

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
pos_parser(char *str, int *i) {
	assert(str != NULL);
	size_t size = strsplit_size(str, "_");
	if(size!=2){
		return NULL;
	}
	char **parts = calloc(size, sizeof(char *));
	assert(parts != NULL);
	strsplit(str, size, parts, "_");   
	if(parts[0]==NULL || parts[1]==NULL){
		free(parts);
		return NULL;
	}
	*i = atoi(parts[1]);
	char *exon;
	exon = strdup(parts[0]);	
	exon[strlen(parts[0])] = '\0'; // double make sure
	free(parts);	
	if(exon != NULL){
		return exon;		
	}
	return NULL;
}

int
strsplit_size (const char *str, const char *delimiter) {
  /* First count parts number*/
  char *pch;
  int i = 0;
  char *tmp = strdup(str);
  pch = strtok(tmp, delimiter);
  /* if there is no delim in the string*/
  i ++;
  while (pch) {
    pch = strtok(NULL, delimiter);
    if (NULL == pch) break;
	i++;
  }
  free(tmp);
  free(pch);
  return i;
}

int
strsplit (const char *str, size_t size, char *parts[], const char *delimiter) {	
	
	if(size==1){
		parts[0] = strdup(str);
		return 0;
	}
  
  char *pch;
  int i = 0;
  char *tmp = strdup(str);
  pch = strtok(tmp, delimiter);
  parts[i++] = strdup(pch);

  while (pch) {
    pch = strtok(NULL, delimiter);
    if (NULL == pch) break;
    parts[i++] = strdup(pch);
  }

  free(tmp);
  free(pch);
  return 0;
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

void 
MPM_display(struct MPM *s){
	if(s==NULL){
		fprintf(stderr, "Error: display an empty MPM");
		exit(-1);
	}else{
		printf("%s\t%d\t%s\t%d\t%d\n", s->READ_NAME, s->READ_POS, s->EXON_NAME, s->EXON_POS, s->LENGTH);		
	}
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