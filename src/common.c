#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include <assert.h>
#include "uthash.h"
#include "kseq.h"
#include "common.h"
KSEQ_INIT(gzFile, gzread)  

/* Reverse complement of DNA-seq */
char*
rev_com(char *s){
	/* Reverse complement of given DNA sequence*/
	assert(s != NULL);
	int n, c, d;
	n=strlen(s);
	char* r;
	r = (char*)malloc(n * sizeof(char));
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

struct kmer_uthash*
find_kmer(char* quary_kmer, struct kmer_uthash *tb) {
    struct kmer_uthash *s;
    HASH_FIND_STR(tb, quary_kmer, s);  /* s: output pointer */
    return s;
}

struct fasta_uthash*
find_fasta(char* quary_name, struct fasta_uthash *tb) {
    struct fasta_uthash *s;
    HASH_FIND_STR(tb, quary_name, s);  /* s: output pointer */
    return s;
}

void 
kmer_uthash_destroy(struct kmer_uthash *table) {
	/*free the kmer_hash table*/
  struct kmer_uthash *cur, *tmp;
  HASH_ITER(hh, table, cur, tmp) {
      HASH_DEL(table, cur);  /* delete it (users advances to next) */
      free(cur);            /* free it */
    }
}

void
kmer_uthash_display(struct kmer_uthash *_kmer_ht) {	
   	struct kmer_uthash *cur, *tmp;
	HASH_ITER(hh, _kmer_ht, cur, tmp) {
		if(cur == NULL)
			exit(-1);
		printf("kmer=%s\tcount=%d\n", cur->kmer, cur->count);
		int i;
		for(i=0; i < cur->count; i++){
			printf("%s\t", cur->pos[i]);
		}
		printf("\n");
	}	
}

void
fasta_uthash_display(struct fasta_uthash *_fasta_ht) {	
   	struct fasta_uthash *cur, *tmp;
	HASH_ITER(hh, _fasta_ht, cur, tmp) {
		if(cur == NULL)
			exit(-1);
		printf(">%s\n%s\n", cur->name, cur->seq);
	}	
}

void 
fasta_uthash_destroy(struct fasta_uthash *table) {
	/*free the kmer_hash table*/
  struct fasta_uthash *cur, *tmp;
  HASH_ITER(hh, table, cur, tmp) {
      HASH_DEL(table, cur);  /* delete it (users advances to next) */
      free(cur);            /* free it */
    }
}


char* 
pos_parser(char *str, int *i) {
	assert(str != NULL);
	size_t size = strsplit_size(str, "_");
	if(size!=2){
		return NULL;
	}
	char **parts = calloc(size, sizeof(char *));
	assert(parts);
	strsplit(str, size, parts, "_");   
	if(parts[0]==NULL || parts[1]==NULL){
		free(parts);
		return NULL;
	}
	*i = atoi(parts[1]);
	char *exon = (char*)malloc(500 * sizeof(char));
	strncpy(exon, parts[0], strlen(parts[0]));	
	free(parts);
	return exon;
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
	int n;
	n=strlen(s);
	char* r;
	r = (char*)malloc(n * sizeof(char));
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

struct kmer_uthash *kmer_uthash_load(char *fname, int *k){	
	/* load kmer_htable*/
	struct kmer_uthash *htable = NULL;
	gzFile fp;  
	kseq_t *seq; 
	int l;
	fp = gzopen(fname, "r");
	if (fp == NULL) {
	  fprintf(stderr, "Can't open input file %s!\n", fname);
	  exit(1);
	}
	seq = kseq_init(fp); // STEP 3: initialize seq  
	if (seq == NULL){
		return NULL;
	}
	while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence 
		/* add kmer */
		char *kmer = seq->name.s;
		*k = strlen(kmer);
		struct kmer_uthash *s;
		s = (struct kmer_uthash*)malloc(sizeof(struct kmer_uthash));
		strncpy(s->kmer, kmer, *k);
		s->kmer[*k] = '\0'; /* just in case*/
		
		/* add count */
		int count = atoi(seq->comment.s);
		s->count = count;
		
		/* add pos */
		char *pos = seq->seq.s;				
		s->pos = malloc((s->count) * sizeof(char*));
		/* split a string by delim */
		char *token;	    
		/* get the first token */
	    token = strtok(pos, "|");				
		/* walk through other tokens */
		int i = 0;
	    while(token != NULL) 
	    {
			s->pos[i] = malloc((strlen(token)+1) * sizeof(char));
			/*duplicate a string*/
			s->pos[i] = strdup(token);
			token = strtok(NULL, "|");
			i ++;
	    }		
		HASH_ADD_STR(htable, kmer, s);
	}	
	gzclose(fp); // STEP 6: close the file handler  
	kseq_destroy(seq);
	return(htable);
}

struct fasta_uthash* fasta_uthash_load(char *fname){
	gzFile fp;
	kseq_t *seq;
	int l;
	struct fasta_uthash *res = NULL;
	fp = gzopen(fname, "r");
	if(fp == NULL){
		return NULL;
	}
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0){
		struct fasta_uthash *s;
		s = (struct fasta_uthash*)malloc(sizeof(struct fasta_uthash));
		/* if we kseq_destroy(seq);, we need to duplicate the string!!!*/
		s->name = strdup(seq->name.s);
		s->seq = strdup(strToUpper(seq->seq.s));
		if(seq->comment.l) s->comment = seq->comment.s;
		HASH_ADD_STR(res, name, s);
	}	
	kseq_destroy(seq);
	gzclose(fp);
	return res;
}