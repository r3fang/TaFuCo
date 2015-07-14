#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h> 
#include <assert.h>
#include "kseq.h" 
#include "uthash.h"

#define MAX_K 100

struct kmer_uthash {
    char kmer[MAX_K];                /* key */
	int count;
	char **pos;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

struct fasta_uthash {
    char* name;                /* key */
	char* seq;
    UT_hash_handle hh;         /* makes this structure hashable */
};

struct MPM
{
	char* READ_NAME;
	char* EXON_NAME;	
	int READ_POS;
	int EXON_POS;
	int LENGTH;
};

struct fasta_uthash *fasta_uthash_load(char *);
struct kmer_uthash *kmer_uthash_load(char *, int *);
void MPM_display(struct MPM *);
char* rev_com(char *s);
char* small_dna_str(char*, char*);
int score_DNA_str(char *s);
struct kmer_uthash *find_kmer(char* , struct kmer_uthash*);
struct fasta_uthash *find_fasta(char* , struct fasta_uthash*);
void kmer_uthash_display(struct kmer_uthash*);	
void fasta_uthash_display(struct fasta_uthash*);
void fasta_uthash_destroy(struct fasta_uthash*);
void kmer_uthash_destroy(struct kmer_uthash*);
char* concat(char*, char*);
char* strToUpper(char*);
int strsplit (const char *str, size_t size, char *parts[], const char *delimiter);
int strsplit_size (const char *str, const char *delimiter);
char *pos_parser(char *str, int *i);

#endif