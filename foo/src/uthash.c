#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include "uthash.h"
#include "utlist.h"
#include "utstring.h"
#include "utarray.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)  

struct kmer_uthash {
    char kmer[15];                /* key */
    int pos;                    
    UT_hash_handle hh;         /* makes this structure hashable */
};

char* concat(char *s1, char *s2)
{
    char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the zero-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
	strcat(result, s2);
	return result;
}

void add_kmer(struct kmer_uthash **table, char kmer[15], int pos) {
	/* You really need to pass a pointer to the hash pointer: **table*/
	struct kmer_uthash *s;
	HASH_FIND_STR(*table, kmer, s);
	if (s==NULL){
		s = (struct kmer_uthash*)malloc(sizeof(struct kmer_uthash));
		strncpy(s->kmer, kmer, 15);
		s->pos = pos;
	}
	HASH_ADD_STR(*table, kmer, s);
}
	
int test(char *fasta_file) {
	struct kmer_uthash *table = NULL;
	gzFile fp;  
	kseq_t *seqs;  
	int l;
	int k=15;
	printf("%s\n", fasta_file);
	fp = gzopen(fasta_file, "r");
	seqs = kseq_init(fp); // STEP 3: initialize seq  
	while ((l = kseq_read(seqs)) >= 0) { // STEP 4: read sequence 
		char *seq = seqs->seq.s;
		int i;
		if (seqs->name.s==NULL)
			return NULL;
		char *name = seqs->name.s;
		char kmer[k];
		for(i=0; i < strlen(seq)-k+1; i++){
			memcpy(kmer, &seq[i], k);
			kmer[k] = '\0';
			/* convert i to string */
			char i_str[100];
			sprintf(i_str, "%d", i);
			add_kmer(&table, kmer, i); /* You really need to pass a pointer to the hash pointer: */
			//printf("%s\t%s\n", kmer, concat(concat(name, "."), i_str));
		}
		break;
	}  
	struct kmer_uthash *s, *tmp;
	HASH_ITER(hh, table, s, tmp) {
	    printf("kmer %s: pos %d\n", s->kmer, s->pos);
	}
	
	kseq_destroy(seqs);
	gzclose(fp);
	return 0;
}
