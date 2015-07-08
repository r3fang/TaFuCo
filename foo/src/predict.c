#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <zlib.h>  
#include "uthash.h"
#include "utlist.h"
#include "utstring.h"
#include "utarray.h"
#include "kseq.h"
#include "common.h"

KSEQ_INIT(gzFile, gzread)  

int predict_main(char *fasta_file){
	/* main function*/
	printf("%s\n", fasta_file);
	gzFile fp;  
	kseq_t *seq;  
	int l;
	fp = gzopen(fasta_file, "r");
	if (fp == NULL) {
	  fprintf(stderr, "Can't open input file %s!\n", fasta_file);
	  exit(1);
	}
	seq = kseq_init(fp); // STEP 3: initialize seq  
	while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence 
		printf("name: %s\n", seq->name.s); 
		printf("seq: %s\n", seq->seq.s); 
	}
	kseq_destroy(seq); // STEP 5: destroy seq  
	gzclose(fp); // STEP 6: close the file handler  
	char *index_file = concat(fasta_file, ".index");
	printf("%s\n", index_file);
	return 0;
}