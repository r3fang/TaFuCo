#include <stdio.h>  
#include <stdlib.h>  
#include <string.h> 
#include <zlib.h>  
#include <assert.h>
#include <math.h>
#include <regex.h>
#include "kseq.h"
#include "utils.h"
#include "junction.h"
#include "alignment.h"


int main(int argc, char *argv[]) {
	char *junction_str = "TGAAGCCTGGAGGTGAAACC";
	char * str2 = "CTGGTGGAAACTAGTTATTTATGCCATGTGGAGAGCCAGTGAGATAGATAGATAGTCTGTTTGTTTTGAGGACTTGGAAAGTTGTTCCTATGAAGCCTGGAGGTGAAACCCTCAAAGTCCGCTACTGGCCTCGGGACAGTTGGCCCGTGGGGCTGCCCTACGTGGAAATCCGGGGTGATGACAAGGACTGCTGAGACCTG";
	char *fq1 = "sample_data/reads_simulated_R1.fq";
	char *fq2 = "sample_data/reads_simulated_R2.fq";
	opt_t *opt = init_opt();
	int mismatch = 2;
	gzFile fp1, fp2;
	int l1, l2;
	kseq_t *seq1, *seq2;
	char *_read1, *_read2;
	solution_t *sol1, *sol2;
	
	if((fp1  = gzopen(fq1, "r")) == NULL)   die("[%s] fail to read fastq files\n",  __func__);
	if((fp2  = gzopen(fq2, "r")) == NULL)   die("[%s] fail to read fastq files\n",  __func__);	
	if((seq1 = kseq_init(fp1))   == NULL)   die("[%s] fail to read fastq files\n",  __func__);
	if((seq2 = kseq_init(fp2))   == NULL)   die("[%s] fail to read fastq files\n",  __func__);	
	
	while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >= 0 ) {
		char *_read1 = rev_com(seq1->seq.s); // reverse complement of read1
		char *_read2 = strdup(seq2->seq.s);		
		if(_read1 == NULL || _read2 == NULL) die("[%s] fail to get _read1 and _read2\n", __func__);
		if(strcmp(seq1->name.s, seq2->name.s) != 0) die("[%s] read pair not matched\n", __func__);		
		if((min_mismatch(_read1, junction_str)) <= mismatch ){
			printf("%s\n%s\n", _read1, str2);
			sol1 = align_with_no_jump(_read1, str2, opt);
			printf("%s\n%s\n", sol1->s1, sol1->s2);
			solution_destory(sol1);
		}
		if((min_mismatch(_read2, junction_str)) <= mismatch ){
			printf("%s\n%s\n", _read2, str2);
			sol2 = align_with_no_jump(_read2, str2, opt);
			printf("%s\n%s\n", sol2->s1, sol2->s2);
			solution_destory(sol2);
		}		
		if(_read1) free(_read1);
		if(_read2) free(_read2);
	}
	kseq_destroy(seq1);
	kseq_destroy(seq2);	
	gzclose(fp1);
	gzclose(fp2);
	destory_opt(opt);
	return 0;
}
