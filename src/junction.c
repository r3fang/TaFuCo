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
	char* str2 = "GAAGCCTGAGTTGGCAGAGGAGGAGCCAGAAGAGAACTGCCCTTTCTGCACTGGTGGAAACTAGTTATTTATGCCATGTGGAGAGCCAGTGAGATAGATAGATAGTCTGTTTGTTTTGAGGACTTGGAAAGTTGTTCCTATGAAGCCTGGAGGTGAAACCCTCAAAGTCCGCTACTGGCCTCGGGACAGTTGGCCCGTGGGGCTGCCCTACGTGGAAATCCGGGGTGATGACAAGGACTGCTGAGACCTGCCAGCCTCTTGACTATCCAACCACCCAACAACTTGCCTGTGTCTTGAGAAG";
	char *fq1 = "sample_data/reads_simulated_R1.fq";
	char *fq2 = "sample_data/reads_simulated_R2.fq";
	opt_t *opt = init_opt();

	int mismatch = 2;
	gzFile fp1, fp2;
	register int l1, l2;
	register kseq_t *seq1, *seq2;
	register char *_read1, *_read2;
	register solution_t *sol1, *sol2;
	if((fp1  = gzopen(fq1, "r")) == NULL)   die("[%s] fail to read fastq files\n",  __func__);
	if((fp2  = gzopen(fq2, "r")) == NULL)   die("[%s] fail to read fastq files\n",  __func__);	
	if((seq1 = kseq_init(fp1))   == NULL)   die("[%s] fail to read fastq files\n",  __func__);
	if((seq2 = kseq_init(fp2))   == NULL)   die("[%s] fail to read fastq files\n",  __func__);	
	
	while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >= 0 ) {
		_read1 = rev_com(seq1->seq.s); // reverse complement of read1
		_read2 = seq2->seq.s;		
		if(_read1 == NULL || _read2 == NULL) die("[%s] fail to get _read1 and _read2\n", __func__);
		if(strcmp(seq1->name.s, seq2->name.s) != 0) die("[%s] read pair not matched\n", __func__);		
		
		if((min_mismatch(_read1, junction_str)) <= mismatch ){
			sol1 = align_with_no_jump(_read1, str2, opt);
			printf("prob=%f\n", sol1->prob);
			free(sol1);
		}
		if((min_mismatch(_read2, junction_str)) <= mismatch){
			sol2 = align_with_no_jump(_read2, str2, opt);
			printf("prob=%f\n", sol2->prob);
			free(sol2);
		}
	}
	kseq_destroy(seq1);
	kseq_destroy(seq2);	
	gzclose(fp1);
	gzclose(fp2);
	destory_opt(opt);
	return 0;
}
