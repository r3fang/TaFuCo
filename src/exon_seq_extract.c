/*--------------------------------------------------------------------*/
/* exon_seq_extract.c                                                 */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Predict Gene Fusion by given fastq files.                          */
/*--------------------------------------------------------------------*/
#include <stdio.h>  
#include <stdlib.h>  
#include <string.h> 
#include <zlib.h>  
#include <assert.h>
#include <math.h>
#include <regex.h>
#include "kseq.h"
#include "kstring.h"
#include "utils.h"
#include "fasta_uthash.h"

#define EXON_HALF_FLANK                     300

int main(int argc, char *argv[]) {
	struct fasta_uthash *HG19_HT       = NULL;
	struct fasta_uthash *s_fasta, *cur_fasta, *FASTA_HT = NULL;
	str_ctr *s_ctr, *ctr = NULL, *gene_name_ctr = NULL;
	char  *line = NULL;
	size_t len = 0;
	ssize_t read;
	char  **fields = NULL;
	int i, j, num;
	register char *gname = NULL;
	register char *category=NULL;
	register char *chrom = NULL;
	register int start, end;
	register char *strand;
	struct fasta_uthash *s;
	char *seq;
	char exon_idx[10];
	char  *fname = "/bioinfoSD/users/rfang/output";
	
	FILE *fp0 = fopen("sample_data/genes.name.txt", "r");
	if(fp0==NULL) die("[%s] can't open %s", __func__, "sample_data/genes.name.txt"); 
	while ((read = getline(&line, &len, fp0)) != -1) {
		if((fields = strsplit(line, 0, &num))==NULL) continue; // get rid of \n 
		str_ctr_add(&gene_name_ctr, fields[0]);		
	}
	fclose(fp0);
	
	if((HG19_HT=fasta_uthash_load("/bioinfoSD/users/rfang/data/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa")) == NULL) die("main: fasta_uthash_load fails\n");	
	FILE *fp = fopen(fname, "r");
	if(fp==NULL) die("[%s] can't open %s", __func__, fname);
	while ((read = getline(&line, &len, fp)) != -1) {
		// get information of exons
		if((fields = strsplit(line, 0, &num))==NULL) continue;
		if(num < 7) continue;
		if((chrom = fields[0])==NULL) continue;
		if((category = fields[2])==NULL) continue;
		if((start = atoi(fields[3]) - EXON_HALF_FLANK)<0) continue;
		if((end = atoi(fields[4]) + EXON_HALF_FLANK)<0) continue;
		if((strand = fields[5])==NULL) continue;
		if((gname = fields[6])==NULL) continue;
		if(strcmp(category, "exon")!=0) continue;
		if((find_str_ctr(gene_name_ctr, gname)) == NULL) continue; // only for targetted genes
		// counting exon index of gene
		str_ctr_add(&ctr, gname);
		// get sequence
		len =  end - start;
		if((s = find_fasta(HG19_HT, chrom))==NULL) continue;
		seq = mycalloc(len + 2, char);
		memset(seq, '\0',len + 2);
		end = min(end, strlen(s->seq));
		start = max(start, 0);
		memcpy(seq, &s->seq[start], len);
		if(strcmp(strand, "-") == 0) seq = rev_com(seq);
		s_ctr = find_ctr(ctr, gname);
		// add to FASTA_HT
		sprintf(exon_idx, "%d", s_ctr->SIZE);
		printf("%s\n", concat(concat(gname, "."), exon_idx));
		if(seq)  free(seq);
	}
	fclose(fp);
	if(gname) free(gname);
	if(HG19_HT) fasta_uthash_destroy(&HG19_HT);
	if (line) free(line);
	if(strand)   free(strand);
	if(category) free(category);
	if(chrom)    free(chrom);
	return 0;
}