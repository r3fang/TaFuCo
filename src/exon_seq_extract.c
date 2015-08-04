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
static struct fasta_uthash *FASTA_HT      = NULL;

int main(int argc, char *argv[]) {
	char  *fname = "sample_data/exons.bed";
	if((FASTA_HT=fasta_uthash_load("/Users/rongxinfang/Documents/genome/hg19/hg19.fa")) == NULL) die("main: fasta_uthash_load fails\n");	
	char  *line = NULL;
	char  **fields = NULL;
	int i, j, num;
	size_t len = 0;
	ssize_t read;
	register char *gname = NULL;
	register char *category=NULL;
	register char *chrom = NULL;
	register int start, end;
	register char *strand;
	struct fasta_uthash *s;
	char *seq;
	FILE *fp = fopen(fname, "r");
	if(fp==NULL) die("[%s] can't open %s", __func__, fname);
	while ((read = getline(&line, &len, fp)) != -1) {
		if((fields = strsplit(line, 0, &num))==NULL) continue;
		gname = fields[0];
		category = fields[1];
		chrom = fields[2];
		start = atoi(fields[3]) - EXON_HALF_FLANK;
		end = atoi(fields[4]) + EXON_HALF_FLANK;
		strand = fields[5];
		printf("%s\t%s\t%s\t%d\t%d\t%s\n", gname, category, chrom, start, end, strand);
		len =  end - start;
		if((s = find_fasta(FASTA_HT, chrom))==NULL) continue;
		seq = mycalloc(len + 2, char);
		memset(seq, '\0',len + 2);
		memcpy(seq, &s->seq[start], len);
		if(strcmp(strand, "-") == 0) seq = rev_com(seq);
		printf("%s\n", seq);
		if(seq)  free(seq);
	}
	fclose(fp);
	if(FASTA_HT) fasta_uthash_destroy(&FASTA_HT);
	if (line) free(line);
	if(gname) free(gname);
	if(strand)   free(strand);
	if(category) free(category);
	if(chrom)    free(chrom);
	return 0;
}