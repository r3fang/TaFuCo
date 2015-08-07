#include "name2fasta.h"

struct fasta_uthash *extract_exon_seq(char* fname, char *fname_db, struct fasta_uthash *HG19_HT){
	if(fname == NULL || fname_db == NULL || HG19_HT == NULL) return NULL;
	struct fasta_uthash *s_fasta, *cur_fasta, *ret_fasta = NULL;
	str_ctr *s_ctr, *ctr = NULL, *gene_name_ctr = NULL;
	char  *line = NULL;
	size_t len = 0;
	int l;
	ssize_t read;
	char  **fields = NULL;
	int i, j, num;
	register char *gname = NULL;
	register char *category=NULL;
	register char *chrom = NULL;
	register int start, end;
	register char *strand;
	register char *exon_name = NULL;
	struct fasta_uthash *s;
	char *seq;
	char exon_idx[50];
	
	FILE * fp0 = fopen (fname, "r");
	if(fp0==NULL) die("[%s] can't open %s", __func__, fname); 

	while ((read = getline(&line, &len, fp0)) != -1) {
		if((fields = strsplit(line, 0, &num))==NULL) continue; // get rid of \n 
		str_ctr_add(&gene_name_ctr, strToUpper(fields[0]));		
	}
	fclose(fp0);

	FILE *fp = fopen(fname_db, "r");
	if(fp==NULL) die("[%s] can't open %s", __func__, fname_db);
	str_ctr *ctr_s, *ctr_tmp;
	while ((read = getline(&line, &len, fp)) != -1) {
		// get information of exons
		if((fields = strsplit(line, 0, &num))==NULL) continue;
		if(num < 7) continue;
		if((chrom = fields[0])==NULL) continue;
		if((category = fields[2])==NULL) continue;
		if((start = atoi(fields[3]))<0) continue;
		if((end = atoi(fields[4]))<0) continue;
		if((strand = fields[5])==NULL) continue;
		if((gname = strToUpper(fields[6]))==NULL) continue;
		if(strcmp(category, "exon")!=0) continue;
		if((ctr_s=find_str_ctr(gene_name_ctr, gname)) == NULL) continue; // only for targetted genes
		ctr_s->SIZE++;
		if((end - start)<=0) continue;
		// counting exon index of gene
		str_ctr_add(&ctr, gname);
		//// get sequence
		if((s = find_fasta(HG19_HT, chrom))==NULL) continue;
		l =  end - start;
		s_ctr = find_ctr(ctr, gname);
		//// add to FASTA_HT
		sprintf(exon_idx, "%zu", s_ctr->SIZE);
		exon_name = concat(concat(gname, "."), exon_idx);
		if((s_fasta = find_fasta(ret_fasta, exon_name)) == NULL){
			s_fasta = mycalloc(1, struct fasta_uthash);
			s_fasta->name = exon_name;
			s_fasta->chrom = chrom;
			s_fasta->start = start;
			s_fasta->end = end;
			s_fasta->seq = mycalloc(l+1, char);
			memset(s_fasta->seq, '\0',l+1);	
			memcpy(s_fasta->seq, &s->seq[start], l);
			if(strcmp(strand, "-") == 0) s_fasta->seq = rev_com(s_fasta->seq);	
			s_fasta->l = l;
			HASH_ADD_STR(ret_fasta, name, s_fasta);
		}
	}
	HASH_ITER(hh, gene_name_ctr, ctr_s, ctr_tmp) {
		if(ctr_s->SIZE == 1) fprintf(stderr,"%s not found\n", ctr_s->KEY);
	}
	fclose(fp);
	if(gname) free(gname);
	if (line) free(line);
	if(strand)   free(strand);
	if(category) free(category);
	return ret_fasta;
}

int name2fasta_usage(){
	fprintf(stderr, "\n");
			fprintf(stderr, "Usage:   tfc name2fasta [options] <gname.txt> <in.fa.gz> <exons.fa> \n\n");
			fprintf(stderr, "Details: name2fasta is to extract genomic sequence of gene candiates\n\n");
			fprintf(stderr, "Options: -g          organism; 0 for human and 1 for mouse\n\n");
			fprintf(stderr, "Inputs:  gname.txt   plain txt file contains names of genes candiates\n");
			fprintf(stderr, "         in.fa       fasta file contains the entire genome sequence [hg19.fa]\n");
			fprintf(stderr, "         exon.fa     output fasta files that contains sequences of targeted genes\n");
			return 1;
}

int name2fasta(int argc, char *argv[]) {
	int c, i;
	srand48(11);
	char *gene_name, *gff_name, *oname, *iname;
	int orgsm = -1;
	while ((c = getopt(argc, argv, "g:c")) >= 0) {
				switch (c) {
				case 'g': orgsm = atoi(optarg); break;
				default: return 1;
		}
	}

	if (optind + 3 > argc) return name2fasta_usage();
	gene_name = argv[optind];
	iname = argv[optind+1];
	oname = argv[optind+2];

	if(orgsm==-1) {
		fprintf(stderr, "[%s] organism missing\n", __func__);
		return 1;
	}
		
	if(orgsm != 0 && orgsm != 1) {
		fprintf(stderr, "[%s] unrecognized organism, must be 0 or 1\n", __func__);
		return 1;
	}
	
	if(orgsm==0) gff_name = "data/hg.bed";
	if(orgsm==1) gff_name = "data/mm.bed";	
	
	struct fasta_uthash *GENO_HT = NULL;
	struct fasta_uthash *EXON_HT = NULL;
	fprintf(stderr, "[%s] loading reference genome sequences ... \n",__func__);
	
	if((GENO_HT = fasta_uthash_load(iname)) == NULL) die("[%s] can't load reference genome %s", __func__, iname);	
	printf("%s\t%s\n", gene_name, gff_name);
	
	fprintf(stderr, "[%s] extracting targeted gene sequences ... \n",__func__);
	if((EXON_HT = extract_exon_seq(gene_name, gff_name, GENO_HT))==NULL) die("[%s] can't extract exon sequences of %s", __func__, gene_name);

	fprintf(stderr, "[%s] writing down sequences ... \n",__func__);
	if((fasta_uthash_write(EXON_HT, oname))!=0) die("[%s] can't write down to %s", __func__, oname);

	fprintf(stderr, "[%s] cleaning up ... \n", __func__);	

	if(EXON_HT)   fasta_uthash_destroy(&EXON_HT);
	if(GENO_HT)   fasta_uthash_destroy(&GENO_HT);    
	return 0;
}
