#include "name2fasta.h"

static char *name_trim(char *s, char delim){
	printf("%s\n", s);
	if(s==NULL || strlen(s) <= 3) return NULL;
	int num, i;
	char**fields_tmp = NULL;
	char* gname;
	printf("%s\n", s);
	fields_tmp = strsplit(s, delim, &num);
	printf("%s\n", s);
	
	if(num<2){
		//for(i=0; i<num; i++){if(fields_tmp[i]) free(fields_tmp[i]);} //free(fields_tmp);
		return NULL;
	}
	gname = strdup(fields_tmp[0]);
	for(i=0; i<num; i++){
	if(fields_tmp[i]) free(fields_tmp[i]);		
	} 
	free(fields_tmp);
	return gname;
}

static fasta_t *extract_exon_seq(char* fname, char *fname_db, fasta_t *HG19_HT, char* genr){
	//if(fname==NULL || fname_db==NULL || HG19_HT==NULL || genr==NULL) return NULL;
	fasta_t *s_fasta, *cur_fasta, *ret_fasta = NULL;
	str_ctr *s_ctr, *ctr = NULL, *gene_name_ctr = NULL;
	char  *line = NULL;
	size_t len = 0;
	int l;
	ssize_t read;
	char  **fields = NULL;
	
	int i, j, num;
	register char *category=NULL;
	register char *chrom = NULL;
	register char *strand = NULL;
	register char *name = NULL;	
	register char *gene_id = NULL;	
	register char *gene_name = NULL;	
	register char *transcript_id = NULL;	
	register char *tss_id = NULL;		
	register int start, end;
	fasta_t *s;
	char *seq;
	char idx[50];
	
	FILE * fp0 = fopen (fname, "r");
	if(fp0==NULL) die("[%s] can't open %s", __func__, fname); 
	while ((read = getline(&line, &len, fp0)) != -1) {
		if((fields = strsplit(line, 0, &num))==NULL) continue; // get rid of \n 
		str_ctr_add(&gene_name_ctr, strToUpper(fields[0]));		
		for(i=0; i<num; i++) free(fields[i]);
	}
	fclose(fp0);
	
	FILE *fp = fopen(fname_db, "r");
	if(fp==NULL) die("[%s] can't open %s", __func__, fname_db);
	str_ctr *ctr_s, *ctr_tmp;
	while ((read = getline(&line, &len, fp)) != -1) {
		fields=NULL;
		category=chrom=strand=name=gene_id=gene_name=transcript_id=tss_id=NULL;
		if(strlen(line)<20) goto CONTINUE;
		// get information of exons
		printf("%s\n", line);
		if((fields = strsplit(line, 0, &num))==NULL) goto CONTINUE;
		if((chrom = fields[0])==NULL) goto CONTINUE;
		if((category = fields[2])==NULL) goto CONTINUE;
		if((strand = fields[6])==NULL) goto CONTINUE;	
		if((gene_id = name_trim(fields[9], '"'))==NULL) goto CONTINUE;
		if((gene_name = name_trim(fields[11], '"'))==NULL) goto CONTINUE;
		if((transcript_id = name_trim(fields[13], '"'))==NULL) goto CONTINUE;
		if((tss_id = name_trim(fields[15], '"'))==NULL) goto CONTINUE;
		start = atoi(fields[3]);
		end = atoi(fields[4]);
		//printf("%s %s %s %s %s %s %s %d %d\n", chrom, category, strand, gene_id, gene_name, transcript_id, tss_id, start, end);
		//fields_tmp = strsplit(fields[17], '"', &j);
		//if(j!=2){for(i=0; i<j; i++) free(fields_tmp[i]); continue;}
		//gene_name = fields_tmp[0];
		//if(strcmp(category, genr)==0 && (ctr_s=find_str_ctr(gene_name_ctr, gname))!=NULL && ((end - start)>0)){
		//	ctr_s->SIZE++;
		//	str_ctr_add(&ctr, gname);
		//	if((s = find_fasta(HG19_HT, chrom))==NULL) continue;
		//	l =  end - start;
		//	s_ctr = find_ctr(ctr, gname);
		//	sprintf(idx, "%d", s_ctr->SIZE);
		//	name = concat(concat(gname, "."), idx);
		//	if((s_fasta = find_fasta(ret_fasta, name)) == NULL){
		//		s_fasta = mycalloc(1, fasta_t);
		//		s_fasta->name = name;
		//		s_fasta->chrom = chrom;
		//		s_fasta->start = start;
		//		s_fasta->end = end;
		//		s_fasta->seq = mycalloc(l+1, char);
		//		memset(s_fasta->seq, '\0',l+1);	
		//		memcpy(s_fasta->seq, &s->seq[start], l);
		//		if(strcmp(strand, "-") == 0) s_fasta->seq = rev_com(s_fasta->seq);	
		//		s_fasta->l = l;
		//		HASH_ADD_STR(ret_fasta, name, s_fasta);
		//	}
		//}
		CONTINUE:
			if(fields){for(i=0; i<num; i++) free(fields[i]); free(fields);}
			if(gene_name)  free(gene_name);
			if(gene_id)  free(gene_id);
			if(transcript_id)  free(transcript_id);
			if(tss_id)  free(tss_id);
			continue;
		
		if(fields){for(i=0; i<num; i++) free(fields[i]); free(fields);}
		if(gene_name)  free(gene_name);
		if(gene_id)  free(gene_id);
		if(transcript_id)  free(transcript_id);
		if(tss_id)  free(tss_id);
	}
	fclose(fp);
	printf("SFAKSFHKASLFHASKLFHA;\n");
	//HASH_ITER(hh, gene_name_ctr, ctr_s, ctr_tmp) {
	//	if(ctr_s->SIZE == 1) fprintf(stderr,"%s not found\n", ctr_s->KEY);
	//}
	return ret_fasta;
}

int name2fasta_usage(){
	fprintf(stderr, "\n");
			fprintf(stderr, "Usage:   tfc name2fasta [options] <gname.txt> <gencode.gtf> <in.fa.gz> <exons.fa> \n\n");
			fprintf(stderr, "Details: name2fasta is to extract genomic sequence of gene candiates\n\n");
			fprintf(stderr, "Options: -g          'exon' or 'transcript' \n\n");
			fprintf(stderr, "Inputs:  .txt        plain txt file contains names of gene candiates e.g. [genes.txt]\n");
			fprintf(stderr, "         .gtf        gft file that contains gene annotation\n");
			fprintf(stderr, "         .fa         fasta file contains the whole genome sequence   e.g. [hg19.fa.gz]\n");
			fprintf(stderr, "         .fa         output fasta files contains extracted seq of targeted genes\n");
			return 1;
}

int name2fasta(int argc, char *argv[]) {
	int c, i;
	srand48(11);
	char *gene_name, *gff_name, *oname, *iname, *genr;
	while ((c = getopt(argc, argv, "g:c")) >= 0) {
				switch (c) {
				case 'g': genr = optarg; break;
				default: return 1;
		}
	}
	
	if(strcmp(genr, "transcript")!=0 && strcmp(genr, "exon")!=0){
		fprintf(stderr, "-g unrecognized gener 'exon' or 'transcript'\n");
		return -1;
	}
	if (optind + 4 > argc) return name2fasta_usage();
	gene_name = argv[optind];
	gff_name = argv[optind+1];
	iname = argv[optind+2];
	oname = argv[optind+3];
		
	fasta_t *GENO_HT = NULL;
	fasta_t *EXON_HT = NULL;
	
	//fprintf(stderr, "[%s] loading reference genome sequences ... \n",__func__);
	//if((GENO_HT = fasta_read(iname)) == NULL) die("[%s] can't load reference genome %s", __func__, iname);	
	
	fprintf(stderr, "[%s] extracting targeted gene sequences ... \n",__func__);
	if((EXON_HT = extract_exon_seq(gene_name, gff_name, GENO_HT, genr))==NULL) die("[%s] can't extract exon sequences of %s", __func__, gene_name);
	
	//fprintf(stderr, "[%s] writing down sequences ... \n",__func__);
	//if((fasta_write(EXON_HT, oname))!=0) die("[%s] can't write down to %s", __func__, oname);
    //
	//fprintf(stderr, "[%s] cleaning up ... \n", __func__);	
    
	if(EXON_HT)   fasta_destroy(&EXON_HT);
	if(GENO_HT)   fasta_destroy(&GENO_HT);    
	return 0;
}