#include "name2fasta.h"

static char *name_trim(char *s, char delim){
	if(s==NULL || strlen(s) <= 3) return NULL;
	int num, i;
	char**fields_tmp = NULL;
	char* gname;
	fields_tmp = strsplit(s, delim, &num);	
	if(num<2){
		return NULL;
	}
	gname = strdup(fields_tmp[0]);
	for(i=0; i<num; i++){
	if(fields_tmp[i]) free(fields_tmp[i]);		
	} 
	free(fields_tmp);
	return gname;
}

static fasta_t *extract_exon_seq(char* fname, char *fname_db, fasta_t *HG19_HT){
	if(fname==NULL || fname_db==NULL || HG19_HT==NULL) return NULL;
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
		name=category=chrom=strand=name=gene_id=gene_name=transcript_id=tss_id=NULL;
		if(strlen(line)<15) goto CONTINUE;
		// get information of exons
		if((fields = strsplit(line, 0, &num))==NULL) goto CONTINUE;
		if(num<15) goto CONTINUE;
		
		if((chrom = fields[0])   ==NULL) goto CONTINUE;
		if((category = fields[2])==NULL) goto CONTINUE;
		if((strand = fields[6])  ==NULL) goto CONTINUE;		
		for(i=0; i<num; i++){
			if(strcmp(fields[i], "gene_id")==0)             gene_id=strToUpper(name_trim(fields[i+1], '"'));
			if(strcmp(fields[i], "gene_name")==0)         gene_name=strToUpper(name_trim(fields[i+1], '"'));			
			if(strcmp(fields[i], "transcript_id")==0) transcript_id=name_trim(fields[i+1], '"');			
			if(strcmp(fields[i], "tss_id")==0)               tss_id=name_trim(fields[i+1], '"');			
		}
		
		if(gene_id==NULL || gene_name==NULL || transcript_id==NULL || tss_id==NULL) goto CONTINUE;
		start = atoi(fields[3]);
		end   = atoi(fields[4]);
		
		if(strcmp(category, "exon")==0 && ((ctr_s=find_str_ctr(gene_name_ctr, gene_name))!=NULL) && ((end - start)>0)){
			ctr_s->SIZE++;
			str_ctr_add(&ctr, gene_name);
			if((s = find_fasta(HG19_HT, chrom))==NULL) goto CONTINUE;
			l =  end - start;
			s_ctr = find_ctr(ctr, gene_name);
			name = join(7, chrom, ".", fields[3], ".", fields[4], ".", gene_name);
			if((s_fasta = find_fasta(ret_fasta, name)) == NULL){
				s_fasta = fasta_init();
				s_fasta->name = strdup(name);
				s_fasta->idx = s_ctr->SIZE;
				s_fasta->chrom = strdup(chrom);
				s_fasta->gene_name = strdup(gene_name);
				s_fasta->gene_id = strdup(gene_id);

				s_fasta->transcript_id = mycalloc(1, char*);
				s_fasta->tss_id = mycalloc(1, char*);

				s_fasta->transcript_id[0] =strdup(transcript_id);
				s_fasta->tss_id[0] = strdup(tss_id);
				
				s_fasta->transcript_num = 1;
				s_fasta->tss_num = 1;
				
				s_fasta->strand = strdup(strand);
				s_fasta->start = start;
				s_fasta->end = end;
				s_fasta->seq = mycalloc(l+1, char);
				memset(s_fasta->seq, '\0',l+1);	
				memcpy(s_fasta->seq, &s->seq[start], l);
				if(strcmp(strand, "-") == 0) s_fasta->seq = rev_com(s_fasta->seq);	
				s_fasta->l = l;
				HASH_ADD_STR(ret_fasta, name, s_fasta);
			}else{
				s_fasta->transcript_num++;
				s_fasta->tss_num++;
				s_fasta->transcript_id = realloc(s_fasta->transcript_id, s_fasta->transcript_num * sizeof(*s_fasta->transcript_id));
				s_fasta->tss_id = realloc(s_fasta->tss_id, s_fasta->tss_num * sizeof(*s_fasta->tss_id));
				s_fasta->transcript_id[s_fasta->transcript_num-1] = strdup(transcript_id);
				s_fasta->tss_id[s_fasta->tss_num-1] = strdup(tss_id);
			}
		}
		/* skip current item */
		CONTINUE:
			if(fields){for(i=0; i<num; i++) free(fields[i]); free(fields);}
			if(gene_name)      free(gene_name);
			if(gene_id)        free(gene_id);
			if(transcript_id)  free(transcript_id);
			if(tss_id)         free(tss_id);
			if(name)           free(name);
			continue;
		/* free all items */
		if(fields){for(i=0; i<num; i++) free(fields[i]); free(fields);}
		if(gene_name)      free(gene_name);
		if(gene_id)        free(gene_id);
		if(transcript_id)  free(transcript_id);
		if(tss_id)         free(tss_id);
		if(name)           free(name);
	}
	fclose(fp);
	
	HASH_ITER(hh, gene_name_ctr, ctr_s, ctr_tmp) {
		if(ctr_s->SIZE == 1) fprintf(stderr,"%s not found\n", ctr_s->KEY);
	}
	return ret_fasta;
}

static int fasta_write_exon(fasta_t *fa, char* fname){
	if(fa==NULL) return -1;
	FILE *fp = fopen(fname, "w");
	fasta_t *s;
	int i;
	if(fp==NULL) die("[%s] can't open %s", __func__, fname);	
	for(s=fa; s!=NULL; s=s->hh.next){
		fprintf(fp, ">%s|%d|%s.%d|\tstrand\t%s\tgene_id\t%s\tgene_id\t%s\t", s->chrom, s->start, s->gene_name, s->idx, s->strand, s->gene_id, s->gene_name);
		s->transcript_id = str_arr_uniq(s->transcript_id, &(s->transcript_num));
		s->tss_id = str_arr_uniq(s->tss_id, &(s->tss_num));
		fprintf(fp, "transcript_id\t");
		if(s->transcript_num<1){
			fprintf(fp, "|\t", s->transcript_id[i]);
		}else{
			for(i=0; i<s->transcript_num; i++) fprintf(fp, "%s|", s->transcript_id[i]);			
		}
		fprintf(fp, "\t");
		fprintf(fp, "tss_id\t");		
		if(s->tss_num<1){
			fprintf(fp, "|\t", s->tss_id[i]);
		}else{			
			for(i=0; i<s->tss_num; i++) fprintf(fp, "%s|", s->tss_id[i]);
		}
		fprintf(fp, "\n");
		fprintf(fp, "%s\n", s->seq);
	}
	fclose(fp);
	return 0;
}

int name2fasta_usage(){
	fprintf(stderr, "\n");
			fprintf(stderr, "Usage:   tfc name2fasta [options] <gname.txt> <gencode.gtf> <in.fa.gz> <exons.fa> \n\n");
			fprintf(stderr, "Details: name2fasta is to extract genomic sequence of gene candiates\n\n");
			fprintf(stderr, "Options: -g          'exon' or 'transcript' \n\n");
			fprintf(stderr, "Inputs:  .txt        plain txt file contains names of gene candiates e.g. [genes.txt]\n");
			fprintf(stderr, "         .gtf        gft file that contains gene annotation [genes.gtf]\n");
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
	
	if (optind + 4 > argc) return name2fasta_usage();
	gene_name = argv[optind];
	gff_name = argv[optind+1];
	iname = argv[optind+2];
	oname = argv[optind+3];

	if(strcmp(genr, "transcript")!=0 && strcmp(genr, "exon")!=0){
		fprintf(stderr, "-g unrecognized gener 'exon' or 'transcript'\n");
		return -1;
	}
		
	fasta_t *GENO_HT = NULL;
	fasta_t *EXON_HT = NULL;
	
	fprintf(stderr, "[%s] loading reference genome sequences ... \n",__func__);
	if((GENO_HT = fasta_read(iname)) == NULL) die("[%s] can't load reference genome %s", __func__, iname);	
	
	fprintf(stderr, "[%s] extracting targeted gene sequences ... \n",__func__);
	if((EXON_HT = extract_exon_seq(gene_name, gff_name, GENO_HT))==NULL) die("[%s] can't extract exon sequences of %s", __func__, gene_name);
	
	fprintf(stderr, "[%s] writing down sequences ... \n",__func__);
	if((fasta_write_exon(EXON_HT, oname))!=0) die("[%s] can't write down to %s", __func__, oname);

	fprintf(stderr, "[%s] cleaning up ... \n", __func__);	
    
	if(EXON_HT)   fasta_destroy(&EXON_HT);
	if(GENO_HT)   fasta_destroy(&GENO_HT);    
	return 0;
}