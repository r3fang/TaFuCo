/*--------------------------------------------------------------------*/
/* main.c                                                             */
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
#include "kmer_uthash.h"
#include "BAG_uthash.h"
#include "fasta_uthash.h"
#include "utils.h"
#include "alignment.h"
#include "junction.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.30-r15"
#endif

#ifndef EXON_FILE
#define EXON_FILE "sample_data/genes.bed"
#endif


/* error code */
#define PR_ERR_NONE		     	 		0 		    // error code
#define MAX_ALLOWED_K                   50
#define EXON_HALF_FLANK                 0

/*--------------------------------------------------------------------*/
/*Global paramters.*/
static struct fasta_uthash *GENO_HT     = NULL;  // reference genome "hg19"
static struct kmer_uthash  *KMER_HT     = NULL;  // kmer hash table
static struct fasta_uthash *EXON_HT     = NULL;  // extracted exon sequences
static struct BAG_uthash   *BAGR_HT     = NULL;  // Breakend Associated Graph
static        junction_t   *JUNC_HT     = NULL;  // Identified Junction sites
static   solution_pair_t   *SOLU_HT     = NULL;  // Identified Junction sites

static char* concat_exons(char*, struct fasta_uthash *, struct kmer_uthash *, int, char *, char* ,  char **, char **, int *);
static int find_junction_one_edge(struct BAG_uthash*, struct fasta_uthash *, opt_t *, junction_t **);
static junction_t *junction_construct(struct BAG_uthash *, struct fasta_uthash *, opt_t *);
static struct fasta_uthash *extract_exon_seq(char*, struct fasta_uthash *);
static struct kmer_uthash  *kmer_uthash_construct(struct fasta_uthash *, int);
static struct BAG_uthash   *BAG_uthash_construct(struct kmer_uthash *, char*, char*, int, int);
static int junction_rediscover_unit(junction_t *, opt_t *, solution_pair_t **);

static inline int
find_junction_one_edge(struct BAG_uthash *eg, struct fasta_uthash *fasta_u, opt_t *opt, junction_t **ret){
	int _k = opt->k;
	int num;
	char* gname1 = strsplit(eg->edge, '_', &num)[0];
	char* gname2 = strsplit(eg->edge, '_', &num)[1];
	register int i, j;
	char* idx;
	register struct kmer_uthash *s_kmer;
	junction_t *m, *n;
	register solution_t *a, *b;
	register char* _read1, *_read2, *str2;
	char *chrom1, *chrom2;
	chrom2 = chrom1 = NULL;
	int start1, start2;
	int junction;
	char *ename1, *ename2;
	int strlen2;
	for(i=0; i<eg->weight; i++){
		_read1 = strsplit(eg->evidence[i], '_', &num)[0];
		_read2 = strsplit(eg->evidence[i], '_', &num)[1];	
		if((str2 =  concat_exons(_read1, fasta_u, KMER_HT, _k, gname1, gname2, &ename1, &ename2, &junction))==NULL) continue;
		
		a = align(_read1, str2, junction, opt);
		b = align(_read2, str2, junction, opt);
		if(a->jump == true && a->prob >= opt->min_align_score){
			idx = idx2str(concat(concat(ename1, "."), ename2), a->jump_start, a->jump_end);
			HASH_FIND_STR(*ret, idx, m);
			if(m==NULL){ // this junction not in ret
				m = junction_init(opt->seed_len);				
				m->idx    = idx;
				m->exon1  = ename1;
				m->exon2  = ename2;				
				m->hits  = 1;
				m->likehood = 10*log(a->prob); 				
				// 20bp junction flanking sequence 
				memcpy( m->s, &str2[a->jump_start-opt->seed_len/2-1], opt->seed_len/2);
				memcpy( &m->s[opt->seed_len/2], &str2[a->jump_end], opt->seed_len/2);
				// junction flanking sequence 
				strlen2 = a->jump_start + strlen(str2)-a->jump_end+1;				
				m->transcript = mycalloc(strlen2, char);
				memset(m->transcript, '\0', strlen2);
				memcpy( m->transcript, str2, a->jump_start);
				memcpy( &m->transcript[a->jump_start], &str2[a->jump_end+1], strlen(str2)-a->jump_end);				
				HASH_ADD_STR(*ret, idx, m);
			}else{
				m->hits ++;
				m->likehood += 10*log(a->prob); 
			}
		}
		if(b->jump == true && b->prob >= opt->min_align_score){
			idx = idx2str(concat(concat(ename1, "."), ename2), b->jump_start, b->jump_end);
			HASH_FIND_STR(*ret, idx, m);
			if(m==NULL){ // this junction not in ret
				m = junction_init(opt->seed_len);				
				m->idx   = idx;
				m->exon1  = ename1;
				m->exon2  = ename2;				
				m->hits  = 1;
				m->likehood = 10*log(b->prob); 				
				memcpy( m->s, &str2[b->jump_start-opt->seed_len/2-1], opt->seed_len/2);
				memcpy( &m->s[opt->seed_len/2], &str2[b->jump_end], opt->seed_len/2);
				strlen2 = b->jump_start + strlen(str2)-b->jump_end+1;
				m->transcript = mycalloc(strlen2, char);
				memset(m->transcript, '\0', strlen2);
				memcpy( m->transcript, str2, b->jump_start);
				memcpy( &m->transcript[b->jump_start], &str2[b->jump_end], strlen(str2)-b->jump_end);				
				HASH_ADD_STR(*ret, idx, m);
			}else{
				m->hits ++;
				m->likehood += 10*log(b->prob); 
			}
		}	
	}
	// delete those junctions with hits < MIN_HITS
	HASH_ITER(hh, *ret, m, n){
		if(m != NULL){
			m->likehood = m->likehood/m->hits;
			if(m->hits < opt->min_hits){
				HASH_DEL(*ret,m);
				free(m);
			}
		}
	}
	if(ename1)         free(ename1);
	if(ename2)         free(ename2);
	if(str2)           free(str2);
	if(_read1)         free(_read1);
	if(_read2)         free(_read2);
	if(a)              solution_destory(a);
	if(b)              solution_destory(b);
	if(idx)            free(idx);
	if(gname1)         free(gname1);     
	if(gname2)         free(gname2);
	return 0;
}
/*
 * generate junction sites from breakend associated graph 
 */
static junction_t 
*junction_construct(struct BAG_uthash *tb, struct fasta_uthash *fa, opt_t *opt){
	if(tb == NULL || fa == NULL || opt==NULL) die("[%s] input error", __func__);	
	struct BAG_uthash *edge, *tmp_bag;
	register int i;
	junction_t *cur_junction, *tmp_junction, *ret = NULL;
	HASH_ITER(hh, tb, edge, tmp_bag){ // iterate every edge
		find_junction_one_edge(edge, fa, opt, &ret);
	}
	return ret;
}

static char 
*concat_exons(char* _read, struct fasta_uthash *fa_ht, struct kmer_uthash *kmer_ht, int _k, char *gname1, char* gname2, char** ename1, char** ename2, int *junction){
	if(_read == NULL || fa_ht == NULL || kmer_ht==NULL || gname1==NULL || gname2==NULL) return NULL;
	char *str1, *str2;
	*ename1 = *ename2 = str1 = str2 = NULL;
	int order_gene1, order_gene2;
	order_gene1 = order_gene2 = 0;
	int i=0;
	int num_tmp;
	char buff[_k];
	char* gname_cur;
	str_ctr *s_ctr, *tmp_ctr, *exons=NULL;
	struct fasta_uthash *fa_tmp;
	find_all_exons(&exons, kmer_ht, _read, _k);
	if(exons==NULL) return NULL; // no exon found
	HASH_ITER(hh, exons, s_ctr, tmp_ctr) { 
		if(s_ctr->SIZE >= 2){ //denoise
			gname_cur = strsplit(s_ctr->KEY, '.', &num_tmp)[0];
			
			if(strcmp(gname_cur, gname1)==0){
				fa_tmp = find_fasta(fa_ht, s_ctr->KEY);
				if(str1 == NULL){
					str1 = strdup(fa_tmp->seq);
					order_gene1 = i++;
					*ename1 = s_ctr->KEY;
				}else{
					*ename1 = s_ctr->KEY;
					str1 = concat(str1, fa_tmp->seq);
				}				
			}
			if(strcmp(gname_cur, gname2)==0){
				fa_tmp = find_fasta(fa_ht, s_ctr->KEY);
				if(str2 == NULL){
					str2 = strdup(fa_tmp->seq);
					*ename2 = s_ctr->KEY;
					order_gene2 = i++;
				}else{
					str2 = concat(str2, fa_tmp->seq);
				}				
			}
		}
	}
	if(gname_cur) free(gname_cur);
	if(exons)         str_ctr_destory(&exons);
	if(str1 == NULL || str2 == NULL) return NULL;
	char *ret = (order_gene1 <= order_gene2) ? concat(str1, str2) : concat(str2, str1);	
	*junction = (order_gene1 <= order_gene2) ? strlen(str1) : strlen(str2);	
	if(order_gene1 > order_gene2){ // switch ename1 and ename2
		char* tmp = strdup(*ename1);
		*ename1 = strdup(*ename2);
		*ename2 = strdup(tmp);
		if(tmp) free(tmp);
 	}
	return ret;
}

static struct fasta_uthash 
*extract_exon_seq(char* fname, struct fasta_uthash *HG19_HT){
	if(fname == NULL || HG19_HT == NULL) return NULL;
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
	FILE *fp0 = fopen(fname, "r");
	if(fp0==NULL) die("[%s] can't open %s", __func__, fname); 
	while ((read = getline(&line, &len, fp0)) != -1) {
		if((fields = strsplit(line, 0, &num))==NULL) continue; // get rid of \n 
		str_ctr_add(&gene_name_ctr, fields[0]);		
	}
	fclose(fp0);
	
	FILE *fp = fopen(EXON_FILE, "r");
	if(fp==NULL) die("[%s] can't open %s", __func__, EXON_FILE);
	while ((read = getline(&line, &len, fp)) != -1) {
		// get information of exons
		if((fields = strsplit(line, 0, &num))==NULL) continue;
		if(num < 7) continue;
		if((chrom = fields[0])==NULL) continue;
		if((category = fields[2])==NULL) continue;
		if((start = atoi(fields[3]))<0) continue;
		if((end = atoi(fields[4]))<0) continue;
		if((strand = fields[5])==NULL) continue;
		if((gname = fields[6])==NULL) continue;
		if(strcmp(category, "exon")!=0) continue;
		if((find_str_ctr(gene_name_ctr, gname)) == NULL) continue; // only for targetted genes
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
	fclose(fp);
	if(gname) free(gname);
	if (line) free(line);
	if(strand)   free(strand);
	if(category) free(category);
	return ret_fasta;
}

static struct kmer_uthash 
*kmer_uthash_construct(struct fasta_uthash *tb, int k){
	if(tb == NULL || k < 0 || k > MAX_ALLOWED_K) return NULL;
	register char *kmer;
	char *name = NULL;	
	char *seq = NULL;	
	register int i, j;
	struct kmer_uthash  *s_kmer, *tmp_kmer, *ret = NULL;
	struct fasta_uthash *s_fasta, *tmp_fasta;
	HASH_ITER(hh, tb, s_fasta, tmp_fasta) {
		seq = strToUpper(s_fasta->seq);
		name = s_fasta->name;
		if(seq == NULL || name == NULL || strlen(seq) <= k){
			continue;
		}
		for(i=0; i < strlen(seq)-k+1; i++){
			kmer = mycalloc(k+1, char);
			memset(kmer, '\0', k+1);
			strncpy(kmer, strToUpper(s_fasta->seq)+i, k);
			kmer_uthash_insert(&ret, kmer, name); 
		}
	}
	if(kmer) free(kmer);
	if(seq)  free(seq);
	if(name)  free(name);
	kmer_uthash_uniq(&ret);
	return ret;
}

/* 
 * Construct breakend associated graph.             
 * fq_file1 - fastq file name for R1
 * fq_file2 - fastq file name for R2
 * _k       - kmer length for hash table
 * cutoff   - min number of kmer matches between a read again a gene.
 * bag      - BAG_uthash object (breakend associated graph)
 */
static struct BAG_uthash
*BAG_uthash_construct(struct kmer_uthash *kmer_uthash, char* fq1, char* fq2, int _min_match, int _k){
	if(kmer_uthash == NULL || fq1 == NULL || fq2 == NULL) return NULL;
	struct BAG_uthash *bag = NULL;
	int error;
	gzFile fp1, fp2;
	register int l1, l2;
	register kseq_t *seq1, *seq2;
	register char *_read1, *_read2;
	register char *edge_name;
	_read1 = _read2 = edge_name = NULL;
	
	if((fp1 = gzopen(fq1, "r"))==NULL) die("[%s] fail to read fastq files", __func__);
	if((fp2 = gzopen(fq2, "r"))==NULL) die("[%s] fail to read fastq files", __func__);	
	if((seq1 = kseq_init(fp1))==NULL)  die("[%s] fail to read fastq files", __func__);
	if((seq2 = kseq_init(fp2))==NULL)  die("[%s] fail to read fastq files", __func__);

	while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >= 0 ) {
		_read1 = rev_com(seq1->seq.s); // reverse complement of read1
		_read2 = seq2->seq.s;		
		if(_read1 == NULL || _read2 == NULL) continue;
		if(strcmp(seq1->name.s, seq2->name.s) != 0) die("[%s] read pair not matched", __func__);		
		if(strlen(_read1) < _k || strlen(_read2) < _k){continue;}
		str_ctr* gene_counter = NULL;
		str_ctr *s, *tmp;
		find_all_genes(&gene_counter, kmer_uthash, _read1, _k);
		find_all_genes(&gene_counter, kmer_uthash, _read2, _k);
		
		HASH_SORT(gene_counter, str_ctr_sort);
		unsigned int num = HASH_COUNT(gene_counter);
		
		char** hits = mycalloc(num, char*);
		int i=0; if(num > 1){
			HASH_ITER(hh, gene_counter, s, tmp) { 
				if(s->SIZE >= _min_match){hits[i++] = strdup(s->KEY);}
			}			
		}
		int m, n; for(m=0; m < i; m++){for(n=m+1; n < i; n++){
				int rc = strcmp(hits[m], hits[n]);
				if(rc<0)  edge_name = concat(concat(hits[m], "_"), hits[n]);
				if(rc>0)  edge_name = concat(concat(hits[n], "_"), hits[m]);
				if(rc==0) edge_name = NULL;
				if(edge_name!=NULL){
					if(BAG_uthash_add(&bag, edge_name, concat(concat(_read1, "_"), _read2)) != 0) die("BAG_uthash_add fails\n");							
				}
		}}
		if(hits)		 free(hits);
		if(gene_counter) free(gene_counter);
	}
	kseq_destroy(seq1);
	kseq_destroy(seq2);	
	gzclose(fp1);
	gzclose(fp2);
	return bag;
}

static solution_pair_t *junction_rediscover(junction_t *junc, opt_t *opt){
	if(junc==NULL || opt==NULL) return NULL;
	solution_pair_t *res = NULL;
	junction_t *cur_junction, *tmp_junction;
	HASH_ITER(hh, junc, cur_junction, tmp_junction) {
		fprintf(stderr, "junction between %s and %s\n", cur_junction->exon1, cur_junction->exon2);
		if((junction_rediscover_unit(cur_junction, opt, &res))!=0) return NULL;
	}
	return res;
}

static int junction_rediscover_unit(junction_t *junc, opt_t *opt, solution_pair_t **sol_pair){
	if(junc==NULL || opt==NULL) return -1;
	int mismatch = opt->max_mismatch;
	gzFile fp1, fp2;
	int l1, l2;
	kseq_t *seq1, *seq2;
	register char *_read1, *_read2;
	solution_t *sol1, *sol2;
	sol1 = sol2 = NULL;
	solution_pair_t *s_sp, *tmp_sp;
	if((fp1  = gzopen(opt->fq1, "r")) == NULL)   die("[%s] fail to read fastq files\n",  __func__);
	if((fp2  = gzopen(opt->fq2, "r")) == NULL)   die("[%s] fail to read fastq files\n",  __func__);	
	if((seq1 = kseq_init(fp1))   == NULL)   die("[%s] fail to read fastq files\n",  __func__);
	if((seq2 = kseq_init(fp2))   == NULL)   die("[%s] fail to read fastq files\n",  __func__);	
	while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >= 0 ) {
		_read1 = rev_com(seq1->seq.s); // reverse complement of read1
		_read2 = seq2->seq.s;		
		if(_read1 == NULL || _read2 == NULL) die("[%s] fail to get _read1 and _read2\n", __func__);
		if(strcmp(seq1->name.s, seq2->name.s) != 0) die("[%s] read pair not matched\n", __func__);
		if((min_mismatch(_read1, junc->s)) <= mismatch || (min_mismatch(_read2, junc->s)) <= mismatch ){
			if((sol1 = align_exon_jump(_read1, junc->transcript, junc->S1, junc->S2, junc->S1_num, junc->S2_num, opt))==NULL) continue;
			if((sol2 = align_exon_jump(_read2, junc->transcript, junc->S1, junc->S2, junc->S1_num, junc->S2_num, opt))==NULL) continue;
			s_sp = find_solution_pair(*sol_pair, seq1->name.s);
			if(s_sp!=NULL){ // if exists
				if(s_sp->prob < sol1->prob*sol2->prob){
					s_sp->r1 = sol1; s_sp->r2=sol2; 
					s_sp->prob = sol1->prob*sol2->prob; 
					s_sp->junc_name = junc->idx;
				}				
			}else{
				if(sol1->prob >= opt->min_align_score && sol2->prob >= opt->min_align_score){
					s_sp = solution_pair_init();
					s_sp->idx = strdup(seq1->name.s); 
					s_sp->junc_name = junc->idx;
					s_sp->r1 = sol1;			
					s_sp->r2 = sol2;
					s_sp->prob = sol1->prob*sol2->prob;
					HASH_ADD_STR(*sol_pair, idx, s_sp);
				}
			}
		}
	}
	if(sol1)     solution_destory(sol1);
	if(sol2)     solution_destory(sol2);
	if(seq1)     kseq_destroy(seq1);
	if(seq2)     kseq_destroy(seq2);	
	if(fp1)      gzclose(fp1);
	if(fp2)      gzclose(fp2);
	return 0;	
}
static junction_t 
*transcript_construct(junction_t *junc_ht, struct fasta_uthash *exon_ht){
	if(junc_ht==NULL || exon_ht==NULL) return NULL;
	junction_t *cur_junction, *tmp_junction;
	char* gname1, *gname2, *ename1, *ename2;
	gname1 = gname2 = ename1 = ename2 = NULL;
	int enum1, enum2;
	char enum1_str[10], enum2_str[10];
	int num1, num2;
	char** fields1, **fields2;
	struct fasta_uthash *exon1_fa,  *exon2_fa;
	exon1_fa = exon2_fa = NULL;
	char *exon1_seq, *exon2_seq;
	int i, j;
	int *S1, *S2;
	int str1_l, str2_l, str3_l;
	HASH_ITER(hh, junc_ht, cur_junction, tmp_junction) {
		exon1_seq = exon2_seq = NULL;
		S1 = S2 = NULL;
		cur_junction->S1_num = cur_junction->S2_num = 0;
		cur_junction->S1 = cur_junction->S2 = NULL;
		
		fields1 = strsplit(cur_junction->exon1, '.', &num1);
		fields2 = strsplit(cur_junction->exon2, '.', &num2);
		if(num1 != 2 || num2 != 2) continue;
		gname1 = fields1[0]; enum1 = atoi(fields1[1]);
		gname2 = fields2[0]; enum2 = atoi(fields2[1]);
		S1 = mycalloc(enum1+1, int);
		for(i=1; i<enum1; i++){
			sprintf(enum1_str, "%d", i);
			ename1 = concat(concat(gname1, "."), enum1_str);
			if((exon1_fa = find_fasta(exon_ht, ename1))!=NULL){
				exon1_seq = concat(exon1_seq, exon1_fa->seq);
				S1[i-1] = strlen(exon1_seq);
				cur_junction->S1_num = i;
			}
		}
		str1_l = (exon1_seq == NULL) ? 0 : strlen(exon1_seq);
		str2_l = (cur_junction->transcript == NULL) ? 0 : strlen(cur_junction->transcript);
		S2 = mycalloc(100, int); memset(S2, INFINITY, 100);
		S2[0] = str1_l + str2_l;
		j = 1;for(i=(enum2+1); i<INFINITY; i++){
			sprintf(enum2_str, "%d", i);
			ename2 = concat(concat(gname2, "."), enum2_str);
			if((exon2_fa = find_fasta(exon_ht, ename2))==NULL) break;
			exon2_seq = concat(exon2_seq, exon2_fa->seq);
			cur_junction->S2_num = j;
			S2[j++] = str1_l + str2_l + strlen(exon2_seq);
		}
		cur_junction->transcript = concat(concat(exon1_seq, cur_junction->transcript), exon2_seq);
		cur_junction->S1 = S1;
		cur_junction->S2 = S2;
	}	
	if(fields1)       free(fields1);
	if(fields2)       free(fields2);
	if(gname1)        free(gname1);
	if(gname2)        free(gname2);
	return junc_ht;
}

static int tfc_usage(opt_t *opt){
	fprintf(stderr, "\n");
			fprintf(stderr, "Usage:   tfc [options] <target.bed> <genome.fa> <R1.fq> <R2.fq>\n\n");
			fprintf(stderr, "Options: -k INT   kmer length [%d]\n", opt->k);
			fprintf(stderr, "         -n INT   min number kmer matches [%d]\n", opt->min_match);
			fprintf(stderr, "         -w INT   min weight for an edge [%d]\n", opt->min_weight);
			fprintf(stderr, "         -m INT   match score [%d]\n", opt->match);
			fprintf(stderr, "         -u INT   penality for mismatch [%d]\n", opt->mismatch);
			fprintf(stderr, "         -o INT   penality for gap open [%d]\n", opt->gap);
			fprintf(stderr, "         -e INT   penality for gap extension [%d]\n", opt->extension);
			fprintf(stderr, "         -j INT   penality for jump between genes [%d]\n", opt->jump_gene);
			fprintf(stderr, "         -s INT   penality for jump between exons [%d]\n", opt->jump_exon);
			fprintf(stderr, "         -h INT   min number of hits for a junction to be called [%d]\n", opt->min_hits);					
			fprintf(stderr, "         -l INT   length of exact match seed for junction rediscovery [%d]\n", opt->seed_len);					
			fprintf(stderr, "         -x INT   max number of mismatches of seed exact match for junction rediscovery [%d]\n", opt->max_mismatch);					
			fprintf(stderr, "         -a FLOAT min alignment score [%.2f]\n", opt->min_align_score);
			fprintf(stderr, "\n");
			return 1;
}
/*--------------------------------------------------------------------*/
/* main function. */
int main(int argc, char *argv[]) {
	opt_t *opt = init_opt(); // initlize options with default settings
	int c, i;
	srand48(11);
	junction_t *junc_ht;
	while ((c = getopt(argc, argv, "m:w:k:n:u:o:e:g:s")) >= 0) {
				switch (c) {
				case 'n': opt->min_match = atoi(optarg); break;
				case 'w': opt->min_weight = atoi(optarg); break;
				case 'k': opt->k = atoi(optarg); break;
				case 'm': opt->match = atoi(optarg); break;
				case 'u': opt->mismatch = atoi(optarg); break;
				case 'o': opt->gap = atoi(optarg); break;
				case 'e': opt->extension = atoi(optarg); break;
				case 'j': opt->jump_gene = atoi(optarg); break;
				case 's': opt->jump_exon = atoi(optarg); break;
				case 'h': opt->min_hits = atoi(optarg); break;
				case 'l': opt->seed_len = atoi(optarg); break;
				case 'x': opt->max_mismatch = atoi(optarg); break;
				case 'a': opt->min_align_score = atof(optarg); break;
				default: return 1;
		}
	}
	if (optind + 4 > argc) return tfc_usage(opt);
	opt->bed = argv[optind];
	opt->fa = argv[optind+1];
	opt->fq1 = argv[optind+2];
	opt->fq2 = argv[optind+3];
	
	//fprintf(stderr, "[%s] loading reference genome sequences ... \n",__func__);
	//if((GENO_HT = fasta_uthash_load(opt->fa)) == NULL) die("[%s] can't load reference genome %s", __func__, opt->fa);	
	//
	//fprintf(stderr, "[%s] extracting targeted gene sequences ... \n",__func__);
	//if((EXON_HT = extract_exon_seq(opt->bed, GENO_HT))==NULL) die("[%s] can't extract exon sequences of %s", __func__, opt->bed);
		
	fprintf(stderr, "[%s] loading reference exon sequences ... \n",__func__);
	if((EXON_HT = fasta_uthash_load(opt->fa)) == NULL) die("[%s] can't load reference genome %s", __func__, opt->fa);	

	fprintf(stderr, "[%s] indexing exon sequneces ... \n",__func__);
	if((KMER_HT = kmer_uthash_construct(EXON_HT, opt->k))==NULL) die("[%s] can't index exon sequences", __func__); 	

	fprintf(stderr, "[%s] constructing graph ... \n", __func__);
	if((BAGR_HT = BAG_uthash_construct(KMER_HT, opt->fq1, opt->fq2, opt->min_match, opt->k)) == NULL)	die("[%s] can't construct BAG graph", __func__); 	

	fprintf(stderr, "[%s] identifying junction sites from graph ... \n", __func__);
	if((JUNC_HT = junction_construct(BAGR_HT, EXON_HT, opt))==NULL) die("[%s] can't identify junctions", __func__);

	fprintf(stderr, "[%s] construct trnascript ... \n", __func__);    
	if((JUNC_HT = transcript_construct(JUNC_HT, EXON_HT))==NULL) die("[%s] can't construct transcript", __func__);
	
	fprintf(stderr, "[%s] rediscover junctions from the reads ... \n", __func__);    	
	if((SOLU_HT = junction_rediscover(JUNC_HT, opt))==NULL) die("[%s] can't rediscover any junction", __func__);;

	solution_pair_t *cur_sop, *tmp_sop;
	fprintf(stderr, "[%s] scoring junction ... \n", __func__);	
	HASH_ITER(hh, SOLU_HT, cur_sop, tmp_sop) {
		printf("%s\t%f\n", cur_sop->junc_name, cur_sop->prob);
	}	
	
	fprintf(stderr, "[%s] cleaning up ... \n", __func__);	
	if(SOLU_HT)  solution_pair_destory(&SOLU_HT);
	if(JUNC_HT)       junction_destory(&JUNC_HT);
	if(BAGR_HT)     BAG_uthash_destroy(&BAGR_HT);
	if(KMER_HT)    kmer_uthash_destroy(&KMER_HT);
	if(EXON_HT)   fasta_uthash_destroy(&EXON_HT);
	//if(GENO_HT)   fasta_uthash_destroy(&GENO_HT);
    
	fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
	fprintf(stderr, "[%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n");
	return 0;
}
