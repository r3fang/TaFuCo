/*--------------------------------------------------------------------*/
/* predict.c                                                          */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Predict Gene Fusion by given fastq files.                          */
/*--------------------------------------------------------------------*/
#include "predict.h"

static char *concat_exons(char* _read, struct fasta_uthash *fa_ht, struct kmer_uthash *kmer_ht, int _k, char *gname1, char* gname2, char** ename1, char** ename2, int *junction);
static int find_junction_one_edge(bag_t *eg, struct fasta_uthash *fasta_u, opt_t *opt, junction_t **ret);
static int align_to_transcript_unit(junction_t *junc, opt_t *opt, solution_pair_t **sol_pair);

/* 
 * Find all genes uniquely matched with kmers on _read.          
 * hash     - a hash table count number of matches between _read and every gene
 * _read    - inqury read
 * _k       - kmer length
 */
static inline int
find_all_genes(str_ctr **hash, struct kmer_uthash *KMER_HT, char* _read, int _k){
/*--------------------------------------------------------------------*/
	/* check parameters */
	if(_read == NULL || _k < 0) die("find_all_MEKMs: parameter error\n");
/*--------------------------------------------------------------------*/
	/* declare vaiables */
	str_ctr *s;
	int _read_pos = 0;
	int num;
	char* gene = NULL;
	register struct kmer_uthash *s_kmer = NULL; 
	char buff[_k];
/*--------------------------------------------------------------------*/
	while(_read_pos<(strlen(_read)-_k+1)){
		/* copy a kmer of string */
		strncpy(buff, _read + _read_pos, _k); buff[_k] = '\0';	
		if(strlen(buff) != _k) die("find_next_match: buff strncpy fails\n");
		/*------------------------------------------------------------*/
		if((s_kmer=find_kmer(KMER_HT, buff)) == NULL){_read_pos++; continue;} // kmer not in table but not an error
		if(s_kmer->count == 1){ // only count the uniq match 
			//gene = strdup(s_kmer->seq_names[0]);
			gene = strsplit(s_kmer->seq_names[0], '.', &num)[0];
			if(gene == NULL) die("find_next_match: get_exon_name fails\n");
			if(str_ctr_add(hash, gene) != 0) die("find_all_MEKMs: str_ctr_add fails\n");
		}
		_read_pos++;
	}
	return 0;
}


/* 
 * Find all exons uniquely matched with kmers on _read.          
 * hash     - a hash table count number of matches between _read and every gene
 * _read    - inqury read
 * _k       - kmer length
 */
static inline int
find_all_exons(str_ctr **hash, struct kmer_uthash *KMER_HT, char* _read, int _k){
/*--------------------------------------------------------------------*/
	/* check parameters */
	if(_read == NULL || _k < 0) die("find_all_MEKMs: parameter error\n");
/*--------------------------------------------------------------------*/
	/* declare vaiables */
	str_ctr *s;
	int _read_pos = 0;
	char* exon = NULL;
	register struct kmer_uthash *s_kmer = NULL; 
	char buff[_k];
/*--------------------------------------------------------------------*/
	while(_read_pos<(strlen(_read)-_k+1)){
		/* copy a kmer of string */
		strncpy(buff, _read + _read_pos, _k); buff[_k] = '\0';	
		if(strlen(buff) != _k) die("find_next_match: buff strncpy fails\n");
		/*------------------------------------------------------------*/
		if((s_kmer=find_kmer(KMER_HT, buff)) == NULL){_read_pos++; continue;} // kmer not in table but not an error
		if(s_kmer->count == 1){ // only count the uniq match 
			exon = strdup(s_kmer->seq_names[0]);
			if(exon == NULL) die("find_next_match: get_exon_name fails\n");
			if(str_ctr_add(hash, exon) != 0) die("find_all_MEKMs: str_ctr_add fails\n");
		}
		_read_pos++;
	}
	return 0;
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
static bag_t
*bag_construct(struct kmer_uthash *kmer_uthash, char* fq1, char* fq2, int min_kmer_matches, int min_edge_weight, int _k){
	if(kmer_uthash == NULL || fq1 == NULL || fq2 == NULL) return NULL;
	bag_t *bag = NULL;
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
		//if(strcmp(seq1->name.s, seq2->name.s) != 0) die("[%s] read pair not matched", __func__);		
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
				if(s->SIZE >= min_kmer_matches){hits[i++] = strdup(s->KEY);}
			}			
		}
		int m, n; for(m=0; m < i; m++){for(n=m+1; n < i; n++){
				int rc = strcmp(hits[m], hits[n]);
				if(rc<0)  edge_name = concat(concat(hits[m], "_"), hits[n]);
				if(rc>0)  edge_name = concat(concat(hits[n], "_"), hits[m]);
				if(rc==0) edge_name = NULL;
				if(edge_name!=NULL){
					if(bag_add(&bag, edge_name, seq1->name.s, concat(concat(_read1, "_"), _read2)) != 0) die("BAG_uthash_add fails\n");							
				}
		}}
		if(hits)		 free(hits);
		if(gene_counter) free(gene_counter);
	}
	// remove the edges with weight < min_edge_weight	
	register bag_t *cur, *tmp;
	HASH_ITER(hh, bag, cur, tmp) {
		if(cur->weight < min_edge_weight) {HASH_DEL(bag, cur); free(cur);}
    }
	
	kseq_destroy(seq1);
	kseq_destroy(seq2);	
	gzclose(fp1);
	gzclose(fp2);
	return bag;
}


static int
edge_junction_gen(bag_t *eg, struct fasta_uthash *fasta_u, opt_t *opt, junction_t **ret){
	if(eg==NULL || fasta_u==NULL || opt==NULL) return -1;
	int _k = opt->k;
	int num;
	char *gname1, *gname2; 
	register char* _read1, *_read2, *str2;
	char *ename1, *ename2;
	char *chrom1, *chrom2;
	char* idx = NULL;
	chrom2 = chrom1 = NULL;
	ename1 = ename2 = NULL;
	_read1 = _read2 = str2 = NULL;
	gname1 = gname2 = NULL;
	gname1 = strsplit(eg->edge, '_', &num)[0];
	gname2 = strsplit(eg->edge, '_', &num)[1];
	register int i, j;
	register struct kmer_uthash *s_kmer;
	junction_t *m, *n;
	register solution_t *a, *b;
	int start1, start2;
	int junction;
	int strlen2;
	for(i=0; i<eg->weight; i++){
		_read1 = strsplit(eg->evidence[i], '_', &num)[0];
		_read2 = strsplit(eg->evidence[i], '_', &num)[1];	
		if((str2 =  concat_exons(_read1, fasta_u, KMER_HT, _k, gname1, gname2, &ename1, &ename2, &junction))==NULL) continue;
		
		a = align(_read1, str2, junction, opt->match, opt->mismatch, opt->gap, opt->extension, opt->jump_gene);
		b = align(_read2, str2, junction, opt->match, opt->mismatch, opt->gap, opt->extension, opt->jump_gene);

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
				// junction flanking sequence 
				memcpy( m->s, &str2[a->jump_start-opt->seed_len/2-1], opt->seed_len/2);
				memcpy( &m->s[opt->seed_len/2], &str2[a->jump_end], opt->seed_len/2);
				// junction flanking sequence 
				strlen2 = a->jump_start + strlen(str2)-a->jump_end+1;				
				m->transcript = mycalloc(strlen2, char);
				m->junc_pos = a->jump_start;				
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
				m->junc_pos = b->jump_start;
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
static bag_t 
*bag_junction_gen(bag_t *tb, struct fasta_uthash *fa, opt_t *opt){
	if(tb == NULL || fa == NULL || opt==NULL) return NULL;	
	bag_t *edge, *bag_cur, *res=NULL;
	register int i;
	junction_t *junc_cur = NULL;
	for(edge=tb; edge != NULL; edge=edge->hh.next) {
		edge_junction_gen(edge, fa, opt, &junc_cur);
		if((bag_cur=find_edge(res, edge->edge))==NULL){
			bag_cur = bag_init();
			bag_cur->edge = edge->edge;
			bag_cur->weight = edge->weight;
			bag_cur->read_names = edge->read_names;
			bag_cur->evidence = edge->evidence;
			bag_cur->junc = transcript_construct(junc_cur, fa);			
			HASH_ADD_STR(res, edge, bag_cur);								
		}
	}
	return res;
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

static solution_pair_t 
*align_edge_to_transcript(junction_t *junc_ht, bag_t *bag_ht, opt_t *opt){
	if(junc_ht==NULL || bag_ht==NULL || opt==NULL) return NULL;
    bag_t *edge;
   	junction_t *junc_cur;
	char *gene1, *gene2, *edge_name;
	gene1 = gene2 = NULL;
	int num, rc;
	solution_t *sol1, *sol2;
	solution_pair_t *sol_cur, *res = NULL;
	char *read1, *read2;
	read1 = read2 = NULL;
	int i, j;
	// iterate every junction
	for(junc_cur=junc_ht; junc_cur != NULL; junc_cur=junc_cur->hh.next){
		gene1 = strsplit(junc_cur->exon1, '.', &num)[0];
		gene2 = strsplit(junc_cur->exon2, '.', &num)[0];
		rc = strcmp(gene1, gene2);
		if(rc<0)  edge_name = concat(concat(gene1, "_"), gene2);
		if(rc>0)  edge_name = concat(concat(gene1, "_"), gene2);
		if(rc==0) continue;
		if((edge = find_edge(bag_ht, edge_name))==NULL) continue;
		for(i=0; i<edge->weight; i++){
			//printf("%s\t%d\n", edge->edge, i);
			read1 = strsplit(edge->evidence[i], '_', &num)[0];
			if(num != 2) continue;
			read2 = strsplit(edge->evidence[i], '_', &num)[1];
			if((sol1 = align_exon_jump(read1, junc_cur->transcript, junc_cur->S1, junc_cur->S2, junc_cur->S1_num, junc_cur->S2_num, opt->match, opt->mismatch, opt->gap, opt->extension, opt->jump_exon))==NULL) continue;
			if((sol2 = align_exon_jump(read2, junc_cur->transcript, junc_cur->S1, junc_cur->S2, junc_cur->S1_num, junc_cur->S2_num, opt->match, opt->mismatch, opt->gap, opt->extension, opt->jump_exon))==NULL) continue;
			if(sol1->prob < opt->min_align_score || sol2->prob < opt->min_align_score) continue;
			
			sol_cur = find_solution_pair(res, edge->read_names[i]);
			if(sol_cur!=NULL){ // if exists
				if(sol_cur->prob < sol1->prob*sol2->prob){
					sol_cur->r1        = sol1; 
					sol_cur->r2        = sol2; 
					sol_cur->prob      = (sol1->prob)*(sol2->prob); 
					sol_cur->junc_name = junc_cur->idx;
				}				
			}else{
					sol_cur = solution_pair_init();
					sol_cur->idx = strdup(edge->read_names[i]); 
					sol_cur->junc_name = junc_cur->idx;
					sol_cur->r1 = sol1;			
					sol_cur->r2 = sol2;
					sol_cur->prob = sol1->prob*sol2->prob;
					HASH_ADD_STR(res, idx, sol_cur);
			}	
		}
	}
	if(gene2)      free(gene2);
	if(gene1)      free(gene1);
	if(edge_name)  free(edge_name);
	if(read1)      free(read1);
	if(read2)      free(read2);	
	return res;
}

/*
 * Description:
 *------------
 * 1) find subset of pairs that contain 20bp junction string by at most 2 mismatches
 * 2) align those reads to constructed transcript returned by transcript_construct

 * Input: 
 *-------
 * junc_ht        - junction_t object: return by transcript_construct
 * opt            - opt_t object: contains all input parameters

 * Output: 
 *-------
 * solution_pair_t object that contains alignment results of all reads.
 */
static int align_reads_to_transcript(solution_pair_t **res, junction_t *junc, opt_t *opt){
	if(junc==NULL || opt==NULL) return -1;
	junction_t *cur_junction, *tmp_junction;
	HASH_ITER(hh, junc, cur_junction, tmp_junction) {
		fprintf(stderr, "[predict] fusion between %s and %s are being tested ... \n", cur_junction->exon1, cur_junction->exon2);
		if((align_to_transcript_unit(cur_junction, opt, res))!=0) return -1;
	}
	return 0;
}

/*
 * align reads to one junction transcript

 * junc      - one junction identified before
 * opt       - opt_t object
 * *sol_pair - solution_pair_t object that contains alignment solutions for all read pair agains junc

 */
static int align_to_transcript_unit(junction_t *junc, opt_t *opt, solution_pair_t **sol_pair){
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
	if((seq1 = kseq_init(fp1))   == NULL)        die("[%s] fail to read fastq files\n",  __func__);
	if((seq2 = kseq_init(fp2))   == NULL)        die("[%s] fail to read fastq files\n",  __func__);	
	
	while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >= 0 ) {
		_read1 = rev_com(seq1->seq.s); // reverse complement of read1
		_read2 = seq2->seq.s;		
		if(_read1 == NULL || _read2 == NULL) die("[%s] fail to get _read1 and _read2\n", __func__);
		if(strcmp(seq1->name.s, seq2->name.s) != 0) die("[%s] read pair not matched\n", __func__);
		if((min_mismatch(_read1, junc->s)) <= mismatch || (min_mismatch(_read2, junc->s)) <= mismatch ){	
			// alignment with jump state between exons 
			if((sol1 = align_exon_jump(_read1, junc->transcript, junc->S1, junc->S2, junc->S1_num, junc->S2_num, opt->match, opt->mismatch, opt->gap, opt->extension, opt->jump_exon))==NULL) continue;
			if((sol2 = align_exon_jump(_read2, junc->transcript, junc->S1, junc->S2, junc->S1_num, junc->S2_num, opt->match, opt->mismatch, opt->gap, opt->extension, opt->jump_exon))==NULL) continue;
			
			s_sp = find_solution_pair(*sol_pair, seq1->name.s);
			if(s_sp!=NULL){ // if exists
				if(s_sp->prob < sol1->prob*sol2->prob){
					s_sp->r1 = sol1; s_sp->r2=sol2; 
					s_sp->prob = (sol1->prob)*(sol2->prob); 
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
		S2 = mycalloc(100, int); memset(S2, UINT_MAX, 100);
		S2[0] = str1_l + str2_l;
		j = 1;for(i=(enum2+1); i<INFINITY; i++){
			sprintf(enum2_str, "%d", i);
			ename2 = concat(concat(gname2, "."), enum2_str);
			if((exon2_fa = find_fasta(exon_ht, ename2))==NULL) break;
			exon2_seq = concat(exon2_seq, exon2_fa->seq);
			cur_junction->S2_num = j;
			S2[j++] = str1_l + str2_l + strlen(exon2_seq);
		}
		cur_junction->junc_pos = str1_l + cur_junction->junc_pos;
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
/*
 * Description:
 *------------
 * revisit junction sites and score them based on alignment results
 
 * Input: 
 *-------
 * sol        - alignment results returned by align_to_transcript
 * junc       - previously identified junction
 * opt        - opt_t object: contains all input parameters

 * Output: 
 *-------
 * junction_t object that contains identified junctions with one more property -> transcript.
 */
static junction_t *junction_score(solution_pair_t *sol, junction_t *junc, double min_align_score, int junc_str_len){
	if(sol==NULL || junc==NULL) return NULL;
	junction_t *junc_cur1, *junc_cur2, *junc_res = NULL;
	solution_pair_t *sol_cur, *sol_tmp;
	double th = min_align_score * min_align_score;
	// iterate every alignment solution pair
	HASH_ITER(hh, sol, sol_cur, sol_tmp) {
		if(sol_cur->prob < th) continue; // filter read pairs with alignment identity smaller than th
		if((junc_cur1 = find_junction(junc, sol_cur->junc_name))==NULL) continue; // continue if junction not exists
		if((junc_cur2 = find_junction(junc_res, sol_cur->junc_name))==NULL){
			junc_cur2 = junction_init(junc_str_len);
			junc_cur2->idx = junc_cur1->idx;
			junc_cur2->exon1 = junc_cur1->exon1;
			junc_cur2->exon2 = junc_cur1->exon2;			
			junc_cur2->s = junc_cur1->s;			
			junc_cur2->junc_pos = junc_cur1->junc_pos;			
			junc_cur2->transcript = junc_cur1->transcript;			
			junc_cur2->hits = 1;	
			junc_cur2->likehood = 10*log(sol_cur->prob);
			HASH_ADD_STR(junc_res, idx, junc_cur2);
		}else{
			junc_cur2->hits++;	
			junc_cur2->likehood += 10*log(sol_cur->prob);			
		}
	}
	return junc_res;
}

//static int 
//junction_display(junction_t *junc, solution_pair_t *sol){
//	if(junc==NULL || sol==NULL){
//		fprintf(stderr, "[%s] junction is empty \n", __func__);
//		return -1;	
//	} 
//	junction_t *junc_cur, *junc_tmp;
//	int i, j;
//	solution_pair_t *sol_cur;
//	char** tmp = mycalloc(3, char*);
//	
//	HASH_ITER(hh, junc, junc_cur, junc_tmp) {
//		if(junc_cur->transcript == NULL) continue;
//		printf("fusion=%s-%s\thits=%zu\tjunction_pos=%d\tlikelihood=%.2f\n",junc_cur->exon1, junc_cur->exon2, junc_cur->hits, junc_cur->junc_pos, junc_cur->likehood);
//		for(i=0; i<junc_cur->junc_pos; i++) junc_cur->transcript[i] = tolower(junc_cur->transcript[i]);
//		printf_line(junc_cur->transcript, 50);
//		for(sol_cur=sol; sol_cur!=NULL; sol_cur=sol_cur->hh.next) {
//			if(sol_cur->junc_name == NULL) continue;
//			if(strcmp(sol_cur->junc_name, junc_cur->idx)==0){
//				printf(">%s\n", sol_cur->idx);
//				if(sol_cur->r1->s1 == NULL || sol_cur->r1->s2 == NULL) continue;
//				i = 1;
//				tmp[0] = tmp[1] = tmp[2] = mycalloc(51, char);
//				memset(tmp[0], '\0', 50);
//				memset(tmp[1], '\0', 50);
//				memset(tmp[2], '\0', 50);
//				while(i<=strlen(sol_cur->r1->s1)){
//					j = i%50;
//					tmp[0][j-1] = sol_cur->r1->s1[i-1]; 
//					tmp[2][j-1] = sol_cur->r1->s2[i-1]; 
//					if(j == 0){
//						printf("%s\n%s\n", tmp[0], tmp[2]);
//						memset(tmp[0], '\0', 50);
//						memset(tmp[1], '|',  50);
//						memset(tmp[2], '\0', 50);
//					}
//					i++;
//				}
//				printf("%s\n%s\n", tmp[0], tmp[2]);
//				printf("%s\t%d\n%s\n", sol_cur->r1->s2, sol_cur->r1->pos, sol_cur->r1->s1);
//				//printf("%s\t%d\n%s\n", sol_cur->r2->s2, sol_cur->r2->pos, sol_cur->r2->s1);		
//			}
//		 }		
//	}
//	return 0;
//}

static int pred_usage(opt_t *opt){
	fprintf(stderr, "\n");
			fprintf(stderr, "Usage:   tfc predict [options] <exon.fa> <R1.fq> <R2.fq>\n\n");
			fprintf(stderr, "Details: predict gene fusion from RNA-seq data\n\n");
			fprintf(stderr, "Options: -k INT    kmer length for indexing genome [%d]\n", opt->k);
			fprintf(stderr, "         -n INT    min number of uniq kmer matches for a gene-read match [%d]\n", opt->min_kmer_match);
			fprintf(stderr, "         -w INT    min weight for an edge in BAG [%d]\n", opt->min_edge_weight);
			fprintf(stderr, "         -m INT    score for match in alignment [%d]\n", opt->match);
			fprintf(stderr, "         -u INT    score for mismatch in alignment [%d]\n", opt->mismatch);
			fprintf(stderr, "         -o INT    penality for gap open [%d]\n", opt->gap);
			fprintf(stderr, "         -e INT    penality for gap extension [%d]\n", opt->extension);
			fprintf(stderr, "         -j INT    penality for jump between genes [%d]\n", opt->jump_gene);
			fprintf(stderr, "         -s INT    penality for jump between exons [%d]\n", opt->jump_exon);
			fprintf(stderr, "         -a FLOAT  min identity score for alignment [%.2f]\n", opt->min_align_score);
			fprintf(stderr, "         -h INT    min hits for a junction [%d]\n", opt->min_hits);					
			fprintf(stderr, "         -l INT    length for junction string [%d]\n", opt->seed_len);					
			fprintf(stderr, "         -x INT    max mismatches of junction string match [%d]\n", opt->max_mismatch);					
			fprintf(stderr, "\n");
			fprintf(stderr, "Inputs:  exon.fa   fasta file that contains exon sequences of targeted \n");
			fprintf(stderr, "                   genes with no flanking sequence which can be generated: \n");
			fprintf(stderr, "                   tfc name2fasta <genes.txt> <in.fa> <exon.fa> \n");
			fprintf(stderr, "         R1.fq     5'->3' end of pair-end sequencing reads\n");
			fprintf(stderr, "         R2.fq     the other end of pair-end sequencing reads\n");
			return 1;
}

/*--------------------------------------------------------------------*/
/* main function. */
int predict(int argc, char *argv[]) {
	opt_t *opt = opt_init(); // initlize options with default settings
	int c, i;
	//srand48(11);
	junction_t *junc_ht;
	while ((c = getopt(argc, argv, "m:w:k:n:u:o:e:g:s:h:l:x:a")) >= 0) {
				switch (c) {
				case 'n': opt->min_kmer_match = atoi(optarg); break;
				case 'w': opt->min_edge_weight = atoi(optarg); break;
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
	if (optind + 3 > argc) return pred_usage(opt);
	opt->fa  = argv[optind];    // exon.fa
	opt->fq1 = argv[optind+1];  // read1
	opt->fq2 = argv[optind+2];  // read2
	
	fprintf(stderr, "[%s] loading reference exon sequences ... \n",__func__);
	if((EXON_HT = fasta_uthash_load(opt->fa)) == NULL) die("[%s] can't load reference sequences %s", __func__, opt->fa);	

	fprintf(stderr, "[%s] indexing exon sequneces ... \n",__func__);
	if((KMER_HT = kmer_uthash_construct(EXON_HT, opt->k))==NULL) die("[%s] can't index exon sequences", __func__); 	

	fprintf(stderr, "[%s] constructing graph ... \n", __func__);
	if((BAGR_HT = bag_construct(KMER_HT, opt->fq1, opt->fq2, opt->min_kmer_match, opt->min_edge_weight, opt->k)) == NULL) return 0;

	fprintf(stderr, "[%s] deleting edges of low weight in the graph ... \n", __func__);
	if((BAGR_HT = bag_trim(BAGR_HT, opt->min_edge_weight))==NULL) return 0;

	fprintf(stderr, "[%s] deleting duplicate reads for every edge ... \n", __func__);
	if((BAGR_HT = bag_uniq(BAGR_HT))==NULL) return 0;
	
	fprintf(stderr, "[%s] identify fusion junctions and construct fused transcript for every fusion ... \n", __func__);
	if((BAGR_HT = bag_junction_gen(BAGR_HT, EXON_HT, opt))==NULL) return 0;
	
	//fprintf(stderr, "[%s] algin edges to constructed transcript ... \n", __func__);    
	//if((SOLU_HT = align_edge_to_transcript(BAGR_HT, opt))==NULL) die("[%s] can't rediscover any junction", __func__);

		
	//junction_t *s;
	//for(s=JUN0_HT; s != NULL; s=s->hh.next){
	//	printf("exon1=%s\texon2=%s\thits=%d\tstr=%s\n", s->exon1, s->exon2, s->hits, s->s);
	//}
	//{ // if junction string identified
	//	fprintf(stderr, "[%s] construct trnascript ... \n", __func__);    
	//	if((JUN0_HT = transcript_construct(JUN0_HT, EXON_HT))==NULL) die("[%s] can't construct transcript", __func__);		
	//	fprintf(stderr, "[%s] algin edges to constructed transcript ... \n", __func__);    
	//	if((SOLU_HT = align_edge_to_transcript(JUN0_HT, BAGR_HT, opt))==NULL) die("[%s] can't rediscover any junction", __func__);
	//}
	//fprintf(stderr, "[%s] align reads to fused transcript ... \n", __func__);    	
	//if((align_reads_to_transcript(&SOLU_HT, JUN0_HT, opt)) != 0) die("[%s] can't rediscover any junction", __func__);;

	//fprintf(stderr, "[%s] scoring junctions ... \n", __func__);    	
	//if((JUN1_HT = junction_score(SOLU_HT, JUN0_HT, opt->min_align_score, opt->seed_len+1))==NULL) die("[%s] can't rediscover any junction", __func__);;
	//
	//if((junction_display(JUN1_HT, SOLU_HT)) != 0) die("[%s] can't display junctions", __func__);
	
	fprintf(stderr, "[%s] cleaning up ... \n", __func__);	
	if(SOLU_HT)  solution_pair_destory(&SOLU_HT);
	if(BAGR_HT)            bag_distory(&BAGR_HT);
	if(KMER_HT)    kmer_uthash_destroy(&KMER_HT);
	if(EXON_HT)   fasta_uthash_destroy(&EXON_HT);
	return 0;
}

