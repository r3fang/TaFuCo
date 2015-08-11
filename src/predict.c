/*--------------------------------------------------------------------*/
/* predict.c                                                          */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Predict gene fusion between targeted genes.                        */
/*--------------------------------------------------------------------*/

#include "predict.h"


static kmer_t *kmer_index(fasta_t *, int);
static bag_t  *bag_construct(kmer_t *, fasta_t *, char*, char*, int, int, int);
static char *concat_exons(char* _read, fasta_t *fa_ht, kmer_t *kmer_ht, int _k, char *gname1, char* gname2, char** ename1, char** ename2, int *junction, int min_kmer_match);
static int find_junction_one_edge(bag_t *eg, fasta_t *fasta_u, opt_t *opt, junction_t **ret);
static int align_to_transcript_unit(junction_t *junc, opt_t *opt, solution_pair_t **sol_pair);
static int gene_order(char* gname1, char* gname2, char* read1, char* read2, kmer_t *kmer_ht, int k, int min_kmer_match);
static junction_t *transcript_construct_no_junc(char* gname1, char *gname2, fasta_t *fasta_ht);
static junction_t *transcript_construct_junc(junction_t *junc_ht, fasta_t *exon_ht);
static inline int find_all_genes(str_ctr **hash, kmer_t *KMER_HT, char* _read, int _k);

/*
 * Description:
 *------------
 * index input sequences by kmer hash table

 * Input: 
 *-------
 * fa        - fasta_t hash table contains sequences to be indexed
 * k         - length of kmer

 * Output: 
 *-------
 * kmer_t hash table that contains kmer and its occurnace positions on input seq.
 */
static kmer_t 
*kmer_index(fasta_t *tb, int k){
	if(tb == NULL || k <= 0 || k > MAX_KMER_LEN) return NULL;
	register char *kmer;
	char *name = NULL;	
	char *seq = NULL;	
	register int i, j;
	kmer_t  *s_kmer, *tmp_kmer, *ret = NULL;
	fasta_t *fa_cur;
	for(fa_cur=tb; fa_cur!=NULL; fa_cur=fa_cur->hh.next){
		seq = strToUpper(fa_cur->seq);
		name = strdup(fa_cur->name);
		if(seq == NULL || name == NULL || strlen(seq) <= k) continue;
		for(i=0; i < strlen(seq)-k+1; i++){
			kmer = mycalloc(k+1, char);
			memset(kmer, '\0', k+1);
			strncpy(kmer, seq+i, k);
			kmer_add(&ret, kmer, name); 
		}		
	}
	if(kmer) free(kmer);
	if(seq)  free(seq);
	if(name)  free(name);
	kmer_uniq(&ret);
	return ret;
}

/* 
 * Find all genes uniquely matched with kmers on _read.          
 * hash     - a hash table count number of matches between _read and every gene
 * _read    - inqury read
 * _k       - kmer length
 */
static inline int
find_all_genes(str_ctr **hash, kmer_t *kmer_ht, char* _read, int _k){
	/* check parameters */
	if(_read == NULL || kmer_ht == NULL || _k < 0) die("[%s]: parameter error\n", __func__);
	/* declare vaiables */
	str_ctr *s;
	int _read_pos = 0;
	int num;
	char* gene = NULL;
	register kmer_t *s_kmer = NULL; 
	char buff[_k];
/*--------------------------------------------------------------------*/
	while(_read_pos<(strlen(_read)-_k+1)){
		/* copy a kmer of string */
		strncpy(buff, _read + _read_pos, _k); buff[_k] = '\0';	
		if(strlen(buff) != _k) continue;

		if((s_kmer=find_kmer(kmer_ht, buff)) == NULL){_read_pos++; continue;} // kmer not in table but not an error
		if(s_kmer->count == 1){ // only count the uniq match 
			gene = strsplit(s_kmer->seq_names[0], '.', &num)[0];
			if(gene != NULL){ str_ctr_add(hash, gene);}
		}
		_read_pos++;
	}
	return 0;
}

/*
 * Description:
 *------------
 * construct breakend associated graph BAG.

 * Input: 
 *-------
 * kmer_ht            - kmer hash table returned by kmer_index
 * fq1                - 5' to 3' end of read
 * fq2                - the other end of read
 * min_kmer_matches   - min number kmer matches between a gene and read needed 
 * min_edge_weight    - edges in the graph with weight smaller than min_edge_weight will be deleted
 * k                  - length of kmer
 * Output: 
 *-------
 * BAG_uthash object that contains the graph.
 */
static bag_t
*bag_construct(kmer_t *kmer_ht, fasta_t *fasta_ht, char* fq1, char* fq2, int min_kmer_matches, int min_edge_weight, int _k){
	if(kmer_ht==NULL || fq1==NULL || fq2==NULL || fasta_ht==NULL) return NULL;
	bag_t *bag = NULL;
	gzFile fp1, fp2;
	register int l1, l2;
	register kseq_t *seq1, *seq2;
	register char *_read1, *_read2, *edge_name;
	_read1 = _read2 = edge_name = NULL;
	int i, num;

	if((fp1 = gzopen(fq1, "r"))==NULL) die("[%s] fail to read fastq files", __func__);
	if((fp2 = gzopen(fq2, "r"))==NULL) die("[%s] fail to read fastq files", __func__);	
	if((seq1 = kseq_init(fp1))==NULL)  die("[%s] fail to read fastq files", __func__);
	if((seq2 = kseq_init(fp2))==NULL)  die("[%s] fail to read fastq files", __func__);
		
	/* iterate read pair in both fastq files */
	while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >= 0 ) {
		_read1 = rev_com(seq1->seq.s); // reverse complement of read1
		_read2 = seq2->seq.s;		
		if(_read1 == NULL || _read2 == NULL) continue;
		//if(strcmp(seq1->name.s, seq2->name.s) != 0) die("[%s] read pair not matched", __func__);		
		if(strlen(_read1) < _k || strlen(_read2) < _k){continue;}
		/* gene_counter is a string counter which will count the number of times a gene occurs in read */
		str_ctr *gene_counter, *s;
		gene_counter = NULL;
		find_all_genes(&gene_counter, kmer_ht, _read1, _k);
		find_all_genes(&gene_counter, kmer_ht, _read2, _k);
		
		if((num = HASH_COUNT(gene_counter))<2) continue; // if less than two genes identified, pass the rest
		char** hits = mycalloc(num, char*);
		/* filter genes that have matches with kmer less than min_kmer_matches */
		i=0; for(s=gene_counter; s!=NULL; s=s->hh.next){if(s->SIZE >= min_kmer_matches){hits[i++] = strdup(s->KEY);}}
			
		int m, n; for(m=0; m < i; m++){for(n=m+1; n < i; n++){
				int rc = strcmp(hits[m], hits[n]);
				if(rc<0)  edge_name = concat(concat(hits[m], "_"), hits[n]);
				if(rc>0)  edge_name = concat(concat(hits[n], "_"), hits[m]);
				if(rc==0) edge_name = NULL;
				if(edge_name!=NULL) if(bag_add(&bag, edge_name, seq1->name.s, concat(concat(_read1, "_"), _read2)) != 0) die("BAG_uthash_add fails\n");							
		}}
		// clean the mess up
		if(hits)		 for(i--;i>0; i--){if(hits[i]) free(hits[i]);};
		if(gene_counter) str_ctr_destory(&gene_counter);
	}
	// remove the edges with weight < min_edge_weight	
	register bag_t *cur, *tmp;
	HASH_ITER(hh, bag, cur, tmp) {
		if(cur->weight < min_edge_weight) {
			HASH_DEL(bag, cur); 
			free(cur);
		}
	}
	// determine gene order by kmer matches
	int order;
	for(cur=bag; cur!=NULL; cur=cur->hh.next){
		order = 0;
		for(i=0; i<cur->weight; i++){
			_read1 = strsplit(cur->evidence[i], '_', &num)[0]; _read2 = strsplit(cur->evidence[i], '_', &num)[1];	
			order += gene_order(strsplit(cur->edge, '_', &num)[0], strsplit(cur->edge, '_', &num)[1], _read1, _read2, kmer_ht, _k, min_kmer_matches);
		}
		if(order > 0){
			cur->gname1 = strdup(strsplit(cur->edge, '_', &num)[1]); 
			cur->gname2 = strdup(strsplit(cur->edge, '_', &num)[0]);
		}
		if(order < 0){
			cur->gname1 = strdup(strsplit(cur->edge, '_', &num)[0]); 
			cur->gname2 = strdup(strsplit(cur->edge, '_', &num)[1]);
		}
		if(order == 0){ HASH_DEL(bag, cur); free(cur);}		
	}	

	// clean the mess up
	if(_read1)      free(_read1);
	if(_read2)      free(_read2);
	if(edge_name)   free(edge_name);
	kseq_destroy(seq1);
	kseq_destroy(seq2);	
	gzclose(fp1);
	gzclose(fp2);
	return bag;
}


/*
 * determine the order of fusion genes.
 * negative means gene1 in front of gene1 from 5'-3'
 */
static int 
gene_order(char* gname1, char* gname2, char* read1, char* read2, kmer_t *kmer_ht, int k, int min_kmer_match){
	if(gname1==NULL || gname2==NULL || read1==NULL || read2==NULL || kmer_ht==NULL) return 0;
	register int i;
	int *gene1 = mycalloc(strlen(read1)+strlen(read2), int);
	int *gene2 = mycalloc(strlen(read1)+strlen(read2), int);
	int gene1_pos = 0;	
	int gene2_pos = 0;	
	int num;
	char* gname_tmp;
	register kmer_t *kmer_cur;
	char *buff = mycalloc(k+1, char);
	for(i=0; i<strlen(read1)-k+1; i++){
		memset(buff, '\0', k+1);
		strncpy(buff, read1 + i, k); buff[k] = '\0';			
		if((kmer_cur=find_kmer(kmer_ht, buff)) == NULL) continue;
		if(kmer_cur->count == 1){  // uniq match
			gname_tmp = strsplit(kmer_cur->seq_names[0], '.', &num)[0];
			if(strcmp(gname_tmp, gname1)==0) gene1[gene1_pos++] = i;
			if(strcmp(gname_tmp, gname2)==0) gene2[gene2_pos++] = i;
		}
	}
	for(i=0; i<strlen(read2)-k+1; i++){
		memset(buff, '\0', k+1);
		strncpy(buff, read2 + i, k); buff[k] = '\0';			
		if((kmer_cur=find_kmer(kmer_ht, buff)) == NULL) continue;
		if(kmer_cur->count == 1){  // uniq match
			gname_tmp = strsplit(kmer_cur->seq_names[0], '.', &num)[0];
			if(strcmp(gname_tmp, gname1)==0) gene1[gene1_pos++] = i+strlen(read1);
			if(strcmp(gname_tmp, gname2)==0) gene2[gene2_pos++] = i+strlen(read1);
		}
	}
	int t = 0;
	if(gene1_pos >= min_kmer_match && gene2_pos >= min_kmer_match){
		for(i=0; i<min(gene1_pos, gene2_pos); i++){
			if(gene1[i] < gene2[i]) t++;
			if(gene1[i] > gene2[i]) t--;
		}
		return t;
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
find_all_exons(str_ctr **hash, kmer_t *KMER_HT, char* _read, int _k){
/*--------------------------------------------------------------------*/
	/* check parameters */
	if(_read == NULL || _k < 0) die("find_all_MEKMs: parameter error\n");
/*--------------------------------------------------------------------*/
	/* declare vaiables */
	str_ctr *s;
	int _read_pos = 0;
	char* exon = NULL;
	register kmer_t *s_kmer = NULL; 
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

static junction_t
*edge_junction_gen(bag_t *eg, fasta_t *fasta_u, opt_t *opt){
	if(eg==NULL || fasta_u==NULL || opt==NULL) return NULL;
	/* variables */
	int _k = opt->k;
	int num;
	register char* _read1, *_read2, *str2;
	char *ename1, *ename2;
	char *chrom1, *chrom2;
	char* idx = NULL;
	chrom2 = chrom1 = NULL;
	ename1 = ename2 = NULL;
	_read1 = _read2 = str2 = NULL;
	char *gname1, *gname2;
	gname1 = eg->gname1;
	gname2 = eg->gname2;
	if(gname1==NULL || gname2==NULL) return NULL;
	register int i, j;
	register kmer_t *s_kmer; 
	register solution_t *sol1, *sol2;           /* alignment solution for read1 and read2 */
	int start1, start2;
	int junc_pos;                               /* position of junction */
	int strlen2;
	junction_t *m, *n, *ret = NULL;

	for(i=0; i<eg->weight; i++){
		_read1 = strsplit(eg->evidence[i], '_', &num)[0];
		_read2 = strsplit(eg->evidence[i], '_', &num)[1];	
		/* string concatnated by exon sequences of two genes */
		if((str2 =  concat_exons(_read1, fasta_u, KMER_HT, _k, gname1, gname2, &ename1, &ename2, &junc_pos, opt->min_kmer_match))==NULL) continue;
		/* fiting alignment that allows jump state */
		sol1 = align(_read1, str2, junc_pos, opt->match, opt->mismatch, opt->gap, opt->extension, opt->jump_gene);
		sol2 = align(_read2, str2, junc_pos, opt->match, opt->mismatch, opt->gap, opt->extension, opt->jump_gene);

		if(sol1->jump == true && sol1->prob >= opt->min_align_score){
			/* idx = exon1.start.exon2.end (uniq id)*/
			idx = idx2str(concat(concat(ename1, "."), ename2), sol1->jump_start, sol1->jump_end); // idx for junction
			HASH_FIND_STR(ret, idx, m);  
			if(m==NULL){ // add this junction to ret
				m = junction_init(opt->seed_len);				
				m->idx    = idx;
				m->exon1  = ename1;
				m->exon2  = ename2;				
				m->hits  = 1;
				m->likehood = 10*log(sol1->prob); 				
				// junction flanking sequence 
				memcpy( m->s, &str2[sol1->jump_start-opt->seed_len/2-1], opt->seed_len/2);
				memcpy( &m->s[opt->seed_len/2], &str2[sol1->jump_end], opt->seed_len/2);
				// junction flanking sequence 
				strlen2 = sol1->jump_start + strlen(str2) - sol1->jump_end+1;				
				m->transcript = mycalloc(strlen2, char);
				m->junc_pos = sol1->jump_start;				
				memset(m->transcript, '\0', strlen2);
				memcpy( m->transcript, str2, sol1->jump_start);
				memcpy( &m->transcript[sol1->jump_start], &str2[sol1->jump_end+1], strlen(str2) - sol1->jump_end);				
				HASH_ADD_STR(ret, idx, m);
			}else{
				m->hits ++;
				m->likehood += 10*log(sol1->prob); 
			}
		}
		if(sol2->jump == true && sol2->prob >= opt->min_align_score){
			idx = idx2str(concat(concat(ename1, "."), ename2), sol2->jump_start, sol2->jump_end);
			HASH_FIND_STR(ret, idx, m);
			if(m==NULL){ // this junction not in ret
				m = junction_init(opt->seed_len);				
				m->idx   = idx;
				m->exon1  = ename1;
				m->exon2  = ename2;				
				m->hits  = 1;
				m->likehood = 10*log(sol2->prob); 				
				memcpy( m->s, &str2[sol2->jump_start-opt->seed_len/2-1], opt->seed_len/2);
				memcpy( &m->s[opt->seed_len/2], &str2[sol2->jump_end], opt->seed_len/2);
				strlen2 = sol2->jump_start + strlen(str2) - sol2->jump_end+1;
				m->transcript = mycalloc(strlen2, char);
				m->junc_pos = sol2->jump_start;
				memset(m->transcript, '\0', strlen2);
				memcpy( m->transcript, str2, sol2->jump_start);
				memcpy( &m->transcript[sol2->jump_start], &str2[sol2->jump_end], strlen(str2)-sol2->jump_end);				
				HASH_ADD_STR(ret, idx, m);
			}else{
				m->hits ++;
				m->likehood += 10*log(sol2->prob); 
			}
		}	
	}
	// delete those junctions with hits < MIN_HITS
	HASH_ITER(hh, ret, m, n){
		if(m != NULL){
			m->likehood = m->likehood/m->hits;
			if(m->hits < opt->min_hits){
				HASH_DEL(ret,m);
				free(m);
			}
		}
	}
	if(ename1)         free(ename1);
	if(ename2)         free(ename2);
	if(str2)           free(str2);
	if(_read1)         free(_read1);
	if(_read2)         free(_read2);
	if(sol1)           solution_destory(sol1);
	if(sol2)           solution_destory(sol2);
	if(idx)            free(idx);
	return ret;
}

static junction_t 
*transcript_construct_junc(junction_t *junc_ht, fasta_t *exon_ht){
	if(junc_ht==NULL || exon_ht==NULL) return NULL;
	junction_t *cur_junction, *tmp_junction;
	char* gname1, *gname2, *ename1, *ename2;
	gname1 = gname2 = ename1 = ename2 = NULL;
	int enum1, enum2;
	char enum1_str[10], enum2_str[10];
	int num1, num2;
	char** fields1, **fields2;
	fasta_t *exon1_fa,  *exon2_fa;
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

static junction_t
*transcript_construct_no_junc(char* gname1, char *gname2, fasta_t *fasta_ht){
	if(gname1==NULL || gname2==NULL || fasta_ht==NULL) return NULL;	
	junction_t *junc_res = junction_init(0);
	if(gname1==NULL || gname2==NULL || fasta_ht==NULL) return NULL;
	fasta_t *fasta_cur;
	char *res1 = NULL;
	char *res2 = NULL;
	int num;
	int idx = 0;
	char** strs;
	int *S1, *S2;
	S1 = S2 = NULL;
	for(fasta_cur=fasta_ht; fasta_cur!=NULL; fasta_cur=fasta_cur->hh.next){
		S1 = mycalloc(100, int); memset(S1, 0, 100);		
		strs = strsplit(fasta_cur->name, '.', &num);
		if(num!=2) continue;
		if(strcmp(strs[0], gname1) != 0) continue;
		if(atoi(strs[1]) < idx) continue; // make sure the order
		res1 = concat(res1, fasta_cur->seq);
		S1[idx] = strlen(res1);
		idx = atoi(strs[1]);		
	}
	junc_res->S1_num = idx;
	
	idx = 0;
	for(fasta_cur=fasta_ht; fasta_cur!=NULL; fasta_cur=fasta_cur->hh.next){
		S2 = mycalloc(100, int); memset(S2, UINT_MAX, 100);
		strs = strsplit(fasta_cur->name, '.', &num);
		if(num!=2) continue;
		if(strcmp(strs[0], gname2) != 0) continue;
		if(atoi(strs[1]) < idx) continue; // make sure the order
		res2 = concat(res2, fasta_cur->seq);
		S2[idx] = strlen(res2);
		idx = atoi(strs[1]);
	}
	junc_res->S2_num = idx;
	junc_res->idx = concat(concat(gname1, "_"), gname2);
	junc_res->exon1 = NULL;
	junc_res->exon2 = NULL;
	junc_res->s = NULL;	
	junc_res->transcript = concat(res1, res2);	
	junc_res->junc_pos = strlen(res1);
	junc_res->S1 = S1;
	junc_res->S2 = S2;	
	junc_res->hits = 0;
	junc_res->likehood = -INFINITY;
	return junc_res;
}

/*
 * generate junction string of every edge based on supportive reads.
 */
static int
bag_junction_gen(bag_t **bag, fasta_t *fa, opt_t *opt){
	if(*bag==NULL || fa==NULL || opt==NULL) return -1;	
	bag_t *edge, *bag_cur;
	register int i;
	junction_t *junc_cur;
	for(edge=*bag; edge!=NULL; edge=edge->hh.next) {
		if((junc_cur = edge_junction_gen(edge, fa, opt))==NULL){ // no junction detected
			edge->junc_flag = false;
			edge->junc = NULL;
		}else{
			edge->junc_flag = true;
			edge->junc = junc_cur;			
		}		
	}   
	return 0;
}


/*
 * generate junction string of every edge based on supportive reads.
 */
static int
bag_transcript_gen(bag_t **bag, fasta_t *fa, opt_t *opt){
	if(*bag==NULL || fa==NULL || opt==NULL) return -1;	
	bag_t *edge, *bag_cur;
	register int i;
	junction_t *junc_cur;
	for(edge=*bag; edge!=NULL; edge=edge->hh.next) {
		if(edge->junc_flag==true){
			printf("sfsdfsadf\n");
			edge->junc = transcript_construct_junc(edge->junc, fa);
		}else{
			edge->junc = transcript_construct_no_junc(edge->gname1, edge->gname2, fa);
		}
	}   
	return 0;
}
/*
 * construct concatnated exon string based on kmer matches
 */
static char 
*concat_exons(char* _read, fasta_t *fa_ht, kmer_t *kmer_ht, int _k, char *gname1, char* gname2, char** ename1, char** ename2, int *junc_pos, int min_kmer_match){
	if(_read == NULL || fa_ht == NULL || kmer_ht==NULL || gname1==NULL || gname2==NULL) return NULL;
	/* variables */
	char *str1, *str2;
	*ename1 = *ename2 = str1 = str2 = NULL;
	int num_tmp;
	char buff[_k];
	char* gname_cur;
	str_ctr *s_ctr, *tmp_ctr, *exons=NULL;
	fasta_t *fa_tmp;
	/* find all exons that uniquely match with gene by kmer */
	find_all_exons(&exons, kmer_ht, _read, _k);
	if(exons==NULL) return NULL; // no exon found
	HASH_ITER(hh, exons, s_ctr, tmp_ctr) { 
		if(s_ctr->SIZE >= min_kmer_match){ //denoise
			gname_cur = strsplit(s_ctr->KEY, '.', &num_tmp)[0]; // name of gene
			// extract string of gene1
			if(strcmp(gname_cur, gname1)==0){
				fa_tmp = find_fasta(fa_ht, s_ctr->KEY);
				if(str1 == NULL){
					str1 = strdup(fa_tmp->seq);
					*ename1 = s_ctr->KEY;
				}else{
					str1 = concat(str1, fa_tmp->seq);
				}				
			}
			// extract string of gene2
			if(strcmp(gname_cur, gname2)==0){
				fa_tmp = find_fasta(fa_ht, s_ctr->KEY);
				if(str2 == NULL){
					str2 = strdup(fa_tmp->seq);
					*ename2 = s_ctr->KEY;
				}else{
					*ename2 = s_ctr->KEY;
					str2 = concat(str2, fa_tmp->seq);
				}				
			}
		}
	}
	if(gname_cur) free(gname_cur);
	if(exons)         str_ctr_destory(&exons);
	if(str1 == NULL || str2 == NULL) return NULL; // only if two genes identified
	char *ret = concat(str1, str2);
	*junc_pos = strlen(str1);	
	return ret;
}

/*
 * align supportive reads of every edge to edge's constructed transcript
 */
static solution_pair_t 
*align_edge_to_transcript(bag_t *bag_ht, opt_t *opt){
	if(bag_ht==NULL || opt==NULL) return NULL;
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
	// iterate every junction in edge
	for(edge=bag_ht; edge!=NULL; edge=edge->hh.next){
		for(junc_cur=edge->junc; junc_cur != NULL; junc_cur=junc_cur->hh.next){
			gene1 = strsplit(junc_cur->exon1, '.', &num)[0];
			gene2 = strsplit(junc_cur->exon2, '.', &num)[0];
			rc = strcmp(gene1, gene2);
			if(rc<0)  edge_name = concat(concat(gene1, "_"), gene2);
			if(rc>0)  edge_name = concat(concat(gene1, "_"), gene2);
			if(rc==0) continue;
			if((edge = find_edge(bag_ht, edge_name))==NULL) continue;
			for(i=0; i<edge->weight; i++){
				read1 = strsplit(edge->evidence[i], '_', &num)[0];
				if(num != 2) continue;
				read2 = strsplit(edge->evidence[i], '_', &num)[1];
				if(junc_cur->transcript==NULL) break;
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
static int align_reads_to_transcript(solution_pair_t **res, bag_t *bag, opt_t *opt){
	if(bag==NULL || opt==NULL) return -1;
	bag_t *bag_cur;
	junction_t *junc_cur;
	for(bag_cur=bag; bag_cur!=NULL; bag_cur=bag_cur->hh.next){		
		if(bag_cur->junc_flag==false) continue;
		fprintf(stderr, "[predict] junction between %s and %s is being tested ... \n", bag_cur->gname1, bag_cur->gname2);		
		for(junc_cur=bag_cur->junc; junc_cur!=NULL; junc_cur=junc_cur->hh.next){
			if((align_to_transcript_unit(junc_cur, opt, res))!=0) return -1;				
		}
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
		//if(strcmp(seq1->name.s, seq2->name.s) != 0) die("[%s] read pair not matched\n", __func__);
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
			fprintf(stderr, "Options: -k INT    kmer length for indexing [%d]\n", opt->k);
			fprintf(stderr, "         -n INT    min uniq kmer matches for a gene and read match [%d]\n", opt->min_kmer_match);
			fprintf(stderr, "         -w INT    edges in graph of weight smaller than -w will be removed [%d]\n", opt->min_edge_weight);
			fprintf(stderr, "         -m INT    score for match in alignment [%d]\n", opt->match);
			fprintf(stderr, "         -u INT    score for mismatch in alignment [%d]\n", opt->mismatch);
			fprintf(stderr, "         -o INT    penality for gap open [%d]\n", opt->gap);
			fprintf(stderr, "         -e INT    penality for gap extension [%d]\n", opt->extension);
			fprintf(stderr, "         -j INT    penality for jump between genes [%d]\n", opt->jump_gene);
			fprintf(stderr, "         -s INT    penality for jump between exons [%d]\n", opt->jump_exon);
			fprintf(stderr, "         -h INT    min hits for a junction [%d]\n", opt->min_hits);					
			fprintf(stderr, "         -l INT    length for junction string [%d]\n", opt->seed_len);					
			fprintf(stderr, "         -x INT    max mismatches of junction string match [%d]\n", opt->max_mismatch);					
			fprintf(stderr, "         -a FLOAT  min identity score for alignment [%.2f]\n", opt->min_align_score);
			fprintf(stderr, "\n");
			fprintf(stderr, "Inputs:  exon.fa   fasta file that contains exon sequences of \n");
			fprintf(stderr, "                   targeted genes, which can be generated by: \n");
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
				case 'k': opt->k = atoi(optarg); break;	
				case 'n': opt->min_kmer_match = atoi(optarg); break;
				case 'w': opt->min_edge_weight = atoi(optarg); break;
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
	
	if(opt->k < MIN_KMER_LEN || opt->k > MAX_KMER_LEN) die("[%s] -k must be within [%d, %d]", __func__, MIN_KMER_LEN, MAX_KMER_LEN); 	
	if(opt->min_kmer_match < MIN_KMER_MATCH) die("[%s] -n must be within [%d, +INF)", __func__, MIN_KMER_MATCH); 	
	if(opt->min_edge_weight < MIN_EDGE_WEIGHT) die("[%s] -w must be within [%d, +INF)", __func__, MIN_EDGE_WEIGHT); 	
	
	fprintf(stderr, "[%s] loading sequences of targeted genes ... \n",__func__);
	if((EXON_HT = fasta_read(opt->fa)) == NULL) die("[%s] fail to read %s", __func__, opt->fa);	
	
	fprintf(stderr, "[%s] indexing sequneces by kmer hash table ... \n",__func__);
	if((KMER_HT = kmer_index(EXON_HT, opt->k))==NULL) die("[%s] can't index exon sequences", __func__); 	

	fprintf(stderr, "[%s] constructing breakend associated graph ... \n", __func__);
	if((BAGR_HT = bag_construct(KMER_HT, EXON_HT, opt->fq1, opt->fq2, opt->min_kmer_match, opt->min_edge_weight, opt->k)) == NULL) return 0;
	
	
	fprintf(stderr, "[%s] triming graph by removing duplicate supportive read pairs of each edge ... \n", __func__);
	if(bag_uniq(&BAGR_HT)!=0){
		fprintf(stderr, "[%s] fail to remove duplicate supportive reads \n", __func__);
		return -1;		
	}
	bag_display(BAGR_HT);
	//if(BAGR_HT == NULL) return 0;
    //
	//fprintf(stderr, "[%s] triming graph by removing edges of weight smaller than %d... \n", __func__, opt->min_edge_weight);
	//if(bag_trim(&BAGR_HT, opt->min_edge_weight)!=0){
	//	fprintf(stderr, "[%s] fail to trim graph \n", __func__);
	//	return -1;
	//}
	//if(BAGR_HT == NULL) return 0;	
    //
	//fprintf(stderr, "[%s] identifying junctions for every fusion candiates... \n", __func__);
	//if(bag_junction_gen(&BAGR_HT, EXON_HT, opt)!=0){
	//	fprintf(stderr, "[%s] fail to identify junctions\n", __func__);
	//	return -1;	
	//}
	//if(BAGR_HT == NULL) return 0;
    //
	//fprintf(stderr, "[%s] constructing fused transcript for every fusion candiates... \n", __func__);
	//if(bag_transcript_gen(&BAGR_HT, EXON_HT, opt)!=0){
	//	fprintf(stderr, "[%s] fail to construct transcript  \n", __func__);
	//	return -1;	
	//}
	//if(BAGR_HT == NULL) return 0;
	//
	//bag_display(BAGR_HT);
	//fprintf(stderr, "[%s] aligning supportive read pairs to fused transcript ... \n", __func__);    
	//if((SOLU_HT = align_edge_to_transcript(BAGR_HT, opt))==NULL) die("[%s] can't rediscover any junction", __func__);
	//
	//fprintf(stderr, "[%s] testing junctions ... \n", __func__);    	
	//if((align_reads_to_transcript(&SOLU_HT, BAGR_HT, opt)) != 0) die("[%s] can't rediscover any junction", __func__);;
	//
	//solution_pair_t *s;
	//for(s=SOLU_HT; s!=NULL; s=s->hh.next){
	//	printf("%s\n", s->junc_name);
	//}
	//
	fprintf(stderr, "[%s] cleaning up ... \n", __func__);	
	if(EXON_HT)          fasta_destroy(&EXON_HT);
	if(KMER_HT)           kmer_destroy(&KMER_HT);
	if(BAGR_HT)            bag_distory(&BAGR_HT);
	if(SOLU_HT)  solution_pair_destory(&SOLU_HT);
	return 0;
}

