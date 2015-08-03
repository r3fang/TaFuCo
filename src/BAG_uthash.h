/*--------------------------------------------------------------------*/
/* Created Date: 15JULY2015                                           */
/* Author: Rongxin Fang                                               */
/* Contact: r3fang@ucsd.edu                                           */
/* Library for Breakend Associated Graph (BAG).                       */
/* Functions it contain:                                              */
/*--------------------------------------------------------------------*/

#ifndef _BAG_UTHASH_H
#define _BAG_UTHASH_H

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include "utils.h"
#include "kseq.h"
#include "kmer_uthash.h"
/* error code */
#define BA_ERR_NONE		     0 // no error
/*
 * the BAG_uthash structure
 */
struct BAG_uthash {
	char *edge;
	size_t weight;
	char **evidence;    
    UT_hash_handle hh;         /* makes this structure hashable */
};

/*
 * intiate BAG_uthash
 *
 */
static inline struct BAG_uthash 
*BAG_uthash_init() {
	struct BAG_uthash *t = mycalloc(1, struct BAG_uthash);
	t->edge = NULL;
	t->weight = 0;
	t->evidence = mycalloc(1, char*);
	return t;
}

/*
 * destory 
 */
static inline int 
BAG_uthash_destroy(struct BAG_uthash **table) {
	if(*table == NULL) die("BAG_uthash_destroy: parameter error\n");
	/*free the kmer_hash table*/
	register struct BAG_uthash *cur, *tmp;
	HASH_ITER(hh, *table, cur, tmp) {
		if(cur==NULL) die("BAG_uthash_destroy: HASH_ITER fails\n");
		HASH_DEL(*table, cur); 
		free(cur);   
    }
	return BA_ERR_NONE;
}
/*
 * display BAG_uthash table on the screen
 *
 * PARAMETERS:	struct BAG_uthash *
 * RETURN:	error code
 */
static inline int BAG_uthash_display(struct BAG_uthash *graph_ht) {
	if(graph_ht == NULL) die("BAG_uthash_display: parameter error\n");
	/*free the kmer_hash table*/
	register struct BAG_uthash *cur, *tmp;
	HASH_ITER(hh, graph_ht, cur, tmp) {
		if(cur == NULL) die("BAG_uthash_display: HASH_ITER fails\n");
		int i; for(i=0; i < cur->weight; i++){
			printf(">%s\n%s\n", cur->edge, cur->evidence[i]);
		}
	}
	return BA_ERR_NONE;
}

static inline struct BAG_uthash 
*find_edge(struct BAG_uthash *tb, char* quary_name) {
	if(quary_name == NULL) die("[%s] input error", __func__);
	struct BAG_uthash* s = NULL;	
    HASH_FIND_STR(tb, quary_name, s);  /* s: output pointer */
	return s;
}

/*
 * add one edge to graph
 */
static inline int 
BAG_uthash_add(struct BAG_uthash** graph_ht, char* edge_name, char* evidence){
	/* check parameters */
	if(edge_name == NULL || evidence == NULL) die("[%s]: parameter error\n", __func__);
	struct BAG_uthash *s;
	register int n;
	if((s = find_edge(*graph_ht, edge_name)) == NULL){
		s = BAG_uthash_init();
		s->edge = strdup(edge_name);
		s->weight = 1;
		s->evidence[0] = strdup(evidence); /* first and only 1 element*/
		HASH_ADD_STR(*graph_ht, edge, s);								
	}else{
		s->weight++;
		char **tmp = mycalloc(s->weight, char*);
	 	for (n = 0; n < s->weight-1; n++){
	 		tmp[n] = strdup(s->evidence[n]);
	 	}
	 	tmp[n] = strdup(evidence);
	 	free(s->evidence);
	 	s->evidence = tmp;
	}
	return BA_ERR_NONE;
}

static inline int 
cmpstr(void const *a, void const *b) { 
    char const *aa = (char const *)a;
    char const *bb = (char const *)b;
    return strcmp(aa, bb);
}

/*
 * trim edges with evidence less than min_weight
 */
static inline int 
BAG_uthash_trim(struct BAG_uthash** tb, int min_weight){
	if(*tb == NULL || min_weight < 0) die("BAG_uthash_tim: wrong parameters\n");
	register struct BAG_uthash *cur, *tmp;
	HASH_ITER(hh, *tb, cur, tmp) {
		if(cur==NULL) die("BAG_uthash_tim: HASH_ITER fails\n");
		if(cur->weight < min_weight) {HASH_DEL(*tb, cur); free(cur);}
    }
	return BA_ERR_NONE;
}

/*
 * rm duplicate evidence for graph edge
 */
static inline int 
BAG_uthash_uniq(struct BAG_uthash **tb){
	if(*tb == NULL) die("BAG_uthash_uniq: input error");
	struct BAG_uthash *cur, *tmp;
	int i, j;
	int before, del;
	HASH_ITER(hh, *tb, cur, tmp) {
		if(cur == NULL) die("kmer_uthash_uniq: HASH_ITER fails\n");
		qsort(cur->evidence, cur->weight, sizeof(char*), mystrcmp);
		before = cur->weight;
		del = 0;
		if(cur->weight > 1){
			for(i=1; i<cur->weight; i++){
				if(mystrcmp(cur->evidence[i], cur->evidence[i-1]) == 0){
					cur->evidence[i-1] = NULL;
					del++;
				}
			}
		}		
		if(del > 0){ // if any element has been deleted
			char** tmp = mycalloc(before - del , char*);
			j=0; for(i=0; i<cur->weight; i++){
				if(cur->evidence[i] != NULL){
					tmp[j++] = strdup(cur->evidence[i]);
				}
			}
			free(cur->evidence);
			cur->evidence = tmp;		
			cur->weight = before - del;	
		}
    }	
	return BA_ERR_NONE;
}

///*
// * trim edges with evidence less than min_weight
// */
static inline struct BAG_uthash 
*BAG_uthash_load(char* fname){
	if(fname == NULL) die("[%s] input error", __func__);
	struct BAG_uthash *tb = NULL;
	struct BAG_uthash *cur = NULL;
	
	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen(fname, "r");
	if(fp == NULL) die("[%s] fail to open %s", __func__, fname);		

	struct fasta_uthash *s;	
	if((seq = kseq_init(fp))==NULL) die("[%s] kseq_init fails", __func__);

	while ((l = kseq_read(seq)) >= 0){
		if(seq->name.s == NULL || seq->seq.s==NULL)
			continue;
		BAG_uthash_add(&tb, seq->name.s, seq->seq.s);		
	}	
	//if(cur) BAG_uthash_destroy(&cur);
	if(seq) kseq_destroy(seq);
	gzclose(fp);
	return tb;
}

/* 
 * Find all genes uniquely matched with kmers on _read.          
 * hash     - a hash table count number of matches between _read and every gene
 * _read    - inqury read
 * _k       - kmer length
 */
static inline int
find_all_matches(str_ctr **hash, struct kmer_uthash *KMER_HT, char* _read, int _k){
/*--------------------------------------------------------------------*/
	/* check parameters */
	if(_read == NULL || _k < 0) die("find_all_MEKMs: parameter error\n");
/*--------------------------------------------------------------------*/
	/* declare vaiables */
	str_ctr *s;
	int _read_pos = 0;
	char* gene = NULL;
	struct kmer_uthash *s_kmer = NULL; 
	char buff[_k];
/*--------------------------------------------------------------------*/
	while(_read_pos<(strlen(_read)-_k-1)){
		/* copy a kmer of string */
		strncpy(buff, _read + _read_pos, _k); buff[_k] = '\0';	
		if(strlen(buff) != _k) die("find_next_match: buff strncpy fails\n");
		/*------------------------------------------------------------*/
		if((s_kmer=find_kmer(KMER_HT, buff)) == NULL){_read_pos++; continue;} // kmer not in table but not an error
		if(s_kmer->count == 1){ // only count the uniq match 
			gene = strdup(s_kmer->seq_names[0]);
			if(gene == NULL) die("find_next_match: get_exon_name fails\n");
			if(str_ctr_add(hash, gene) != 0) die("find_all_MEKMs: str_ctr_add fails\n");
		}
		_read_pos++;
	}
	return 0;
}

/* 
 * Construct breakend associated graph.             
 * fq_file1 - fastq file name for R1
 * fq_file2 - fastq file name for R2
 * _k       - kmer length for hash table
 * cutoff   - min number of kmer matches between a read again a gene.
 * bag      - BAG_uthash object (breakend associated graph)
 */

static inline int
construct_BAG(struct BAG_uthash **bag, struct kmer_uthash *kmer_uthash, char* fq1, char* fq2, int _min_match, int _k){
	if(bag == NULL || kmer_uthash == NULL) die("[%s] input error\n", __func__);
	
	int error;
	gzFile fp1, fp2;
	register int l1, l2;
	register kseq_t *seq1, *seq2;
	register char *_read1, *_read2;
	register char *edge_name;
	_read1 = _read2 = edge_name = NULL;
	
	if((fp1 = gzopen(fq1, "r"))==NULL) die("construct_BAG: fail to read fastq files\n");
	if((fp2 = gzopen(fq2, "r"))==NULL) die("construct_BAG: fail to read fastq files\n");	
	if((seq1 = kseq_init(fp1))==NULL) die("construct_BAG: fail to read fastq files\n");
	if((seq2 = kseq_init(fp2))==NULL) die("construct_BAG: fail to read fastq files\n");	

	while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >= 0 ) {
		_read1 = rev_com(seq1->seq.s); // reverse complement of read1
		_read2 = seq2->seq.s;		
		if(_read1 == NULL || _read2 == NULL) die("construct_BAG: fail to get _read1 and _read2\n");
		if(strcmp(seq1->name.s, seq2->name.s) != 0) die("construct_BAG: read pair not matched\n");		
		if(strlen(_read1) < _k || strlen(_read2) < _k){continue;}
		str_ctr* gene_counter = NULL;
		str_ctr *s, *tmp;
		find_all_matches(&gene_counter, kmer_uthash, _read1, _k);
		find_all_matches(&gene_counter, kmer_uthash, _read2, _k);
		
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
					if(BAG_uthash_add(bag, edge_name, concat(concat(_read1, "_"), _read2)) != 0) die("BAG_uthash_add fails\n");							
				}
		}}
		if(hits)		 free(hits);
		if(gene_counter) free(gene_counter);
	}
	if(edge_name)   free(edge_name);
	kseq_destroy(seq1);
	kseq_destroy(seq2);	
	gzclose(fp1);
	gzclose(fp2);
	return 0; // no error raised
}

#endif
