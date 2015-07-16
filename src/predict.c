/*--------------------------------------------------------------------*/
/* predict.c                                                          */
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
#include "uthash.h"
#include "kseq.h"
#include "common.h"
#include "kmer_uthash.h"
#include "fasta_uthash.h"
#include "BAG_uthash.h"
KSEQ_INIT(gzFile, gzread);


/* error code */
#define PR_ERR_NONE		     	 0 // no error
#define PR_ERR_PARAM			-1 // bad paraemter
#define PR_ERR_MALLOC			-2 // fail to malloc memory uthash
#define PR_ERR_POSPARSE			-3
#define PR_ERR_FINDHASH			-4
#define PR_ERR_UNDEFINED		-100 // fail to malloc memory uthash

/*
 * error handling
 * all errors are returned as an integer code, and a string
 * amplifying the error is saved in here; this can then be
 * printed
 */
char pr_errbuf[256] = "no error\n";	/* the error message buffer */
#define ERRBUF(str)			(void) strncpy(qe_errbuf, str, sizeof(qe_errbuf))
#define ERRBUF2(str,n)		(void) sprintf(qe_errbuf, str, n)
#define ERRBUF3(str,n,m)	(void) sprintf(qe_errbuf, str, n, m)

/*--------------------------------------------------------------------*/
/*Global paramters.*/
static const char *pcPgmName="predict.c";
static struct kmer_uthash *KMER_HT = NULL;
static struct fasta_uthash *FASTA_HT = NULL;
static struct BAG_uthash *BAG_HT = NULL;

/*--------------------------------------------------------------------*/

/* Find next Maximal Extended Kmer Match (MEKM) on _read at pos_read. */
int
find_next_MEKM(char *_read, int pos_read, int k, int min_match){
	/*------------------------------------------------------------*/
	// parameters decleration
	typedef struct freq{
		int LEN;
		size_t SIZE;
	    UT_hash_handle hh;         /* makes this structure hashable */
	} freq; 
	freq* match_lens = NULL;
	char* buff = NULL;	
	char **matches = NULL;
	int i = 0;
	int num = 0;
	char *exon_tmp = NULL;
	char *exon_max = NULL;
	int max_len = -10;
	freq *s_freq = NULL;
	freq *tmp = NULL;
	int error;
	/*------------------------------------------------------------*/
	/* check parameters */
	if(read == NULL) goto FAIL_PARAM;
	
	/*------------------------------------------------------------*/
	/* copy a kmer of string */
	if((buff = malloc((k+1) * sizeof(char)))==NULL) goto FAIL_MALLOC;
	strncpy(buff, _read + pos_read, k);
	buff[k] = '\0';	
	if(buff == NULL || strlen(buff) != k) goto FAIL_OTHER;
	/*------------------------------------------------------------*/
	struct kmer_uthash *s_kmer;
	if((s_kmer = find_kmer(buff, KMER_HT)) == NULL) goto SUCCESS;
	/*------------------------------------------------------------*/
	if((matches = malloc(s_kmer->count * sizeof(char*)))==NULL) goto FAIL_MALLOC;	
	/*------------------------------------------------------------*/
	/*
	 * iterate all kmer matches and extend to max length
	 */
	for(i=0; i<s_kmer->count; i++){		
		int _pos_exon; // position on exon
		exon_tmp = pos_parser(s_kmer->pos[i], &_pos_exon);
		/* error report*/
		if(exon_tmp == NULL || _pos_exon < 0) goto FAIL_POSPARSE;
		
		struct fasta_uthash *s_fasta;
		if((s_fasta = find_fasta(exon_tmp, FASTA_HT)) == NULL) goto SUCCESS;
		int m = 0;
		/* extending kmer to find MPM */
		while(*(s_fasta->seq + m + _pos_exon) == *(_read+ pos_read + m)){
			m ++;
		}
		if(m >= min_match){
			if(m>max_len){
				max_len=m;
				exon_max = strdup(exon_tmp);
			}
			HASH_FIND_INT(match_lens, &m, s_freq);	
			if(s_freq==NULL){
				if((s_freq=(freq*)malloc(sizeof(freq))) == NULL) goto FAIL_MALLOC;
				s_freq->LEN  = m;
				s_freq->SIZE = 1;				
				HASH_ADD_INT(match_lens, LEN, s_freq);
			}else{
				s_freq->SIZE++;
			}
		}
	}
	/*------------------------------------------------------------*/	
	HASH_FIND_INT(match_lens, &max_len, s_freq);	
	if(s_freq == NULL) goto NO_MATCH; 
	if(s_freq->SIZE==1){
		printf("%s\n", exon_max);
	}
	goto SUCCESS;
	
	FAIL_PARAM:
		ERRBUF("find_next_MEKM: invalid NULL parameter");
		error = PR_ERR_PARAM;
		goto EXIT;
	FAIL_MALLOC:
		ERRBUF("find_next_MEKM: fail to malloc memory");
		error = PR_ERR_MALLOC;
		goto EXIT;
	FAIL_POSPARSE:
		ERRBUF("find_next_MEKM: fail to parse kmer match postion e.g. gene.exon19_101 ");
		error = PR_ERR_POSPARSE;
		goto EXIT;
	FAIL_HASH_FIND:
		ERRBUF("find_next_MEKM: fail to find element in uthahs");
		error = PR_ERR_FINDHASH;
		goto EXIT;	
	FAIL_OTHER:
		ERRBUF("find_next_MEKM: undefined error");
		error = PR_ERR_UNDEFINED;
		goto EXIT;	

	NO_MATCH:
		error = PR_ERR_NONE;
		goto EXIT;
		
	SUCCESS:
		error = PR_ERR_NONE;
		goto EXIT;

	EXIT:		
		if(buff) 		free(buff);
		if(matches) 	free(matches);		
		if(exon_max) 	free(exon_max);	
		if(exon_tmp) 	free(exon_tmp);	
		HASH_ITER(hh, match_lens, s_freq, tmp) {
			HASH_DEL(match_lens, s_freq);
      		free(s_freq);
		}
		return error;
}
/*--------------------------------------------------------------------*/

/* Find all Maximal Extended Kmer Matchs (MEKMs) on _read.            */
size_t 
find_all_MEKMs(char **hits, char* _read, int _k, int min_match){
	int error;
	int _read_pos = 0;
	size_t _i = 0; 
	while(_read_pos<(strlen(_read)-_k-1)){
		char* _exon = NULL;
		find_next_MEKM(_read, _read_pos, _k, min_match);
		//if((error=find_next_MEKM(_read, _read_pos, _k, min_match)) < 0){
		//    fprintf(stderr, "find_all_MEKMs: error=%d\n", error);  
		//	exit(-1);
		//}
		//printf("%d\t_exon=%s\n", _read_pos, _exon);
		//if(_exon) free(_exon);
		_read_pos += 1;
		if (_exon != NULL){
			hits[_i] = _exon;
			_i++;
		}
	}
	return _i;
}

/* Construct breakend associated graph by given fq files.             */
struct BAG_uthash * 
construct_BAG(char *fq_file1, char *fq_file2, int _k, int min_match){
	struct BAG_uthash *graph_ht = NULL;
	gzFile fp1, fp2;
	kseq_t *seq1, *seq2;
	int l1, l2;
	fp1 = gzopen(fq_file1, "r");
	fp2 = gzopen(fq_file2, "r");
	seq1 = kseq_init(fp1);
	seq2 = kseq_init(fp2);
	while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >= 0 ) {
		char *_read1 = rev_com(seq1->seq.s);
		char *_read2 = seq2->seq.s;
		
		if(strcmp(seq1->name.s, seq2->name.s) != 0){
			ERRBUF("construct_BAG: read pair not matched");
			return NULL;
		}			
		
		if(_read1 == NULL || _read2 == NULL)
			continue;    
		
		char **hits1, **hits2;
		if((hits1 = malloc(strlen(_read1) * sizeof(char*)))==NULL){
			ERRBUF("construct_BAG: malloc: no more memory");
			return NULL;
		}
		if((hits2 = malloc(strlen(_read2) * sizeof(char*)))==NULL){
			ERRBUF("construct_BAG: malloc: no more memory");
			return NULL;
		}
		
		if(strlen(_read1) < _k || strlen(_read2) < _k){
			continue;
		}
		size_t num1, num2;		
		num1 = find_all_MEKMs(hits1, _read1, _k, min_match);
		num2 = find_all_MEKMs(hits2, _read2, _k, min_match);
		//// get uniq elements in hits1
		char** hits_uniq1 = malloc(num1 * sizeof(char*));  
		char** hits_uniq2 = malloc(num2 * sizeof(char*));  

		if((hits_uniq1 = malloc(num1 * sizeof(char*)))==NULL){
			ERRBUF("construct_BAG: malloc: no more memory");
			return NULL;
		}
		if((hits_uniq2 = malloc(num2 * sizeof(char*)))==NULL){
			ERRBUF("construct_BAG: malloc: no more memory");
			return NULL;
		}

		size_t size1 = set_str_arr(hits1, hits_uniq1, num1);
		size_t size2 = set_str_arr(hits2, hits_uniq2, num2);
		free(hits1);
		free(hits2);
		int i, j;
		for(i=0; i<size1; i++){
			for(j=0; j<size2; j++){
				char *edge_name;
				size_t len1 = strsplit_size(hits_uniq1[i], ".");
				size_t len2 = strsplit_size(hits_uniq2[j], ".");
				char **parts1 = calloc(len1, sizeof(char *));
				char **parts2 = calloc(len2, sizeof(char *));				
				strsplit(hits_uniq1[i], len1, parts1, ".");  
				strsplit(hits_uniq2[j], len2, parts2, ".");  
				 
				if(parts1[0]==NULL || parts2[0]==NULL){
					free(parts1);
					free(parts2);
					continue;
				}
				
				char* gene1 = strdup(parts1[0]);
				char* gene2 = strdup(parts2[0]);
				
				int rc = strcmp(gene1, gene2);
				if(rc<0)
					edge_name = concat(concat(gene1, "_"), gene2);
				if(rc>0)
					edge_name = concat(concat(gene2, "_"), gene1);
				if(rc==0)
					edge_name = NULL;
				
				if(edge_name!=NULL){
					printf("%s\t%s\n", edge_name, concat(concat(_read1, "_"), _read2));
					BAG_uthash_add(&graph_ht, edge_name, concat(concat(_read1, "_"), _read2));					
				}
				if(edge_name!=NULL){free(edge_name);}
				if(parts1!=NULL){free(parts1);}
				if(parts2!=NULL){free(parts2);}
				if(gene1!=NULL){free(gene1);}
				if(gene2!=NULL){free(gene2);}	
			}
		}
		free(hits_uniq1);
		free(hits_uniq2);
	}
	kseq_destroy(seq1);
	kseq_destroy(seq2);	
	gzclose(fp1);
	gzclose(fp2);
	return graph_ht;
}

/*--------------------------------------------------------------------*/

/* main function. */
int main(int argc, char *argv[]) { 
	if (argc != 5) {  
	        fprintf(stderr, "Usage: %s <in.fa> <read_R1.fq> <read_R2.fq> <int k>\n", argv[0]);  
	        return 1;  
	 }
	char *fasta_file = argv[1];
	char *fq_file1 = argv[2];
	char *fq_file2 = argv[3];
	int min_mtch;		
	if (sscanf (argv[4], "%d", &min_mtch)!=1) { printf ("error - not an integer");}
	/* load kmer hash table in the memory */
	
	int error;
	///* load kmer_uthash table */
	char *index_file = concat(fasta_file, ".index");
	if(index_file == NULL)
		return -1;
	
	//printf("loading kmer uthash table ...\n");
	int k;
	KMER_HT = kmer_uthash_load(index_file, &k);	
	if(KMER_HT == NULL){
		fprintf(stderr, "Fail to load kmer_uthash table\n");
		exit(-1);		
	}
	/* MAX_K is defined in common.h */
	if(k > MAX_K){
		fprintf(stderr, "input k(%d) greater than allowed lenght - 100\n", k);
		exit(-1);		
	}		
	/* load fasta_uthash table */
	FASTA_HT = fasta_uthash_load(fasta_file);
	if(FASTA_HT == NULL){
		fprintf(stderr, "Fail to load fasta_uthash table\n");
		exit(-1);		
	}
	BAG_HT = construct_BAG(fq_file1, fq_file2, k, min_mtch);	
	kmer_uthash_destroy(&KMER_HT);	
	fasta_uthash_destroy(&FASTA_HT);	

	if((error=BAG_uthash_display(BAG_HT))!=0){
		fprintf(stderr, "fails to display BAG_uthash; error=%d\n", error);
		exit(-1);				
	};	
	
	if((error=BAG_uthash_destroy(&BAG_HT))!=0){
		fprintf(stderr, "fails to destory BAG_uthahs with error=%d\n", error);
		exit(-1);		
	}	
	return 0;
}

