/*--------------------------------------------------------------------*/
/* fit_affine_jump.c	                                              */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Date: 07-23-2015                                                   */
/* Pair wise fit alignment with affine gap and jump state.            */
/* This could be used to align RNA-seq reads with intron splicing as  */
/* jump state.                                                        */
/*                                                                    */
/* initilize L(i,j), M(i,j), U(i,j), J(i,j), G1(i,j), G2(i,j):        */
/*--------------------------------------------------------------------*/
/* M(0,0)  = -INF                                                     */
/* M(i,0)  = -INF               (i>0 )                                */
/* M(0,j)  = 0                  (j>0 )                                */
/* L(0,j)  = -INF               (j>=0)                                */
/* L(i,0)  = -INF               (i>=0)                                */
/* U(i,0)  = -INF               (i>0 )                                */
/* U(0,j)  = 0                  (j>=0)                                */
/* G1(0,j) = -INF               (j>=0)                                */
/* G1(i,0) = -INF               (i>=0)                                */
/* G2(0,j) = -INF               (j>=0)                                */
/* G2(i,0) = -INF               (i>=0)                                */
/* J(0,j)  = -INF               (j>=0)                                */
/* J(i,0)  = -INF               (i>=0)                                */
/*                                                                    */
/* reccurrance relations:                                             */
/*--------------------------------------------------------------------*/
/*  M(i,j) = max{M(i-1, j-1) + s(x,y), U(i-1, j-1) + s(x,y),          */
/*              L(i-1, j-1) + s(x,y),                                 */
/*              J(i-1, j-1) + s(x,y) (j>x),                           */
/*              G1(i-1, j-1) + s(x,y) if j in S1,                     */
/*              G2(i-1, j-1) + s(x,y) if j in S2}                     */                                          
/*  U(i,j) = max{M(i-1, j) + GAP, U(i-1, j) + EXTENSION}              */
/*  L(i,j) = max{M(i, j-1) + GAP, L(i, j-1) + EXTENSION}              */
/*--------------------------------------------------------------------*/
/* G1(i,j) = max{M(i-1, j) + JUMP_EXON (j in S1), G1(i-1, j)}        */
/* G2(i,j) = max{M(i-1, j) + JUMP_EXON (j in S2), G2(i-1, j)}        */
/*  J(i,j) = max{M(i, j-1) + JUMP_GENE (j < x), J(i-1, j)}            */

/* Traceback:                                                         */
/*--------------------------------------------------------------------*/
/* start at max(M(m,j_max), L(n, j_max)), Stop at any of i=0 on M/L;  */
/* The rational behind this is no gap allowed to flank s1             */
/*--------------------------------------------------------------------*/
/* NOTE:                                                              */
/* 1. WE ONLY ALLOW STATE CHANGE FROM M(MATCH) TO J(JUMP) AT SPECIFIC */
/* POSITIONS ON S2.                                                   */ 
/* 2. WE ALLOW STATE CHANGE FROM JUMP TO MATCH AT ANYWHERE            */
/* 3. THE READ SHOULD NOT START WITH JUMP AND ALSO SHOULD NOT END     */
/* WITH JUMP STATE.                                                   */ 
/*--------------------------------------------------------------------*/

#include "alignment.h"
int S1[6] = {130, 131, 132, 407, 408, 409};
int S2[6] = {761, 762, 763, 844, 845, 846};
int JUNCTION = 612;

/*
 * fit alignment with affine gap penality
 */
double align(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2, opt_t *opt){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL || opt == NULL) die("align: parameter error\n");
	if(s1->l > s2->l) die("first sequence must be shorter than the second to do fitting alignment"); 
	size_t m   = s1->l + 1; size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	// copy alignment parameter
	double match = opt->m;
	double mismatch = opt->u;
	double gap = opt->o;
	double extension = opt->e;
	double jump_penality = opt->j;
	// initlize leftmost column
	int i, j;
	for(i=0; i<S->m; i++){
		S->M[i][0] = -INFINITY;
		S->U[i][0] = -INFINITY;
		S->L[i][0] = -INFINITY;
		S->J[i][0] = -INFINITY;
		S->G1[i][0] = -INFINITY;
		S->G2[i][0] = -INFINITY;
	}
	// initlize first row
	for(j=0; j<S->n; j++){
		S->M[0][j] = 0.0;
		S->U[0][j] = 0.0;
		S->L[0][j] = -INFINITY;
		S->J[0][j] = -INFINITY;
		S->G1[0][j] = -INFINITY;
		S->G2[0][j] = -INFINITY;
	}
	double delta, tmp_J, tmp_G1, tmp_G2, tmp_M;
	int idx;
	
	// recurrance relation
	for(i=1; i<=s1->l; i++){
		for(j=1; j<=s2->l; j++){
			// MID any state can goto MID
			delta = ((s1->s[i-1] - s2->s[j-1]) == 0) ? MATCH : MISMATCH;
			tmp_J = (j > JUNCTION) ?  S->J[i-1][j-1]+delta : -INFINITY;
			tmp_G1 = (isvalueinarray(j, S1, 6)) ?  S->G1[i-1][j-1]+delta : -INFINITY;
			tmp_G2 = (isvalueinarray(j, S2, 6)) ?  S->G2[i-1][j-1]+delta : -INFINITY;
			idx = max6(&S->M[i][j], S->L[i-1][j-1]+delta, S->M[i-1][j-1]+delta, S->U[i-1][j-1]+delta, tmp_J, tmp_G1, tmp_G2);
			if(idx == 0) S->pointerM[i][j]=LOW;
			if(idx == 1) S->pointerM[i][j]=MID;
			if(idx == 2) S->pointerM[i][j]=UPP;
			if(idx == 3) S->pointerM[i][j]=JUMP;			 				
			if(idx == 4) S->pointerM[i][j]=GENE1;			 				
			if(idx == 5) S->pointerM[i][j]=GENE2;			 				
			// LOW
			idx = max6(&S->L[i][j], S->L[i][j-1]+EXTENSION, S->M[i][j-1]+GAP, -INFINITY, -INFINITY, -INFINITY,  -INFINITY);
			if(idx == 0) S->pointerL[i][j]=LOW;
			if(idx == 1) S->pointerL[i][j]=MID;
			
			// UPP
			idx = max6(&S->U[i][j], S->M[i-1][j]+gap, S->U[i-1][j]+extension,  -INFINITY, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerU[i][j]=MID;
			if(idx == 1) S->pointerU[i][j]=UPP;
			// JUMP 
			tmp_M = (j < JUNCTION) ?  S->M[i][j-1]+JUMP_GENE : -INFINITY;
			idx = max6(&S->J[i][j], tmp_M, S->J[i][j-1], -INFINITY,  -INFINITY,  -INFINITY, -INFINITY);
			if(idx == 0) S->pointerJ[i][j] = MID;			
			if(idx == 1) S->pointerJ[i][j] = JUMP;			
			// GENE1
			tmp_M = (isvalueinarray(j, S1, 6)) ?  S->M[i][j-1]+JUMP_EXON : -INFINITY;
			idx = max6(&S->G1[i][j], tmp_M, S->G1[i][j-1], -INFINITY, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerG1[i][j] = MID;			
			if(idx == 1) S->pointerG1[i][j] = GENE1;
			// GENE2
			tmp_M = (isvalueinarray(j, S2, 6)) ?  S->M[i][j-1]+JUMP_EXON : -INFINITY;
			idx = max6(&S->G2[i][j], tmp_M, S->G2[i][j-1], -INFINITY, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerG2[i][j] = MID;			
			if(idx == 1) S->pointerG2[i][j] = GENE2;
			}
		}
	
	// find trace-back start point
	// NOTE: ALWAYS STARTS TRACING BACK FROM MID OR LOW
	int i_max, j_max;
	double max_score = -INFINITY;
	int max_state;
	i_max = s1->l;
	for(j=0; j<s2->l; j++){
		if(max_score < S->M[i_max][j]){
			max_score = S->M[i_max][j];
			j_max = j;
			max_state = MID;
		}
	}
	for(j=0; j<s2->l; j++){
		if(max_score < S->L[i_max][j]){
			max_score = S->L[i_max][j];
			j_max = j;
			max_state = LOW;
		}
	}
	//trace_back(S, s1, s2, r1, r2, max_state, i_max, j_max);	
	destory_matrix(S);
	return max_score;
}


/* main function. */
int main(int argc, char *argv[]) {
	opt_t *opt = init_opt(); // initlize options with default settings
	int c;
	srand48(11);
	while ((c = getopt(argc, argv, "m:u:o:e:j:s")) >= 0) {
			switch (c) {
			case 'm': opt->m = atoi(optarg); break;
			case 'u': opt->u = atoi(optarg); break;
			case 'o': opt->o = atoi(optarg); break;
			case 'e': opt->e = atoi(optarg); break;
			case 'j': opt->j = atoi(optarg); break;
			case 's': opt->s = true; break;
			default: return 1;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
				fprintf(stderr, "Usage:   alignTools fit [options] <target.fa>\n\n");
				fprintf(stderr, "Options: -m INT   score for a match [%d]\n", opt->m);
				fprintf(stderr, "         -u INT   mismatch penalty [%d]\n", opt->u);
				fprintf(stderr, "         -o INT   gap open penalty [%d]\n", opt->o);
				fprintf(stderr, "         -e INT   gap extension penalty [%d]\n", opt->e);
				fprintf(stderr, "         -j INT   jump penality [%d]\n", opt->j);
				fprintf(stderr, "         -s       weather jump state include\n");
				fprintf(stderr, "\n");
				return 1;
	}
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	printf("%s\n", argv[argc-1]);
	ks1->s = "ACTGTACAACTGAAGGACTGACATGGCAATCCTTAAGAATTTTACCTACAGAATGAATGCACACATATATCATTCCAAATTTGTCAACTTATATAAATATGGTCCCATTTTCACAGTTAATTGGCTCAGTACAGTCACCCATAAGCCACTGTTCTTTAAACAGGAAGCTAAATTTATTTACATAGGAAGCTGCAATTTATTCCCATCTCAATATGGCAATAAATGGAGATATACAAGCTTTGTAAAGAAAAAGAGATCCATACGTACTGGGGAAATACTCATGTGTTTCATGTATTTTTCAAAAGAATTATATTCATAAGGAAACCA";
	ks2->s = "ATCCACTTTTACAATATGTAAAAGGTACTTTTAACTTCCTTTCATTGAACCAGTGTACAACAGTTCACTGTACAACTGAAGGACTGACATGGCAATCCTTAAGAATTTTACCTACAGAATGAATGCACACAatttgacacattttcttagtttcaaaagattatttaaaaaaggaattcagtagattgacttgtaaataaccattgcagattttgaatctgcaaaaatccgtcacattgctgttgggacagattaagataaggctaaaatttttttccaagttattgcaaactgctatggaaagagaaagtaatccaaaatgtataaaatggcccatggacaaatccaaaccacgcaatttttgtaaataaaggtttattgcaatatggccacatctacttactcatgTATATCATTCCAAATTTGTCAACTTATATAAATATGGTCCCATTTTCACAGTTAATTGGCTTCACCAAGTAAGAAAATATGGGTAAAAACACAATTCAAGGTCACTCAAGTTTATCATCCTCGTAAGTAACAACAGCTCTCTATTTGAAGGTATATGGGAATCTCAAGTAGAATATTCAAGACTTTCTTAACAATATGTAAATAGCTCAATAATAGCAAACTTTAAAATGTACATTTTTCTCAGTACAGTCACCCATAAGCCACTGTTCTTTAAACAGGAAGCTAAATTTATTTACATAGGAAGCTGCAATTTATTCCCATCTCAATATGGCAATAAATGGAGATATACAAGTAttcgctatattaggtaagttatatccgataggacaggtttacatttacagttactcaattcatgttaataagcttttttccagCCATGTTTCAACTACTTTGTAAAGAAAAAGAGATCCATACGTACTGGGGAAATACTCATGTGTTTCATGTATTTTTCAAAAGAATTATATTCATAAGGAAACCACCACTTAATTATTCCAATCTCATGCAGGCTCAATGATTTTTTAAGTCAGTATTTTGTCCTATGGAACTAATCCAACCGAGTAACAGATGATTTTAAAAATACGCTCATTTAAATACAAACA";
	ks1->l = strlen(ks1->s);
	ks2->l = strlen(ks2->s);
	
	//kstring_read(argv[argc-1], ks1, ks2, opt);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	if(ks1->l > ks2->l) die("first sequence must be shorter than the second\n");
	kstring_t *r1 = mycalloc(1, kstring_t);
	kstring_t *r2 = mycalloc(1, kstring_t);
	r1->s = mycalloc(ks1->l + ks2->l, char);
	r2->s = mycalloc(ks1->l + ks2->l, char);
	printf("score=%f\n", align(ks1, ks2, r1, r2, opt));
	//printf("%s\n%s\n", r1->s, r2->s);
	//kstring_destory(ks1);
	//kstring_destory(ks2);
	//kstring_destory(r1);
	//kstring_destory(r2);
	//free(opt);
	return 0;
}
