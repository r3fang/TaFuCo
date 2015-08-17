/*--------------------------------------------------------------------*/
/* Created Date: 15JULY2015                                           */
/* Author: Rongxin Fang                                               */
/* Contact: r3fang@ucsd.edu                                           */
/* Extract exon sequences of candiate genes.                          */
/*--------------------------------------------------------------------*/

#ifndef _NAME2FASTA_H
#define _NAME2FASTA_H

#include <stdio.h>  
#include <stdlib.h>  
#include <string.h> 
#include <zlib.h>  
#include <math.h>
#include "kseq.h"
#include "kstring.h"
#include "fasta_uthash.h"
#include "utils.h"

/*
 * Description:
 *------------
 * Extract exon sequences by providing gene's name

 * Input: 
 *-------
 * fname     - name of file that contains genes' name e.g. genes.name.txt
 * fname_db  - name of the file that contains all genes' annotation e.g. data/hg.bed
 * HG19_HT   - fasta_uthash object that contains reference genome loaded by fasta_uthash_load

 * Output: 
 *-------
 * fasta_uthash object that contains extracted sequences.
 */
static fasta_t *extract_exon_seq(char* fname, char *fname_db, fasta_t *HG19_HT);

/*
 * usage info
 */
int name2fasta_usage();

/*
 * main function
 */
int name2fasta(int argc, char *argv[]);

#endif