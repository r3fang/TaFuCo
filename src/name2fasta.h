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
 * usage info
 */
int name2fasta_usage();

/*
 * main function
 */
int name2fasta(int argc, char *argv[]);

#endif