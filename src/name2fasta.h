#ifndef _NAME2FASTA_H
#define _NAME2FASTA_H

#include <stdio.h>  
#include <stdlib.h>  
#include <string.h> 
#include <zlib.h>  
#include <assert.h>
#include <math.h>
#include <regex.h>
#include "kseq.h"
#include "fasta_uthash.h"
#include "utils.h"

struct fasta_uthash *extract_exon_seq(char*, char *, struct fasta_uthash *);
int name2fasta_usage();
int name2fasta(int argc, char *argv[]);

#endif