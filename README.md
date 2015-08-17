##Get Started

```
$ git clone git@github.com:r3fang/tfc.git
$ cd tfc
$ make
$ ./tfc predict exon.fa A431-1-ABGHI_S1_L001_R1_001.fastq.gz A431-1-ABGHI_S1_L001_R2_001.fastq.gz
```


##Introduction

TFC is a super lightwieght, stand-alone, ultrafast, C-implemented, mapping-free and highly sensitive and precise software desgined for fusion detection between candiate genes by RNA-seq data. It consists of two major components: 
 
```
$ ./tfc

Program: tfc (targeted gene fusion calling)
Version: 08.17-r15
Contact: Rongxin Fang <r3fang@ucsd.edu>

Usage:   tfc <command> [options]

Command: name2fasta     extract exon sequences
         predict        predict gene fusion
```

 - name2fasta

```
$./tfc name2fasta

Usage:   tfc name2fasta [options] <gname.txt> <in.fa.gz> <out.fa> 

Details: name2fasta is to extract genomic sequence of gene candiates

Options: -g          organism - 0 for human; 1 for mouse

Inputs:  gname.txt   plain txt file contains names of gene candiates e.g. [genes.txt]
         hg19.fa     fasta file contains the whole genome sequence   e.g. [hg19.fa]
         out.fa      output fasta files contains extracted seq of targeted genes
```

  - predict
	
```
$ ./tfc predict

Usage:   tfc predict [options] <exon.fa> <R1.fq> <R2.fq>

Details: predict gene fusion from RNA-seq data

Options: -k INT    kmer length for indexing genome [15]
         -n INT    min number of uniq kmer matches for a gene-read hit [10]
         -w INT    min weight for an edge in BAG [3]
         
		 -m INT    score for match in alignment [2]
         -u INT    score for mismatch in alignment [-2]
         -o INT    penality for gap open [-5]
         -e INT    penality for gap extension [-1]
         -j INT    penality for jump between genes [-10]
         -s INT    penality for jump between exons [-8]
         -a FLOAT  min identity for an alignment [0.80]
         
		 -h INT    min unique read hits for a junction [3]
		 -l INT    length for junction string [20]         
		 -x INT    max mismatches of junction string match [2]

		 -A INT    alpha
		 -B INT    beta
		 
Inputs:  exon.fa   fasta file that contains exon sequences of targeted 
                   genes with no flanking sequence which can be generated: 
                   tfc name2fasta <genes.txt> <in.fa> <exon.fa> 
         R1.fq     5'->3' end of pair-end sequencing reads
         R2.fq     the other end of pair-end sequencing reads
```

## Workflow

![workflow](https://github.com/r3fang/tfc/blob/master/img/workflow.jpg)

## FAQ

1. How fast is tfc?
2. Does tfc need reads to be mapped in advance?
3. How precise is tfc?
4. Does tfc work for single-end reads?

#### Version
08.17-r15

#### Author
Rongxin Fang (r3fang@eng.ucsd.edu)

