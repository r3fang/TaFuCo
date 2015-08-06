##Get Started

```
$ git clone git@github.com:r3fang/tfc.git
$ cd tfc
$ make
$ ./tfc -predict exon.fa reads1.fq.gz reads2.fq.gz > tfc.out
```


##Introduction

TFC is a lightwieght, stand-alone, ultrafast, C-implemented, mapping-free and highly sensitive software desgined for detection of fusion between candiate genes using Illumina RNA-seq data. It consists of two major components: 
 
 - name2fasta is to extract genomic sequences of genes candiates that only requires user to provide the names of targeted genes.

```
$./tfc name2fasta

Usage:   tfc name2fasta <genes.txt> <genes.gff> <in.fa> <exon.fa>

Inputs:  genes.txt   plain txt file contains names of targeted genes
         genes.gff   gff files contains information of all genes
         in.fa       fasta file contains the entire genome sequence [hg19.fa]
         exon.fa     output fasta files that contains exon sequences of targeted genes
```

  - predict is to predict gene fusion from RNA-seq data.
	
```
$./tfc predict

Usage:   tfc predict [options] <exon.fa> <R1.fq> <R2.fq>

Options: -k INT    kmer length [15]
         -n INT    min number kmer matches [10]
         -w INT    min weight for an edge [3]
         -m INT    match score [2]
         -u INT    penality for mismatch [-2]
         -o INT    penality for gap open [-5]
         -e INT    penality for gap extension [-1]
         -j INT    penality for jump between genes [-10]
         -s INT    penality for jump between exons [-8]
         -h INT    min number of hits for a junction to be called [3]
         -l INT    length of exact match seed for junction rediscovery [20]
         -x INT    max number of mismatches of seed exact match for junction rediscovery [2]
         -a FLOAT  min alignment score [0.80]

Inputs:  exon.fa   fasta file that contains exon sequences of targeted 
                   genes with no flanking sequence which can be generated: 
                   tfc -seq <genes.txt> <genes.gff> <in.fa> <exon.fa> 
         R1.fq     5'->3' end of pair-end sequencing reads
         R2.fq     the other end of pair-end sequencing reads
```

#### Version
0.8.05-r15


#### Author
Rongxin Fang (r3fang@eng.ucsd.edu)

