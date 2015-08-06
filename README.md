# TFC (Targeted Gene Fusion Calling)

TFC is a lightwieght, stand-alone, ultrafast and highly sensitive tool to detect fusion of candiate genes by Illumina high throughput sequencing reads.

## Version
0.8.05-r15

## Get Started

  - Install

```
$ git clone git@github.com:r3fang/tfc.git
$ cd tfc
$ make
$ ./bin/tfc

Program: tfc (targeted gene fusion calling)
Version: 0.7.30-r15
Contact: Rongxin Fang <r3fang@ucsd.edu>

Usage:   tfc <command> [options]

Command: -seq        extract exon sequences by gene name
         -pred       predict gene fusion
```

  - seq (extract exon sequences of gene name)


```
$./bin/tfc -seq

Usage:   tfc -seq <genes.txt> <genes.gff> <in.fa> <exon.fa>

Inputs:  genes.txt   plain txt file contains names of targeted genes
         genes.gff   gff files contains information of all genes
         in.fa       fasta file contains the entire genome sequence [hg19.fa]
         exon.fa     output fasta files that contains exon sequences of targeted genes
```

  - pred (predict gene fusion)


```
$./bin/tfc -pred

Usage:   tfc -pred [options] <exon.fa> <R1.fq> <R2.fq>

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

## Author
Rongxin Fang (r3fang@eng.ucsd.edu)

