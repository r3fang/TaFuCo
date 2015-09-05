[<img src="https://github.com/r3fang/TaFuCo/blob/master/img/sp5.jpg" width="110px" height="30px">](https://www.illumina.com/)
##Get Started     
```
$ git clone https://github.com/r3fang/TaFuCo.git
$ cd TaFuCo
$ make
$ ./TaFuCo rapid A431-1-ABGHI_S1_L001_R1_001.fastq.gz A431-1-ABGHI_S1_L001_R2_001.fastq.gz
```

##Introduction

**TaFuCo** is a precise, fast, C-implemented, lightweight, stand-alone and mapping-free Bioinformatics software designed for **targeted fusion detection** from RNA-seq data. TaFuCo has two modes, **rapid** and **predict**, **rapid** allows user to have a quick run against a list of predefined gene candidates with default parameter settings. **predict** is more flexible and allows user to decide the gene candidates and parameters, but it needs more input files from the user (e.g. hg19.fa, genes.gtf). In brief, **rapid** is easier to use and **predict** is more flexible.

```
$ ./TaFuCo 

Program: TaFuCo (targeted gene fusion calling)
Version: 08.29-r15
Contact: Rongxin Fang <r3fang@ucsd.edu>

Usage:   TaFuCo <command> [options]

Command: rapid          predict gene fusions in rapid mode
         predict        predict gene fusions in predict mode
         name2fasta     extract DNA sequences of targeted genes
```

- **rapid** (predict fusions in rapid mode)

```
$ ./TaFuCo rapid

Usage:   TaFuCo rapid <R1.fq> <R2.fq>

Details: predict fusions in a rapid mode

Inputs:  R1.fq     5'->3' end of pair-end sequencing reads
         R2.fq     the other end of sequencing reads
```

- **predict** (predict fusions in predict mode).

```
$ ./TaFuCo predict

Usage:   TaFuCo predict [options] <gname.txt> <genes.gtf> <in.fa> <R1.fq> <R2.fq>

Details: predict gene fusion from pair-end RNA-seq data

Options:
   -- Graph:
         -k INT    kmer length for indexing in.fa [15]
         -n INT    min unique kmer matches for a hit between gene and pair [10]
         -w INT    edges in graph of weight smaller than -w will be removed [4]
   -- Alignment:
         -m INT    score for match [2]
         -u INT    penality for mismatch[-2]
         -o INT    penality for gap open [-5]
         -e INT    penality for gap extension [-1]
         -j INT    penality for jump between genes [-10]
         -s INT    penality for jump between exons [-8]
         -a FLOAT  min identity score for alignment [0.80]
   -- Junction:
         -h INT    min hits for a junction [3]
         -l INT    length for junction string [20]
         -x INT    max mismatches allowed for junction string match [2]
   -- Fusion:
         -A INT    weight for junction containing reads [3]
         -p FLOAT  p-value cutoff for fusions [0.05]

Inputs:  gname.txt plain txt file that contains name of gene candidates
         genes.gtf gtf file that contains gene annotation
         in.fa     fasta file that contains reference genome
         R1.fq     5'->3' end of pair-end sequencing reads
         R2.fq     the other end of sequencing reads
```
## Workflow

![workflow](https://github.com/r3fang/TaFuCo/blob/master/img/workflow.jpg)

### A Full Example for Rapid Mode
```
$ ./TaFuCo rapid A431-1-ABGHI_S1_L001_R1_001.fastq.gz A431-1-ABGHI_S1_L001_R2_001.fastq.gz
```
### A Full Example for Predict Mode
```
$ sort -k5,5n genes.gtf > genes.sorted.gtf
$ ./TaFuCo predict data/genes.txt data/genes.sorted.gtf hg19.fa A431-1-ABGHI_S1_L001_R1_001.fastq.gz A431-1-ABGHI_S1_L001_R2_001.fastq.gz
```

## FAQ

 1. **How fast is TaFuCo?**     
 On average, **~5min** per million pairs using a single x86_64 32-bit 2000 MHz GenuineIntel processor.   
 We tested TaFuCo (rapid mode) on 43 real RNA-seq data against 506 genes candidates. On average, TaFuCo spends ~5min per million pairs. However, the running time is not absolutely linear to the number of reads. We found most of the time has been spent on the alignment for the step *fusion refinement* and *junction refinement*, therefore, the more fusions in the sample identified, the longer TaFuCo usually runs. 

 2. **What's the maximum memory requirement for TaFuCo?**   
 **1GB** would be enough for **rapid** mode predicting against 1,000 genes, **predict** will take over more memory since it reads the whole reference genome into RAM.     
 The majority (~90%) of the memory occupied by TaFuCo-rapid is used for storing the kmer hash table indexed from reference sequences. Thus, the more genes are being tested, theoretically the more memory will be taken over. Based on our simulations, predicting against ~1000 genes with k=15 always takes less than **1GB** memory in rapid mode, which means TaFuCo can safely be used on most of today's PCs.  

 3. **How precise is TaFuCo?**  
 **~0.85** and **~0.99** for sensitivity and specificity on the simulated data.     
 We randomly constructed 50 fused transcripts and simulated Illumina pair-end sequencing reads from constructed transcripts using [art](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) in paired-end mode with parameters setting as `-p -l 75 -ss HS25 -f 30 -m 200 -s 10` and run ```./TaFuCo rapid``` against simulated reads, then calculate sensitivity and specificity. Repeat above process for 200 time and get the average sensitivity and specificity.

 4. **What are the genes being tested in the rapid mode?**  
 We have a list of 506 genes that have been found fused with other genes in the literature, those genes serve as a default gene set if user does not really know what genes should be tested. 

 5. **How many genes should I test each time?**  
 Let's make it clear, the more genes being tested, the more memory TaFuCo will grasp. The relationship between number of tested genes (N) and memory usage is linear(R=0.99) when N<3000 as shown in the table below.
 
 | # of tested genes  | vmsize(GB) | vmpeak(GB) | vmrss(GB) | vmhwm(GB) |   
 |:------------------:|:-------------:|:-----:|:------:|:-------:|
 | 500  | 0.29	| 0.29	|0.28	| 0.28 |
 | 1000 | 0.57	| 0.57	|0.56	| 0.56 |
 | 1500 | 0.89	| 0.89	|0.88	| 0.88 |
 | 2000 | 1.19	| 1.19	|1.19	| 1.19 |
 | 2500 | 1.54	| 1.59	|1.53	| 1.58 |
 | 3000 | 1.89	| 1.89	|1.88	| 1.88 |
 
 6. **How is the score of fusion calculated?**   
 Let Si, Sj be the set of reads that match with gene_i and gene_j respectively. S1_ij is a set of reads that support fusion e_ij and overlap with its junction. S2_ij is a set of reads that support fusion e_ij but not overlapping with its junction. f(s) indicates the alignment score of s against the real transcript of the fusion. Then score equals the product of alignment score of reads that support the fusion normalized by sequencing depth.
<p align="center">
  <img src="https://github.com/r3fang/TaFuCo/blob/master/img/likelihood.jpg" width="400px" height="250px">
</p>
 
 7. **What's the null model for p-value?**   
 We extracted normal transcripts of targeted genes and simulated pair-end reads from the normal transcripts. Then run TaFuCo against the simulated data and calculate the score (defined above) for every gene pair. Repeat this for 200 times and get the distribution of score of every gene pair as the null model. 

 8. **How does TaFuCo guarantee specificity without comparing sequencing reads against regions outside targeted genes?**   
 we have several strict criteria to filter out reads that are likely to come from regions outside targeted intervals. For instance, both ends of a pair are aligned against the constructed transcript and those pairs of any end not being successfully aligned will be discarded. Also, any pair with too large or too small insertion size will be filtered out. The likelihood of fusion will be normalized by sequencing depth of the two genes before p-value is calculated.

 9.  **Is there anything I should be very careful about for `./TaFuCo name2fasta`?**    
 Yes, genes.gtf needs to be sorted by its 5th column, the end position of the feature.

 10. **Is there anything I should be very careful about for `./TaFuCo predict`?**  
 Yes, before running `TaFuCo predict [options] <exon.fa> <R1.fq> <R2.fq>`, user has to make sure R1.fq and R2.fq (RNA-seq) are in the right order that R2.fq must be identical to the psoitive strand of reference genome.         

####Version     
09.04-r15

####Author     
Rongxin Fang    
r3fang@eng.ucsd.edu
