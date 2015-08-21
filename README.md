##Get Started     
```
$ git clone git@github.com:r3fang/tfc.git
$ cd tfc
$ make
$ ./tfc predict exon.fa.gz A431-1-ABGHI_S1_L001_R1_001.fastq.gz A431-1-ABGHI_S1_L001_R2_001.fastq.gz
```

##Introduction

**TFC** is a super *lightweight*, *stand-alone*, *ultrafast*, *C-implemented*, *mapping-free* and *precise* Bioinformatics software desgined for **fast fusion detection** of targeted genes from RNA-seq data. It consists of two major components:
 
```
$ ./tfc 

Program: tfc (targeted fusion calling)
Version: 08.19-r15
Contact: Rongxin Fang <r3fang@ucsd.edu>

Usage:   tfc <command> [options]

Command: name2fasta     extract DNA sequences
         predict        predict gene fusions
```

- **name2fasta** (extract *exon/transcript/CDS* sequences of targeted genes)
 
```
$./tfc name2fasta

Usage:   tfc name2fasta [options] <gname.txt> <genes.gtf> <in.fa> <out.fa> 

Details: name2fasta is to extract genomic sequence of gene candiates

Options: -g               'exon' or 'transcript' or 'CDS' 

Inputs:  gname.txt        .txt file contains the names of gene candiates
         genes.gtf        .gft file contains gene annotations sorted by 5th column
         in.fa            .fa file contains the whole genome sequence e.g. [hg19.fa]
         out.fa           .fa files contains output sequences
```

- **predict** (predict fusions between targeted genes).

```
$ ./tfc predict

Usage:   tfc predict [options] <exon.fa> <R1.fq.gz> <R2.fq.gz>

Details: predict gene fusion from pair-end RNA-seq data

Options:
   -- Graph:
         -k INT    kmer length for indexing in.fa [15]
         -n INT    min unique kmer matches for a hit between gene and pair [10]
         -w INT    edges in graph of weight smaller than -w will be removed [3]
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

Inputs:  exon.fa   .fasta file that contains exon sequences of 
                   targeted genes, which can be generated by: 
                   tfc name2fasta <gname.txt> <genes.gtf> <in.fa> <exon.fa>  
         R1.fq.gz  5'->3' end of pair-end sequencing reads
         R2.fq.gz  the other end of sequencing reads
```
## Workflow

![workflow](https://github.com/r3fang/tfc/blob/master/img/workflow.jpg)

## A Full Example
```
$ sort -k5,5n genes.gtf > genes.sorted.gtf
$ ./tfc name2fasta genes.txt genes.sorted.gtf hg19.fa.gz exon.fa
$ zcat A431-1-ABGHI_S1_L001_R1_001.fastq.gz | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip > A431-1-ABGHI_S1_L001_R1_001.sorted.fastq.gz
$ zcat A431-1-ABGHI_S1_L001_R2_001.fastq.gz | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip > A431-1-ABGHI_S1_L001_R2_001.sorted.fastq.gz
$ ./tfc predict exon.fa A431-1-ABGHI_S1_L001_R1_001.sorted.fastq.gz A431-1-ABGHI_S1_L001_R2_001.sorted.fastq.gz
```
## FAQ

 1. **How fast is TFC?**     
 **~6min** per million read pairs using one CPU core.     
 TFC is 100% implemented in C. We tested TFC on 43 real RNA-seq data with various number of reads ranging from 0.9m to 4m against 506 targeted genes. On average, TFC has ~6min run per million reads for ~500 targeted genes.   
 
 2. **What's the maximum memory requirement for TFC?**   
 **1GB** would be the up limit for most of the cases.   
 The majority (~90%) of the memory occupied by TFC is used for storing the kmer hash table indexed from reference sequences. Thus, the more genes are being tested, the more memory will probably be needed (it also depends on the complexity of the sequences). Based on our simulations, predicting on ~500 genes with k=15 always takes less than **1GB** memory, which means you can definately run TFC on most of today's PC.

 3. **How precise is TFC?**  
 **~0.85** and **~0.99** for sensitivity and specificity on our simulated data.     
 We randomly generated 50 fused transcripts and simulated illumina pair-end sequencing reads from constructed transcripts using [art](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) in paired-end read simulation mode with parameters setting `-l 75 -ss HS25 -f 30 -m 200 -s 10` and run TFC against *paired_reads1.fq* and *paired_reads2.fq* then caculate sensitivity and specificity. Repeat above process for 200 time.

 4. **How is likelihood of fusion being calculated?**   
 In very brief, likelihood equals the product of alignment score of reads that support the fusion normalized by sequencing depth.   
 In detail, let *e(i,j)* indicates the fusion between *gene(i)* and *gene(j)* and *s(i,j)* and *junc(i,j)* be the transcript string and junction site of *e(i,j)*. Let *f(x, y)* be the alignment function between quary string *x* and reference *y*, for any *x* and *y* (*f(x,y)* is always between [0,1]). Let **S(i)** and **S(j)** be the subset of read pairs that aligned to *gene(i)* and *gene(j)* respectively. **S1(i,j)** is the subset of read pairs that support *e(i,j)* and overlapped with *junc(i,j)* and **S2(i,j)** be the subset of read pairs also also support *e(i,j)* but not overlaped with *junc(i,j)*. Liklihood of *e(i,j)* can be calculated by      
 ![equation](https://github.com/r3fang/tfc/blob/master/img/Tex2Img_1440195992.jpg)    
 in which ![equation](https://github.com/r3fang/tfc/blob/master/img/Tex2Img_1440196064.jpg)

 5. **What's the null model for calculating p-value?**   
 In brief, p-value is the probability of observing the likelihood by null model.   
 We extracted all transcripts of targeted genes and simulate pair-end reads by art. Then run `./tfc predict` against simulated data and calculate likelihood for all gene pairs. Repeat this for 200 times and get the distribution of likelihood of every gene pair. 

 6. **How does TFC guarantee specificity without comparing sequencing reads against regions outside targeted genes?**   
 we have several strict criteria to filter out read pairs that are likely to come from regions outside targeted loci. For instance, both ends of a pair are aligned to the constructed transcript and those pairs of any end not being aligned with score smaller than `-a FLOAT min identity score for alignment` will be discarded. Also, any pair with too large or too small insertion size will be filtered out. 

 7. **Does tfc work for single-end reads?**   
 Unfornately, TFC only works for pair-end sequencing data now, but having it run for single-end read is a feature we would love to add in the near future.

 8. **Does TFC support parallel computing?**    
 No. We realize TFC is fast enough, but this is a feature we would love to add in the near future.

 9.  **Is there anything I should be very careful about for `./tfc name2fasta`?**    
 Yes, genes.gtf needs to be sorted by its 5th column as shown above. 

 10. **Is there anything I should be very careful about for `./tfc predict`?**  
 3 things.    

- First, exon.fa has to be in the following format, in which *SORT1.1* indicates this is the first exon of gene *SORT1*. exon.fa can be generated by **name2fasta**     
 *\>SORT1.1	chr1|109852188	strand	- gene_id	SORT1	transcript_id	NM_002959|NM_001205228|	tss_id	TSS12777|TSS22486|*     
 *ATCCAGTT...TTAACACAC*    
 *\>SORT1.2        chr1|109856883  strand  - gene_id SORT1   transcript_id   NM_002959|NM_001205228| tss_id  TSS12777|TSS22486|*  
 *TACACAC...TTTTTTTTTAA*       
- Second, when you run `tfc predict [options] <exon.fa> <R1.fq> <R2.fq>`, R1.fq and R2.fq (RNA-seq) must be in the right order that R2.fq must be identical to the psoitive strand of reference genome.         
- Third, name of reads has to be paired up in R1.fq and R2.fq, sort them based on read name if necessary.

####Author   

Rongxin Fang
r3fang@eng.ucsd.edu
