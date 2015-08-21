##Get Started

```
$ git clone git@github.com:r3fang/tfc.git
$ cd tfc
$ make
$ ./tfc predict exon.fa.gz A431-1-ABGHI_S1_L001_R1_001.fastq.gz A431-1-ABGHI_S1_L001_R2_001.fastq.gz
```

##Introduction

**TFC** is a super *lightweight*, *stand-alone*, *ultrafast*, *C-implemented*, *mapping-free* and *sensitive* Bioinformatics software for **fusion detection** between candidate genes from RNA-seq data. It consists of two major components:
 
```
$ ./tfc 

Program: tfc (targeted fusion calling)
Version: 08.19-r15
Contact: Rongxin Fang <r3fang@ucsd.edu>

Usage:   tfc <command> [options]

Command: name2fasta     extract DNA sequences
         predict        predict gene fusions
```

- **name2fasta** 
  
> extract *exon/transcript/CDS* sequences of targeted genes. Before running it, .gtf file has to be sorted based on the 4th column by  
`sort -k5,5n genes.gtf > genes.sorted.gtf`
 
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

- **predict** 
  
> predict fusions between targeted genes. **IMPORTANT** before running it, you need to make sure R1.fq and R2.fq have their read's name matched up. Sort R1.fq and R2.fq based on id if necessary 
```
$ zcat R1.fq.gz | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip > R1.sorted.fq.gz
$ zcat R2.fq.gz | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip > R2.sorted.fq.gz
```

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

## FAQ

 1. **How fast is TFC?**     
 On average, ~6min for 1 million read pairs.     
 TFC is 100% implemented in C. We tested TCF on 43 real RNA-seq data with various number of reads ranging from 0.9m to 4m against 506 targeted genes. On average, TFC has ~6min run per million reads.   
 
  |Sample         | Reads Number   | Running Time |
  |:-------------:| :-------------:| :-------------:|
  |1       | 10M            | 6min         |
  |2  | 5M             | 5min         |
  |3              | 0.4M           | 5min         |
  |4              | 5M             | 5min         |
  |...            | ...            | ...          |
  |41             | 5M             | 5min         |
  |42             | 5M             | 5min         |
  |43             | 5M             | 5min         |  
  |**average**        | **1.7M**           | **8min**         |  
  
 2. **What's the maximum memory requirement for TFC?**   
 **1GB** would be the up limit for most of the cases.   
 The majority (~90%) of the memory occupied by TFC is used for storing the kmer hash table indexed from reference sequences. Thus, the more genes are being tested, the more memory will probably be needed (it also depends on the complexity of the sequences). Based on our simulations, predicting on ~500 genes with k=15 always takes less than **1GB** memory, which means you can definately run TFC on most of today's PC.
 3. **Does TFC depend on any third-party software?**   
 No. TFC is compeletely stand-alone.
 4. **How precise is TFC?**  
 **0.85+-0.04** for sensitivity and **0.99+-0.005** for specificity based on our simulations.     
 We randomly generated 50 fused transcripts and simulated illumina pair-end sequencing reads from fused transcripts using [art](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) in paired-end read simulation mode with parameters `-l 75 -ss HS25 -f 30 -m 200 -s 10` and run TFC on *paired_reads1.fq* and *paired_reads2.fq* then caculate Sensitivity and Specificity. Repeat above process for 100 time. 
 5. **How does TFC guarantee specificity without comparing sequencing reads against regions outside targeted genes?**   
 we have several strict criteria to filter out read pairs that are likely to come from regions outside targeted loci. For instance, both ends of a pair are aligned to the constructed transcript and those pairs of any end not being aligned with a fair score will be discarded. Also, any pair with too large or too small insertion size will be filtered out. 
 6. **Does tfc work for single-end reads?**   
 Unfornately, TFC only works for pair-end sequencing data now, but having it run for single-end read is a feature we would love to add in the near future.
 7. **Is there anything I should be very careful about for `./tfc predict`?**  
 3 things.    
 
 - First, exon.fa has to be in the following format, in which *SORT1.1* indicates this is the first exon of gene *SORT1*. exon.fa can be generated by **name2fasta**     
 *\>SORT1.1	chr1|109852188	strand	- gene_id	SORT1	transcript_id	NM_002959|NM_001205228|	tss_id	TSS12777|TSS22486|*     
 *ATCCAGTT...TTAACACAC*    
 *\>SORT1.2        chr1|109856883  strand  - gene_id SORT1   transcript_id   NM_002959|NM_001205228| tss_id  TSS12777|TSS22486|*  
 *TACACAC...TTTTTTTTTAA*       
 - Second, When you run `tfc predict [options] <exon.fa> <R1.fq> <R2.fq>`, R1.fq and R2.fq (RNA-seq) must be in the right order that R2.fq must be identical to the psoitive strand of reference genome.      
 - Third, name of reads has to be paired up in R1.fq and R2.fq, sort them based on read name if necessary.           
 
 8. **Does TFC support parallel computing?**    
 No. We realize TFC is fast enough, but this is a feature we would love to add in the near future. 

#### Version
08.19-r15

#### Author
Rongxin Fang (r3fang@eng.ucsd.edu)

