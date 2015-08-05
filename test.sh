#/bin/bash

GENOME="/Users/rongxinfang/Documents/genome/hg19/hg19.fa"
EXON="sample_data/exon.fa"
FQ1="sample_data/reads_simulated_R1.fq"
FQ2="sample_data/reads_simulated_R2.fq"
GENE="sample_data/genes.name.txt"
TFF="sample_data/genes.bed"
TFC="bin/tfc"

# generate exon sequneces if not exist
if [ ! -f $EXON ]
then
	$TFC -seq $GENE $TFF $GENOME $EXON
fi
# predict junction
$TFC -pred $EXON $FQ1 $FQ2