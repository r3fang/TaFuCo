#!/bin/bash
# 07-06-2015
# Rongxin Fang

# FILEs
GENOME="/bioinfoSD/users/rfang/data/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.7.x/genome.fa"
GTF="/bioinfoSD/users/rfang/data/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
EXONS_NAME="exons"
GENE_NAME="genes.name.txt"

cd /bioinfoSD/users/rfang/tfc/sample_data
	
# Step 1. Extract the exons by given genes' name
#if [ -e $EXONS_NAME.bed ]
#then
#  echo "$EXONS_NAME.bed already exists, delete it first."
#  exit
#fi
################################################################
while read name
do
	grep "$name" $GTF |\
	awk -v myvar2=$name -v myvar=\"$name\"\; '{if($3 == "exon" && $12 == myvar) printf("%s\t%s\t%s\t%s\t.\t%s\n", $1, $4, $5, myvar2, $7)}' |\
	sort -k 1,1 -k2,2n -k3 - | uniq |\
	python ../bin/extract_exon.py - >> $EXONS_NAME.bed
done < $GENE_NAME
#
# Step 2. Extract exon sequences
 
 fastaFromBed -s -name -fi $GENOME -bed $EXONS_NAME.bed -fo $EXONS_NAME.fa
 rm $EXONS_NAME.bed
# Step 3. concate exons
 
# fastaFromBed -s -name -fi $GENOME -bed $EXONS_NAME.bed -fo $EXONS_NAME.fa

# Step 3. Index the exon sequences
#python2.7 $INDEX_PY -d 10 -k 18 -r $EXONS_NAME.fa -o $EXONS_NAME.index

# Step 3. Construct graph and identify supporting reads
# python2.7 $GRAPH_PY -d 10 -i genes.index.gz -r genes.fa -f $END1 -b $END2 -o $GRAPH.fa

# Step 4. Construct graph and identify supporting reads
# cat $GRAPH.fa | python2.7 $RMDUP_PY - > $GRAPH.uniq.fa 