#/bin/bash

GENOME="sample_data/exons.fa"
FQ1="sample_data/reads_simulated_R1.fq"
FQ2="sample_data/reads_simulated_R2.fq"
GENE="sample_data/genes.name.txt"
TFC="bin/tfc"

cd /bioinfoSD/users/rfang/tfc
$TFC $GENE $GENOME $FQ1 $FQ2  
