# qsub -pe threaded 20 -N predict_1 -m abe test.sh
#gcc src/predict.c src/common.c  src/utils.c -o bin/predict -lz
#gcc src/index.c src/common.c src/fasta_uthash.c src/kmer_uthash.c -o bin/index -lz
gcc -g -O2 src/alignment.c src/utils.c src/kstring.c -o bin/alignment -lz

#cd sample_data
#../bin/index exons.fa.gz 20
#../bin/predict exons.fa.gz A431-1-ABGHI_S1_L001_R1_001.fastq.gz A431-1-ABGHI_S1_L001_R2_001.fastq.gz 1
#../bin/predict exons.fa.gz U2OS-1-ABGHI_S29_L001_R1_001.fastq.gz U2OS-1-ABGHI_S29_L001_R2_001.fastq.gz 1 5
