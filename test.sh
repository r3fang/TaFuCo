# qsub -pe threaded 20 -N predict_1 -m abe test.sh
gcc src/predict.c src/common.c  src/utils.c -o bin/predict -lz

#gcc src/index.c src/common.c src/fasta_uthash.c src/kmer_uthash.c -o bin/index -lz

#cd sample_data
#../bin/index exons.fa.gz 20
#../bin/predict exons.fa.gz A431-1-ABGHI_S1_L001_R1_001.fastq.gz A431-1-ABGHI_S1_L001_R2_001.fastq.gz 0 0
cat results/A431.bag.fa | python bin/rmDup.py > results/A431.bag.uniq.fa

