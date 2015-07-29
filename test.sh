# qsub -pe threaded 20 -N predict_1 -m abe test.sh
gcc src/predict.c src/kstring.c -o bin/predict -lz
#gcc -g -O2 src/alignment.c src/utils.c src/kstring.c -o bin/alignment -lz

cd sample_data
../bin/predict exons.fa U2OS-1-ABGHI_S29_L001_R1_001.fastq.gz U2OS-1-ABGHI_S29_L001_R2_001.fastq.gz 10 2
