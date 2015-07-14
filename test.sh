gcc src/predict.c src/common.c -o bin/predict -lz
# gcc src/index.c src/common.c -o bin/index -lz

cd sample_data
#../bin/index exons.fa.gz 20
../bin/predict exons.fa.gz U2OS-1-ABGHI_S29_L001_R1_001.fastq.gz U2OS-1-ABGHI_S29_L001_R2_001.fastq.gz
