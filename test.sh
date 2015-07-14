gcc src/predict.c src/common.c src/index.c -o bin/predict -lz
cd sample_data
../bin/predict exons.fa.gz U2OS-1-ABGHI_S29_L001_R1_001.fastq.gz U2OS-1-ABGHI_S29_L001_R2_001.fastq.gz