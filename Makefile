all:
		#$(CC) -g -O2 src/main.c     src/kstring.c -o bin/tfc      -lz  -lm
		#$(CC) -g -O2 src/junction.c src/kstring.c -o bin/junction -lz  -lm
		$(CC) -g -O2 src/exon_seq_extract.c src/kstring.c -o bin/exon_seq_extract -lz  -lm
		
clean:
		rm -f bin/*.dSYM