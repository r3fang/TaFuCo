all:
		#$(CC) -g -O2 src/main.c src/kstring.c -o bin/tfc -lz
		$(CC) -g -O2 src/alignment.c src/kstring.c -o bin/alignment -lz
clean:
		rm -f bin/*.dSYM