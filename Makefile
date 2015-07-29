all:
		$(CC) -g -O2 src/main.c src/kstring.c -o bin/tfc -lz
clean:
		rm -f bin/*.dSYM