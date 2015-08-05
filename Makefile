all:
		$(CC) -g -O2 src/main.c src/predict.c src/kstring.c -o bin/tfc      -lz  -lm
clean:
		rm -f bin/*.dSYM