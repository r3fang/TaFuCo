all:
		$(CC) -g -O2 src/main.c src/predict.c src/kstring.c -o tfc -lz  -lm