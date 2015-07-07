#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int FastaReader(char *fasta_file)
{
	gzFile fp;
	kseq_t *seq;
	int l;
//	if (argc == 1) {
//		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
//		return 1;
//	}
	fp = gzopen(fasta_file, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		printf("name: %s\n", seq->name.s);
		if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		printf("seq: %s\n", seq->seq.s);
		if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
	}
	printf("return value: %d\n", l);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
