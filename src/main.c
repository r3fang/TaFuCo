#include <stdio.h>
#include <string.h>
#include "kstring.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.30-r15"
#endif

int main_exon_seq(int argc, char *argv[]);
int main_prefict(int argc, char *argv[]);

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: tfc (targeted gene fusion calling)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Rongxin Fang <r3fang@ucsd.edu>\n\n");
	fprintf(stderr, "Usage:   tfc <command> [options]\n\n");
	fprintf(stderr, "Command: name2fasta     extract exon sequences by gene name\n");
	fprintf(stderr, "         predict        predict gene fusion\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int i, ret;
	double t_real;
	kstring_t pg = {0,0,0};
	ksprintf(&pg, "@PG\tID:tfc\tPN:tfc\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	if (argc < 2) return usage();
	else if (strcmp(argv[1], "name2fasta") == 0) ret = main_exon_seq(argc-1, argv+1);
	else if (strcmp(argv[1], "predict") == 0) ret = main_prefict(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	if (ret == 0) {
		fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		fprintf(stderr, "[%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n");
	}
	return ret;
}