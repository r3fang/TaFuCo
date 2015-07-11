#include <Python.h>
#include <zlib.h>
#include <stdio.h>
#include <assert.h>
#include "kseq.h"
#include "uthash.h"
#include "common.h"

KSEQ_INIT(gzFile, gzread)

static PyObject *foo_try(PyObject *self, PyObject *args)
{
	char str0[] = "CTCF.exon3_176";
	char *str =  malloc((strlen(str0)+1) * sizeof(char));
	strncpy(str, str0, strlen(str0));
	/* robust strsplit */	
	int i;
	//char* exon = pos_parser(str, &i);
	//printf("%d\n", i);
	char* exon = pos_parser(str, &i);
	if(exon==NULL){
		return Py_BuildValue("");
	}
	printf("exon=%s\tpos=%d\n", exon, i);
	int num = 10;
	return Py_BuildValue("");
}

static PyObject *foo_predict(PyObject *self, PyObject *args)
{
	/* Index the reference genome. */
	char *fasta_file;
	char *fastq_file;
	int k;
	
	/* Parsing input paramters */
	if (!PyArg_ParseTuple(args, "ssi", &fasta_file, &fastq_file, &k)){
			return NULL;
	}
	if (fasta_file == NULL){
		return NULL;
	}
	predict_main(fasta_file, fastq_file, k);
	return Py_BuildValue("");
}

static PyObject *foo_index(PyObject *self, PyObject *args)
{
	/* Index the reference genome. */
	char *fasta_file;
	int k;
	/* Parsing input paramters */
	if (!PyArg_ParseTuple(args, "si", &fasta_file, &k)){
		return NULL;
	}
	if (fasta_file == NULL){
		return NULL;
	}
	index_main(fasta_file, k);
	return Py_BuildValue("");
}

static PyObject *foo_FastaReader(PyObject *self, PyObject *args){
	char *fasta_file;
	if (!PyArg_ParseTuple(args, "s", &fasta_file)){
		return NULL;
	}
	/* Initlize variable*/
	gzFile fp; kseq_t *seq; int l;
	fp = gzopen(fasta_file, "r");
	if(!fp){
		return NULL;
	}
	seq = kseq_init(fp);
	/* Create a python list to contain the (name, seq)*/
	PyObject *pylist = PyList_New(0);
	while ((l = kseq_read(seq)) >= 0) {
		PyObject *tmp = Py_BuildValue("(s,s)", seq->name.s, seq->seq.s);
		PyList_Append(pylist, tmp);
		Py_DECREF(tmp);
	}
	return pylist;
}

static PyObject *foo_ReverseComplement(PyObject *self, PyObject *args){
	/* Reverse complement of given DNA sequence*/
	char *s;
	if (!PyArg_ParseTuple(args, "s", &s)){
		return NULL;
	}
	int n, c, d;
	n=strlen(s);
	char* r;
	r = (char*)malloc(n * sizeof(char));
	for (c = n - 1, d = 0; c >= 0; c--, d++){
		switch(toupper(s[c])){
			case 'A':
				r[d] = 'T';
				break;
			case 'T':
				r[d] = 'A';
				break;
			case 'C':
				r[d] = 'G';
				break;
			case 'G':
				r[d] = 'C';
				break;
			default:
				r[d] = s[d];
				break;
		}
	}
	r[n] = '\0';
	return Py_BuildValue("s", r);
}

static PyObject *foo_kmer_match(PyObject *self, PyObject *args){
	char *ref;
	char *quary;
	if (!PyArg_ParseTuple(args, "ss", &quary, &ref)){
		return NULL;
	}
	if(strlen(quary) > strlen(ref))
		return NULL;
	int ref_len = strlen(ref);
	int quary_len = strlen(quary);
	int i;
	char* buff;
	buff = (char*)malloc(quary_len * sizeof(char));
	for(i=0; i<ref_len-quary_len; i++){
		memcpy(buff, &ref[i], quary_len);
		buff[quary_len] = '\0';
		if(strcmp(buff, quary) == 0){
			free(buff);
			return Py_BuildValue("i", i);			
		}
	}
	return Py_BuildValue("i", -1);
}

/* docstring for functions */
static char FastaReader_docs[] = 
	"Fasta parser.";
static char ReverseComplement_docs[] = 
	"Reverse complement of given DNA sequence.";
static char kmer_match_docs[] = 
	"If given kmer occurs in ref seq.";
static char index_docs[] = 
	"Index reference DNA sequence.";
static char predict_docs[] = 
	"Predict gene fusion.";
static char try_docs[] = 
	"try.";
static char foo_docs[] = 
	"A collections of non-sense functions.";

/* Method Mapping Fucntion of the module */
static PyMethodDef foo_funcs[] = {
	{"FastaReader", (PyCFunction)foo_FastaReader, METH_VARARGS, FastaReader_docs},
	{"ReverseComplement", (PyCFunction)foo_ReverseComplement, METH_VARARGS, ReverseComplement_docs},	
	{"kmer_match", (PyCFunction)foo_kmer_match, METH_VARARGS, kmer_match_docs},	
	{"index", (PyCFunction)foo_index, METH_VARARGS, index_docs},	
	{"predict", (PyCFunction)foo_predict, METH_VARARGS, predict_docs},	
	{"TRY", (PyCFunction)foo_try, METH_VARARGS, try_docs},	
	{NULL}
};

/* Initialization Function */
PyMODINIT_FUNC initfoo() {
	Py_InitModule3("foo", foo_funcs,
					foo_docs);
}
