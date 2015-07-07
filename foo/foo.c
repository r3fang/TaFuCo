#include <Python.h>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

/* C function */
static PyObject *foo_add(PyObject *self, PyObject *args){
	int a, b;
	/* Parsing input paramters */
	if (!PyArg_ParseTuple(args, "ii", &a, &b)){
		return NULL;
	}
	return Py_BuildValue("i", a + b);
}

static PyObject *foo_strlen(PyObject *self, PyObject *args){
	char *s;
	/* Parsing input paramters */
	if (!PyArg_ParseTuple(args, "s", &s)){
		return NULL;
	}
	return Py_BuildValue("i", strlen(s));
}

static PyObject *foo_FastaReader(PyObject *self, PyObject *args){
	char *fasta_file;
	if (!PyArg_ParseTuple(args, "s", &fasta_file)){
		return NULL;
	}
	
	gzFile fp;
	kseq_t *seq;
	int l;
	PyObject *list = PyList_New(100);
	fp = gzopen(fasta_file, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		return Py_BuildValue("{s:s}", seq->name.s, seq->seq.s);
		printf("name: %s\n", seq->name.s);
		if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		printf("seq: %s\n", seq->seq.s);
		if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
	}
}
	
/* docstring for functions */
static char add_docs[] = 
	"Sum of two integers.";
static char strlen_docs[] = 
	"Length of a string.";
static char FastaReader_docs[] = 
	"Fasta parser.";

static char foo_docs[] = 
	"A collections of non-sense functions.";

/* Method Mapping Fucntion of the module */
static PyMethodDef foo_funcs[] = {
	{"add", (PyCFunction)foo_add, METH_VARARGS, add_docs},
	{"strlen", (PyCFunction)foo_strlen, METH_VARARGS, strlen_docs},
	{"FastaReader", (PyCFunction)foo_FastaReader, METH_VARARGS, FastaReader_docs},
	{NULL}
};

/* Initialization Function */
PyMODINIT_FUNC initfoo() {
	Py_InitModule3("foo", foo_funcs,
					foo_docs);
}
