#include <Python.h>

/* C function */
static PyObject *foo_add(PyObject *self, PyObject *args){
	int a, b;
	/* Parsing input paramters */
	if (!PyArg_ParseTuple(args, "ii", &a, &b)){
		return NULL;
	}
	return Py_BuildValue("i", a + b);
}

/* docstring for functions */
static char add_docs[] = 
	"Sum of two integers.";
static char foo_docs[] = 
	"A collections of non-sense functions.";

/* Method Mapping Fucntion of the module */
static PyMethodDef foo_funcs[] = {
	{"add", (PyCFunction)foo_add, 
	METH_VARARGS, add_docs},
	{NULL}
};

/* Initialization Function */
PyMODINIT_FUNC initfoo() {
	Py_InitModule3("foo", foo_funcs,
					foo_docs);
}
