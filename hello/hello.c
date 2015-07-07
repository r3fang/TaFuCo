#include <Python.h>

/* The interface*/
static PyObject* helloworld(PyObject* self)
{
	return Py_BuildValue("s", "Hello, Python extension!!");
}

/* docstrings of the function*/
static char helloworld_docs[] =
	"helloworld(): Any message you want to put here!!\n";

/* Method Mapping Table */
static PyMethodDef helloworld_funcs[] = {
	{"helloworld", (PyCFunction)helloworld,
	METH_NOARGS, helloworld_docs},
	{NULL}
};

/* Initilization */
PyMODINIT_FUNC inithelloworld()
{
	Py_InitModule3("helloworld", helloworld_funcs,
				   "Extension module example!");
}
