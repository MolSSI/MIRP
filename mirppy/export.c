#include <Python.h>

/* In export_boys.c */
PyObject * export_mirp_boys_mp(PyObject *self, PyObject *args);
PyObject * export_mirp_boys_interval(PyObject *self, PyObject *args);
PyObject * export_mirp_boys_double(PyObject *self, PyObject *args);

/* In export_eri.c */
PyObject * export_mirp_single_eri_double(PyObject *self, PyObject *args);
PyObject * export_mirp_prim_eri_double(PyObject *self, PyObject *args);
PyObject * export_mirp_eri_double(PyObject *self, PyObject *args);


static PyMethodDef mirp_methods[] = {
    {"boys_mp", export_mirp_boys_mp, METH_VARARGS,
     "Test function"},
    {"boys_interval", export_mirp_boys_interval, METH_VARARGS,
     "Test function"},
    {"boys_double", export_mirp_boys_double, METH_VARARGS,
     "Test function"},
    {"single_eri_double", export_mirp_single_eri_double, METH_VARARGS,
     "Test function"},
    {"prim_eri_double", export_mirp_prim_eri_double, METH_VARARGS,
     "Test function"},
    {"eri_double", export_mirp_eri_double, METH_VARARGS,
     "Test function"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


static struct PyModuleDef mirp_module = {
   PyModuleDef_HEAD_INIT,
   "mirp",   /* name of module */
   "hello",    /* module documentation, may be NULL */
   -1,         /* size of per-interpreter state of the module,
                  or -1 if the module keeps state in global variables. */
   mirp_methods,
   NULL,
   NULL,
   NULL,
   NULL
};



PyMODINIT_FUNC
PyInit_mirppy(void)
{
    return PyModule_Create(&mirp_module);
}
