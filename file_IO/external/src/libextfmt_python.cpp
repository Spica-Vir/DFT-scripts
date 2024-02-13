// A python wrapper for libextfmt
#include <Python.h>
#include "libextfmt.h"

// Wrapper function for CUBEtoXSF_3D
static PyObject* CUBEtoXSF_3D(PyObject* self, PyObject* args) {
    const char* cubeFileName;
    const char* xsfFileName;
    const char* datatype;

    if (!PyArg_ParseTuple(args, "sss", &cubeFileName, &xsfFileName, &datatype)) {
        return nullptr;
    }

    int out = CUBEtoXSF_3D(cubeFileName, xsfFileName, datatype);
    PyObject* outpy = PyLong_FromLong(out);

    return outpy;
}

// Module method definitions
static PyMethodDef libextfmtMethods[] = {
    {"CUBEtoXSF_3D", CUBEtoXSF_3D, METH_VARARGS, "Convert 3D data format CUBE to XSF"},
    {nullptr, nullptr, 0, nullptr} // definition end indicator
};

// Module definition
static struct PyModuleDef libextfmtModule = {
    PyModuleDef_HEAD_INIT,
    "libextfmt", // python import module name
    nullptr, // description
    -1,
    libextfmtMethods
};

// Module initialization function
PyMODINIT_FUNC PyInit_libextfmt (void) {
    return PyModule_Create(&libextfmtModule);
}
