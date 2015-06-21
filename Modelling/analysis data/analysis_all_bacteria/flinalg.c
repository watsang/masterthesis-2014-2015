// Top-level module docstring
static char flinalg_module_doc[] =
"flinalg.c \n\
  Pointer-based (no BLAS) fast linear algebra routines. \n\
\n\
The routines in this module are UNSAFE.  They do no checks on their \n\
inputs, and are meant to be used only by code which is 100% sure that it \n\
will pass them the proper inputs from Python.";

// C code begins
#include "Python.h"
#include "numpy/arrayobject.h"

// For the intel compiler
//#ifdef __INTEL_COMPILER
//#include "mkl_cblas.h"
//#endif

// Regular blas
//#include "cblas.h"

// Our own exception
static PyObject *flinalgError;

//--------------------------------------------------------------------------
static char doc_innerproduct[] =
"inner_add(a,b,out): UNSAFE but fast inner product for double precision.\n\
\n\
Computes: \n\
  out = innerproduct(a,b) \n\
\n\
Inputs: \n\
  - a,b: numpy arrays of double-precision elements.  They must be of \n\
         compatible dimensions for the inner product to be computed, as \n\
         no checks of any kind are made.\n\
\n\
WARNING:\n\
This routine does NO CHECKS of ANY KIND!!!  It will happily segfault or \n\
give incorrect results if the inputs are not contiguous, not doubles, of \n\
the wrong dimensionality, etc.  It is strictly meant to be used for speed\n\
at the cost of safety.  For debugging, use Numpy's innerproduct() instead.";

static PyObject *array_innerproduct(PyObject *dummy, PyObject *args) {

    PyObject *op1, *op2;
    PyArrayObject *ap1, *ap2, *ret=NULL;
    double *ip1, *ip2, *op;
    npy_intp dimensions[NPY_MAXDIMS];
    double tmp;
    int nd;
    int i, j, l, i1, i2, n1, n2;

    if (!PyArg_ParseTuple(args, "OO", &op1, &op2)) return NULL;

    // This routine will be happy to segfault if the input is anything but a
    // Numpy array, since we don't check anything, we just cast the inputs
    // directly as array objects.
    ap1 = (PyArrayObject *)(op1);
    ap2 = (PyArrayObject *)(op2);

    // At least check for contiguity, which if not satisfied, can lead to
    // really bizarre and hard-to-track errors.
    if (!PyArray_ISCONTIGUOUS(ap1) || !PyArray_ISCONTIGUOUS(ap2)) {
        PyErr_SetString(flinalgError,
                        "flinalg routines only work with contiguous arrays");
        goto fail;
    }
    
    // Determine input array sizes
    l  = ap1->dimensions[ap1->nd-1];
    n1 = PyArray_SIZE(ap1)/l;
    n2 = PyArray_SIZE(ap2)/l;

    // Compute dimensionality for output
    nd = ap1->nd+ap2->nd-2;
    j = 0;
    for(i=0; i < ap1->nd-1; i++) {
        dimensions[j++] = ap1->dimensions[i];
    }
    for(i=0; i < ap2->nd-1; i++) {
        dimensions[j++] = ap2->dimensions[i];
    }
    // Allocate output
    // We simply *assume* that both inputs have the same type (and it better
    // be double).
    ret = (PyArrayObject *)PyArray_SimpleNew(nd,dimensions,
                                             ap1->descr->type_num);
    if (ret == NULL) goto fail;

    // Get pointers to input and ouptut memory zones and compute dot product
    op  = (double*)ret->data;
    ip1 = (double*)ap1->data;

    for(i1=0; i1<n1; i1++) {
        ip2 = (double*)ap2->data;
        for(i2=0; i2<n2; i2++) {
            // It's better to use a local temp than to dereference the *op
            // pointer on every pass (I measured it).
            // Compute dot product with an in-place add
            tmp = 0.0;
            for(i=0;i<l;i++) {
                tmp += ip1[i]*ip2[i];
            }
            // Accumulate result into output array value            
            *op = tmp;
            // Advance pointers for next row and output array element
            ip2 += l;
            op++;
        }
        ip1 += l;
    }

    return (PyObject *)ret;
    
  fail:
    Py_XDECREF(ret);
    return NULL;
}

//--------------------------------------------------------------------------
static char doc_inner_add[] =
"inner_add(a,b,out): UNSAFE inner product with in-place addition.\n\
\n\
Computes: \n\
  out += innerproduct(a,b) \n\
\n\
Inputs: \n\
  - a,b: numpy arrays of double-precision elements.  They must be of \n\
         compatible dimensions for the inner product to be computed, as \n\
         no checks of any kind are made.\n\
  - out: numpy array where output is accumulated.\n\
\n\
WARNING:\n\
This routine does NO CHECKS of ANY KIND!!!  It will happily segfault or \n\
give incorrect results if the inputs are not contiguous, not doubles, of \n\
the wrong dimensionality, etc.  It is strictly meant to be used for speed\n\
at the cost of safety.  For debugging, use Numpy's innerproduct() instead.";

static PyObject *array_inner_add(PyObject *dummy, PyObject *args) {

    PyObject *op1, *op2, *op3;
    PyArrayObject *ap1, *ap2, *ret;
    double *ip1, *ip2, *op;
    double tmp;
    int i, l, i1, i2, n1, n2;

    if (!PyArg_ParseTuple(args, "OOO", &op1, &op2,&op3)) return NULL;

    // This routine will segfault if the input is anything but a Numpy array,
    // since we don't check anything.
    ap1 = (PyArrayObject *)(op1);
    ap2 = (PyArrayObject *)(op2);
    ret = (PyArrayObject *)(op3);

    // At least check for contiguity, which if not satisfied, can lead to
    // really bizarre and hard-to-track errors.
    if (!PyArray_ISCONTIGUOUS(ap1) || !PyArray_ISCONTIGUOUS(ap2) ||
        !PyArray_ISCONTIGUOUS(ret) ) {
        PyErr_SetString(flinalgError,
                        "flinalg routines only work with contiguous arrays");
        return NULL;
    }
    
    // Determine input array sizes
    l  = ap1->dimensions[ap1->nd-1];
    n1 = PyArray_SIZE(ap1)/l;
    n2 = PyArray_SIZE(ap2)/l;

    /*  // debugging
    printf("\n");
    printf("n1,n2 : %d %d\n",n1,n2);
    //printf("elsize: %d \n",os);
    printf("l     : %d\n",l);
    printf("is1   : %d\n",ap1->strides[ap1->nd-1]);
    printf("is2   : %d\n",ap2->strides[ap2->nd-1]);
    */

    // Get pointers to input and ouptut memory zones and compute dot product
    op  = (double*)ret->data;
    ip1 = (double*)ap1->data;

    for(i1=0; i1<n1; i1++) {
        ip2 = (double*)ap2->data;
        for(i2=0; i2<n2; i2++) {
            // It's better to use a local temp than to dereference the *op
            // pointer on every pass (I measured it).
            // Compute dot product with an in-place add
            tmp = 0.0;
            for(i=0;i<l;i++) {
                tmp += ip1[i]*ip2[i];
            }
            // Add result into output array value            
            *op += tmp;
            // Advance pointers for next row and output array element
            ip2 += l;
            op++;
        }
        ip1 += l;
    }

    /*
    // blas version.  MUCH slower.
    for(i1=0; i1<n1; i1++) {
        ip2 = (double*)ap2->data;
        for(i2=0; i2<n2; i2++) {
            // It's better to use a local temp than to dereference the *op
            // pointer on every pass (I measured it).
            // Compute dot product with an in-place add
            // Add result into output array value            
            *op += cblas_ddot(l,ip1,1,ip2,1);

            // Advance pointers for next row and output array element
            ip2 += l;
            op++;
        }
        ip1 += l;
    }
    */

    // Return None, since we modify the output in-place.
    Py_INCREF(Py_None);
    return Py_None;
}

//--------------------------------------------------------------------------
/* registration table -- the list of methods defined in the module */
static struct PyMethodDef flinalg_module_methods[] = {
    {"inner", (PyCFunction)array_innerproduct,METH_VARARGS,doc_innerproduct},
    {"inner_add", (PyCFunction)array_inner_add,METH_VARARGS,doc_inner_add},

    // Never remove the sentinel below
    {NULL, (PyCFunction)NULL, 0, NULL}          /* sentinel */
};

//--------------------------------------------------------------------------
/* Initialization function for the module (*must* be called initArray) */
PyMODINIT_FUNC
initflinalg(void) {
    PyObject *m;

    /* Create the module and add the functions */
    m = Py_InitModule4("flinalg",
                       flinalg_module_methods,
                       flinalg_module_doc,
                       (PyObject*)NULL,
                       PYTHON_API_VERSION);

    // Add our own exception
    flinalgError = PyErr_NewException("flinalg.error", NULL, NULL);
    Py_INCREF(flinalgError);
    PyModule_AddObject(m, "error", flinalgError);
    
    /* Import the array object */
    import_array();

    /* Check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module flinalg");
}
//*************************** <END flinalg.c> *****************************
