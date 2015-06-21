import scipy.linalg.blas
cdef extern from "f2pyptr.h":
    void *f2py_pointer(object) except NULL
ctypedef int dgemm_t(
  char *transa, char *transb,
	int *m, int *n, int *k,
	double *alpha,
	double *a, int *lda,
	double *b, int *ldb,
	double *beta,
	double *c, int *ldc)
# Since Scipy >= 0.12.0
cdef dgemm_t *dgemm = <dgemm_t*>f2py_pointer(scipy.linalg.blas.dgemm._cpointer)
#
cdef dot(double[::1,:] a, double[::1,:] b, double[::1,:] c,
         char* transa='N', char* transb='N',
         double alpha=1., double beta=0.):
    cdef int m, n, k, lda, ldb, ldc
    cdef char* N='N'
    if(transa==N):
        m = a.shape[0]
        k = a.shape[1]
        n = b.shape[1]
    else:
        m = a.shape[1]
        k = a.shape[0]
        n = b.shape[1]
    lda = a.shape[0]
    ldb = b.shape[0]
    ldc = c.shape[0]
    dgemm(transa, transb, &m, &n, &k, &alpha, &a[0,0], &lda, &b[0,0], &ldb,
          &beta, &c[0,0], &ldc)