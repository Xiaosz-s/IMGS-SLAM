#ifndef EIGEN_PARTIALLU_LAPACK_H
#define EIGEN_PARTIALLU_LAPACK_H

namespace Eigen { 

namespace internal {

#define EIGEN_LAPACKE_LU_PARTPIV(EIGTYPE, LAPACKE_TYPE, LAPACKE_PREFIX) \
template<int StorageOrder> \
struct partial_lu_impl<EIGTYPE, StorageOrder, lapack_int> \
{ \
   \
  static lapack_int blocked_lu(Index rows, Index cols, EIGTYPE* lu_data, Index luStride, lapack_int* row_transpositions, lapack_int& nb_transpositions, lapack_int maxBlockSize=256) \
  { \
    EIGEN_UNUSED_VARIABLE(maxBlockSize);\
    lapack_int matrix_order, first_zero_pivot; \
    lapack_int m, n, lda, *ipiv, info; \
    EIGTYPE* a; \
 \
    matrix_order = StorageOrder==RowMajor ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR; \
    lda = convert_index<lapack_int>(luStride); \
    a = lu_data; \
    ipiv = row_transpositions; \
    m = convert_index<lapack_int>(rows); \
    n = convert_index<lapack_int>(cols); \
    nb_transpositions = 0; \
\
    info = LAPACKE_##LAPACKE_PREFIX##getrf( matrix_order, m, n, (LAPACKE_TYPE*)a, lda, ipiv ); \
\
    for(int i=0;i<m;i++) { ipiv[i]--; if (ipiv[i]!=i) nb_transpositions++; } \
\
    eigen_assert(info >= 0); \
 \
\
    first_zero_pivot = info; \
    return first_zero_pivot; \
  } \
};

EIGEN_LAPACKE_LU_PARTPIV(double, double, d)
EIGEN_LAPACKE_LU_PARTPIV(float, float, s)
EIGEN_LAPACKE_LU_PARTPIV(dcomplex, lapack_complex_double, z)
EIGEN_LAPACKE_LU_PARTPIV(scomplex, lapack_complex_float,  c)

} 

} 

#endif 
