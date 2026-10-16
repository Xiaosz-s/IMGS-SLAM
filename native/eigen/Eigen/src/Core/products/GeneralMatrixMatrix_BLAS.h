#ifndef EIGEN_GENERAL_MATRIX_MATRIX_BLAS_H
#define EIGEN_GENERAL_MATRIX_MATRIX_BLAS_H

namespace Eigen { 

namespace internal {

#define GEMM_SPECIALIZATION(EIGTYPE, EIGPREFIX, BLASTYPE, BLASFUNC) \
template< \
  typename Index, \
  int LhsStorageOrder, bool ConjugateLhs, \
  int RhsStorageOrder, bool ConjugateRhs> \
struct general_matrix_matrix_product<Index,EIGTYPE,LhsStorageOrder,ConjugateLhs,EIGTYPE,RhsStorageOrder,ConjugateRhs,ColMajor,1> \
{ \
typedef gebp_traits<EIGTYPE,EIGTYPE> Traits; \
\
static void run(Index rows, Index cols, Index depth, \
  const EIGTYPE* _lhs, Index lhsStride, \
  const EIGTYPE* _rhs, Index rhsStride, \
  EIGTYPE* res, Index resIncr, Index resStride, \
  EIGTYPE alpha, \
  level3_blocking<EIGTYPE, EIGTYPE>& , \
  GemmParallelInfo<Index>* ) \
{ \
  using std::conj; \
\
  EIGEN_ONLY_USED_FOR_DEBUG(resIncr); \
  eigen_assert(resIncr == 1); \
  char transa, transb; \
  BlasIndex m, n, k, lda, ldb, ldc; \
  const EIGTYPE *a, *b; \
  EIGTYPE beta(1); \
  MatrixX##EIGPREFIX a_tmp, b_tmp; \
\
 \
  transa = (LhsStorageOrder==RowMajor) ? ((ConjugateLhs) ? 'C' : 'T') : 'N'; \
  transb = (RhsStorageOrder==RowMajor) ? ((ConjugateRhs) ? 'C' : 'T') : 'N'; \
\
 \
  m = convert_index<BlasIndex>(rows);  \
  n = convert_index<BlasIndex>(cols);  \
  k = convert_index<BlasIndex>(depth); \
\
 \
  lda = convert_index<BlasIndex>(lhsStride); \
  ldb = convert_index<BlasIndex>(rhsStride); \
  ldc = convert_index<BlasIndex>(resStride); \
\
 \
  if ((LhsStorageOrder==ColMajor) && (ConjugateLhs)) { \
    Map<const MatrixX##EIGPREFIX, 0, OuterStride<> > lhs(_lhs,m,k,OuterStride<>(lhsStride)); \
    a_tmp = lhs.conjugate(); \
    a = a_tmp.data(); \
    lda = convert_index<BlasIndex>(a_tmp.outerStride()); \
  } else a = _lhs; \
\
  if ((RhsStorageOrder==ColMajor) && (ConjugateRhs)) { \
    Map<const MatrixX##EIGPREFIX, 0, OuterStride<> > rhs(_rhs,k,n,OuterStride<>(rhsStride)); \
    b_tmp = rhs.conjugate(); \
    b = b_tmp.data(); \
    ldb = convert_index<BlasIndex>(b_tmp.outerStride()); \
  } else b = _rhs; \
\
  BLASFUNC(&transa, &transb, &m, &n, &k, (const BLASTYPE*)&numext::real_ref(alpha), (const BLASTYPE*)a, &lda, (const BLASTYPE*)b, &ldb, (const BLASTYPE*)&numext::real_ref(beta), (BLASTYPE*)res, &ldc); \
}};

#ifdef EIGEN_USE_MKL
GEMM_SPECIALIZATION(double,   d,  double, dgemm)
GEMM_SPECIALIZATION(float,    f,  float,  sgemm)
GEMM_SPECIALIZATION(dcomplex, cd, MKL_Complex16, zgemm)
GEMM_SPECIALIZATION(scomplex, cf, MKL_Complex8,  cgemm)
#else
GEMM_SPECIALIZATION(double,   d,  double, dgemm_)
GEMM_SPECIALIZATION(float,    f,  float,  sgemm_)
GEMM_SPECIALIZATION(dcomplex, cd, double, zgemm_)
GEMM_SPECIALIZATION(scomplex, cf, float,  cgemm_)
#endif

} 

} 

#endif 
