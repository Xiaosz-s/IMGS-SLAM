#ifndef EIGEN_LLT_LAPACKE_H
#define EIGEN_LLT_LAPACKE_H

namespace Eigen { 

namespace internal {

template<typename Scalar> struct lapacke_llt;

#define EIGEN_LAPACKE_LLT(EIGTYPE, BLASTYPE, LAPACKE_PREFIX) \
template<> struct lapacke_llt<EIGTYPE> \
{ \
  template<typename MatrixType> \
  static inline Index potrf(MatrixType& m, char uplo) \
  { \
    lapack_int matrix_order; \
    lapack_int size, lda, info, StorageOrder; \
    EIGTYPE* a; \
    eigen_assert(m.rows()==m.cols()); \
     \
    size = convert_index<lapack_int>(m.rows()); \
    StorageOrder = MatrixType::Flags&RowMajorBit?RowMajor:ColMajor; \
    matrix_order = StorageOrder==RowMajor ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR; \
    a = &(m.coeffRef(0,0)); \
    lda = convert_index<lapack_int>(m.outerStride()); \
\
    info = LAPACKE_##LAPACKE_PREFIX##potrf( matrix_order, uplo, size, (BLASTYPE*)a, lda ); \
    info = (info==0) ? -1 : info>0 ? info-1 : size; \
    return info; \
  } \
}; \
template<> struct llt_inplace<EIGTYPE, Lower> \
{ \
  template<typename MatrixType> \
  static Index blocked(MatrixType& m) \
  { \
    return lapacke_llt<EIGTYPE>::potrf(m, 'L'); \
  } \
  template<typename MatrixType, typename VectorType> \
  static Index rankUpdate(MatrixType& mat, const VectorType& vec, const typename MatrixType::RealScalar& sigma) \
  { return Eigen::internal::llt_rank_update_lower(mat, vec, sigma); } \
}; \
template<> struct llt_inplace<EIGTYPE, Upper> \
{ \
  template<typename MatrixType> \
  static Index blocked(MatrixType& m) \
  { \
    return lapacke_llt<EIGTYPE>::potrf(m, 'U'); \
  } \
  template<typename MatrixType, typename VectorType> \
  static Index rankUpdate(MatrixType& mat, const VectorType& vec, const typename MatrixType::RealScalar& sigma) \
  { \
    Transpose<MatrixType> matt(mat); \
    return llt_inplace<EIGTYPE, Lower>::rankUpdate(matt, vec.conjugate(), sigma); \
  } \
};

EIGEN_LAPACKE_LLT(double, double, d)
EIGEN_LAPACKE_LLT(float, float, s)
EIGEN_LAPACKE_LLT(dcomplex, lapack_complex_double, z)
EIGEN_LAPACKE_LLT(scomplex, lapack_complex_float, c)

} 

} 

#endif 
