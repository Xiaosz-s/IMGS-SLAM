#ifndef EIGEN_COLPIVOTINGHOUSEHOLDERQR_LAPACKE_H
#define EIGEN_COLPIVOTINGHOUSEHOLDERQR_LAPACKE_H

namespace Eigen { 

#define EIGEN_LAPACKE_QR_COLPIV(EIGTYPE, LAPACKE_TYPE, LAPACKE_PREFIX, EIGCOLROW, LAPACKE_COLROW) \
template<> template<typename InputType> inline \
ColPivHouseholderQR<Matrix<EIGTYPE, Dynamic, Dynamic, EIGCOLROW, Dynamic, Dynamic> >& \
ColPivHouseholderQR<Matrix<EIGTYPE, Dynamic, Dynamic, EIGCOLROW, Dynamic, Dynamic> >::compute( \
              const EigenBase<InputType>& matrix) \
\
{ \
  using std::abs; \
  typedef Matrix<EIGTYPE, Dynamic, Dynamic, EIGCOLROW, Dynamic, Dynamic> MatrixType; \
  typedef MatrixType::RealScalar RealScalar; \
  Index rows = matrix.rows();\
  Index cols = matrix.cols();\
\
  m_qr = matrix;\
  Index size = m_qr.diagonalSize();\
  m_hCoeffs.resize(size);\
\
  m_colsTranspositions.resize(cols);\
   \
\
  m_nonzero_pivots = 0; \
  m_maxpivot = RealScalar(0);\
  m_colsPermutation.resize(cols); \
  m_colsPermutation.indices().setZero(); \
\
  lapack_int lda = internal::convert_index<lapack_int,Index>(m_qr.outerStride()); \
  lapack_int matrix_order = LAPACKE_COLROW; \
  LAPACKE_##LAPACKE_PREFIX##geqp3( matrix_order, internal::convert_index<lapack_int,Index>(rows), internal::convert_index<lapack_int,Index>(cols), \
                              (LAPACKE_TYPE*)m_qr.data(), lda, (lapack_int*)m_colsPermutation.indices().data(), (LAPACKE_TYPE*)m_hCoeffs.data()); \
  m_isInitialized = true; \
  m_maxpivot=m_qr.diagonal().cwiseAbs().maxCoeff(); \
  m_hCoeffs.adjointInPlace(); \
  RealScalar premultiplied_threshold = abs(m_maxpivot) * threshold(); \
  lapack_int *perm = m_colsPermutation.indices().data(); \
  for(Index i=0;i<size;i++) { \
    m_nonzero_pivots += (abs(m_qr.coeff(i,i)) > premultiplied_threshold);\
  } \
  for(Index i=0;i<cols;i++) perm[i]--;\
\
   \
\
  return *this; \
}

EIGEN_LAPACKE_QR_COLPIV(double,   double,        d, ColMajor, LAPACK_COL_MAJOR)
EIGEN_LAPACKE_QR_COLPIV(float,    float,         s, ColMajor, LAPACK_COL_MAJOR)
EIGEN_LAPACKE_QR_COLPIV(dcomplex, lapack_complex_double, z, ColMajor, LAPACK_COL_MAJOR)
EIGEN_LAPACKE_QR_COLPIV(scomplex, lapack_complex_float,  c, ColMajor, LAPACK_COL_MAJOR)

EIGEN_LAPACKE_QR_COLPIV(double,   double,        d, RowMajor, LAPACK_ROW_MAJOR)
EIGEN_LAPACKE_QR_COLPIV(float,    float,         s, RowMajor, LAPACK_ROW_MAJOR)
EIGEN_LAPACKE_QR_COLPIV(dcomplex, lapack_complex_double, z, RowMajor, LAPACK_ROW_MAJOR)
EIGEN_LAPACKE_QR_COLPIV(scomplex, lapack_complex_float,  c, RowMajor, LAPACK_ROW_MAJOR)

} 

#endif 
