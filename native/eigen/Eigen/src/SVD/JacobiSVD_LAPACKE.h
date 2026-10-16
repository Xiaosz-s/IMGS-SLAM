#ifndef EIGEN_JACOBISVD_LAPACKE_H
#define EIGEN_JACOBISVD_LAPACKE_H

namespace Eigen { 

#define EIGEN_LAPACKE_SVD(EIGTYPE, LAPACKE_TYPE, LAPACKE_RTYPE, LAPACKE_PREFIX, EIGCOLROW, LAPACKE_COLROW) \
template<> inline \
JacobiSVD<Matrix<EIGTYPE, Dynamic, Dynamic, EIGCOLROW, Dynamic, Dynamic>, ColPivHouseholderQRPreconditioner>& \
JacobiSVD<Matrix<EIGTYPE, Dynamic, Dynamic, EIGCOLROW, Dynamic, Dynamic>, ColPivHouseholderQRPreconditioner>::compute(const Matrix<EIGTYPE, Dynamic, Dynamic, EIGCOLROW, Dynamic, Dynamic>& matrix, unsigned int computationOptions) \
{ \
  typedef Matrix<EIGTYPE, Dynamic, Dynamic, EIGCOLROW, Dynamic, Dynamic> MatrixType; \
   \
   \
  allocate(matrix.rows(), matrix.cols(), computationOptions); \
\
   \
  m_nonzeroSingularValues = m_diagSize; \
\
  lapack_int lda = internal::convert_index<lapack_int>(matrix.outerStride()), ldu, ldvt; \
  lapack_int matrix_order = LAPACKE_COLROW; \
  char jobu, jobvt; \
  LAPACKE_TYPE *u, *vt, dummy; \
  jobu  = (m_computeFullU) ? 'A' : (m_computeThinU) ? 'S' : 'N'; \
  jobvt = (m_computeFullV) ? 'A' : (m_computeThinV) ? 'S' : 'N'; \
  if (computeU()) { \
    ldu  = internal::convert_index<lapack_int>(m_matrixU.outerStride()); \
    u    = (LAPACKE_TYPE*)m_matrixU.data(); \
  } else { ldu=1; u=&dummy; }\
  MatrixType localV; \
  lapack_int vt_rows = (m_computeFullV) ? internal::convert_index<lapack_int>(m_cols) : (m_computeThinV) ? internal::convert_index<lapack_int>(m_diagSize) : 1; \
  if (computeV()) { \
    localV.resize(vt_rows, m_cols); \
    ldvt  = internal::convert_index<lapack_int>(localV.outerStride()); \
    vt   = (LAPACKE_TYPE*)localV.data(); \
  } else { ldvt=1; vt=&dummy; }\
  Matrix<LAPACKE_RTYPE, Dynamic, Dynamic> superb; superb.resize(m_diagSize, 1); \
  MatrixType m_temp; m_temp = matrix; \
  LAPACKE_##LAPACKE_PREFIX##gesvd( matrix_order, jobu, jobvt, internal::convert_index<lapack_int>(m_rows), internal::convert_index<lapack_int>(m_cols), (LAPACKE_TYPE*)m_temp.data(), lda, (LAPACKE_RTYPE*)m_singularValues.data(), u, ldu, vt, ldvt, superb.data()); \
  if (computeV()) m_matrixV = localV.adjoint(); \
  \
  m_isInitialized = true; \
  return *this; \
}

EIGEN_LAPACKE_SVD(double,   double,                double, d, ColMajor, LAPACK_COL_MAJOR)
EIGEN_LAPACKE_SVD(float,    float,                 float , s, ColMajor, LAPACK_COL_MAJOR)
EIGEN_LAPACKE_SVD(dcomplex, lapack_complex_double, double, z, ColMajor, LAPACK_COL_MAJOR)
EIGEN_LAPACKE_SVD(scomplex, lapack_complex_float,  float , c, ColMajor, LAPACK_COL_MAJOR)

EIGEN_LAPACKE_SVD(double,   double,                double, d, RowMajor, LAPACK_ROW_MAJOR)
EIGEN_LAPACKE_SVD(float,    float,                 float , s, RowMajor, LAPACK_ROW_MAJOR)
EIGEN_LAPACKE_SVD(dcomplex, lapack_complex_double, double, z, RowMajor, LAPACK_ROW_MAJOR)
EIGEN_LAPACKE_SVD(scomplex, lapack_complex_float,  float , c, RowMajor, LAPACK_ROW_MAJOR)

} 

#endif 
