#ifndef EIGEN_QR_LAPACKE_H
#define EIGEN_QR_LAPACKE_H

namespace Eigen { 

namespace internal {

#define EIGEN_LAPACKE_QR_NOPIV(EIGTYPE, LAPACKE_TYPE, LAPACKE_PREFIX) \
template<typename MatrixQR, typename HCoeffs> \
struct householder_qr_inplace_blocked<MatrixQR, HCoeffs, EIGTYPE, true> \
{ \
  static void run(MatrixQR& mat, HCoeffs& hCoeffs, Index = 32, \
      typename MatrixQR::Scalar* = 0) \
  { \
    lapack_int m = (lapack_int) mat.rows(); \
    lapack_int n = (lapack_int) mat.cols(); \
    lapack_int lda = (lapack_int) mat.outerStride(); \
    lapack_int matrix_order = (MatrixQR::IsRowMajor) ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR; \
    LAPACKE_##LAPACKE_PREFIX##geqrf( matrix_order, m, n, (LAPACKE_TYPE*)mat.data(), lda, (LAPACKE_TYPE*)hCoeffs.data()); \
    hCoeffs.adjointInPlace(); \
  } \
};

EIGEN_LAPACKE_QR_NOPIV(double, double, d)
EIGEN_LAPACKE_QR_NOPIV(float, float, s)
EIGEN_LAPACKE_QR_NOPIV(dcomplex, lapack_complex_double, z)
EIGEN_LAPACKE_QR_NOPIV(scomplex, lapack_complex_float, c)

} 

} 

#endif 
