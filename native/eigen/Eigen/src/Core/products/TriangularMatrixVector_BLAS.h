#ifndef EIGEN_TRIANGULAR_MATRIX_VECTOR_BLAS_H
#define EIGEN_TRIANGULAR_MATRIX_VECTOR_BLAS_H

namespace Eigen { 

namespace internal {

template<typename Index, int Mode, typename LhsScalar, bool ConjLhs, typename RhsScalar, bool ConjRhs, int StorageOrder>
struct triangular_matrix_vector_product_trmv :
  triangular_matrix_vector_product<Index,Mode,LhsScalar,ConjLhs,RhsScalar,ConjRhs,StorageOrder,BuiltIn> {};

#define EIGEN_BLAS_TRMV_SPECIALIZE(Scalar) \
template<typename Index, int Mode, bool ConjLhs, bool ConjRhs> \
struct triangular_matrix_vector_product<Index,Mode,Scalar,ConjLhs,Scalar,ConjRhs,ColMajor,Specialized> { \
 static void run(Index _rows, Index _cols, const Scalar* _lhs, Index lhsStride, \
                                     const Scalar* _rhs, Index rhsIncr, Scalar* _res, Index resIncr, Scalar alpha) { \
      triangular_matrix_vector_product_trmv<Index,Mode,Scalar,ConjLhs,Scalar,ConjRhs,ColMajor>::run( \
        _rows, _cols, _lhs, lhsStride, _rhs, rhsIncr, _res, resIncr, alpha); \
  } \
}; \
template<typename Index, int Mode, bool ConjLhs, bool ConjRhs> \
struct triangular_matrix_vector_product<Index,Mode,Scalar,ConjLhs,Scalar,ConjRhs,RowMajor,Specialized> { \
 static void run(Index _rows, Index _cols, const Scalar* _lhs, Index lhsStride, \
                                     const Scalar* _rhs, Index rhsIncr, Scalar* _res, Index resIncr, Scalar alpha) { \
      triangular_matrix_vector_product_trmv<Index,Mode,Scalar,ConjLhs,Scalar,ConjRhs,RowMajor>::run( \
        _rows, _cols, _lhs, lhsStride, _rhs, rhsIncr, _res, resIncr, alpha); \
  } \
};

EIGEN_BLAS_TRMV_SPECIALIZE(double)
EIGEN_BLAS_TRMV_SPECIALIZE(float)
EIGEN_BLAS_TRMV_SPECIALIZE(dcomplex)
EIGEN_BLAS_TRMV_SPECIALIZE(scomplex)

#define EIGEN_BLAS_TRMV_CM(EIGTYPE, BLASTYPE, EIGPREFIX, BLASPREFIX, BLASPOSTFIX) \
template<typename Index, int Mode, bool ConjLhs, bool ConjRhs> \
struct triangular_matrix_vector_product_trmv<Index,Mode,EIGTYPE,ConjLhs,EIGTYPE,ConjRhs,ColMajor> { \
  enum { \
    IsLower = (Mode&Lower) == Lower, \
    SetDiag = (Mode&(ZeroDiag|UnitDiag)) ? 0 : 1, \
    IsUnitDiag  = (Mode&UnitDiag) ? 1 : 0, \
    IsZeroDiag  = (Mode&ZeroDiag) ? 1 : 0, \
    LowUp = IsLower ? Lower : Upper \
  }; \
 static void run(Index _rows, Index _cols, const EIGTYPE* _lhs, Index lhsStride, \
                 const EIGTYPE* _rhs, Index rhsIncr, EIGTYPE* _res, Index resIncr, EIGTYPE alpha) \
 { \
   if (ConjLhs || IsZeroDiag) { \
     triangular_matrix_vector_product<Index,Mode,EIGTYPE,ConjLhs,EIGTYPE,ConjRhs,ColMajor,BuiltIn>::run( \
       _rows, _cols, _lhs, lhsStride, _rhs, rhsIncr, _res, resIncr, alpha); \
     return; \
   }\
   Index size = (std::min)(_rows,_cols); \
   Index rows = IsLower ? _rows : size; \
   Index cols = IsLower ? size : _cols; \
\
   typedef VectorX##EIGPREFIX VectorRhs; \
   EIGTYPE *x, *y;\
\
 \
   Map<const VectorRhs, 0, InnerStride<> > rhs(_rhs,cols,InnerStride<>(rhsIncr)); \
   VectorRhs x_tmp; \
   if (ConjRhs) x_tmp = rhs.conjugate(); else x_tmp = rhs; \
   x = x_tmp.data(); \
\
\
\
   char trans, uplo, diag; \
   BlasIndex m, n, lda, incx, incy; \
   EIGTYPE const *a; \
   EIGTYPE beta(1); \
\
 \
   n = convert_index<BlasIndex>(size); \
   lda = convert_index<BlasIndex>(lhsStride); \
   incx = 1; \
   incy = convert_index<BlasIndex>(resIncr); \
\
 \
   trans = 'N'; \
   uplo = IsLower ? 'L' : 'U'; \
   diag = IsUnitDiag ? 'U' : 'N'; \
\
 \
   BLASPREFIX##trmv##BLASPOSTFIX(&uplo, &trans, &diag, &n, (const BLASTYPE*)_lhs, &lda, (BLASTYPE*)x, &incx); \
\
 \
   BLASPREFIX##axpy##BLASPOSTFIX(&n, (const BLASTYPE*)&numext::real_ref(alpha),(const BLASTYPE*)x, &incx, (BLASTYPE*)_res, &incy); \
 \
   if (size<(std::max)(rows,cols)) { \
     if (ConjRhs) x_tmp = rhs.conjugate(); else x_tmp = rhs; \
     x = x_tmp.data(); \
     if (size<rows) { \
       y = _res + size*resIncr; \
       a = _lhs + size; \
       m = convert_index<BlasIndex>(rows-size); \
       n = convert_index<BlasIndex>(size); \
     } \
     else { \
       x += size; \
       y = _res; \
       a = _lhs + size*lda; \
       m = convert_index<BlasIndex>(size); \
       n = convert_index<BlasIndex>(cols-size); \
     } \
     BLASPREFIX##gemv##BLASPOSTFIX(&trans, &m, &n, (const BLASTYPE*)&numext::real_ref(alpha), (const BLASTYPE*)a, &lda, (const BLASTYPE*)x, &incx, (const BLASTYPE*)&numext::real_ref(beta), (BLASTYPE*)y, &incy); \
   } \
  } \
};

#ifdef EIGEN_USE_MKL
EIGEN_BLAS_TRMV_CM(double,   double, d,  d,)
EIGEN_BLAS_TRMV_CM(dcomplex, MKL_Complex16, cd, z,)
EIGEN_BLAS_TRMV_CM(float,    float,  f,  s,)
EIGEN_BLAS_TRMV_CM(scomplex, MKL_Complex8,  cf, c,)
#else
EIGEN_BLAS_TRMV_CM(double,   double, d,  d, _)
EIGEN_BLAS_TRMV_CM(dcomplex, double, cd, z, _)
EIGEN_BLAS_TRMV_CM(float,    float,  f,  s, _)
EIGEN_BLAS_TRMV_CM(scomplex, float,  cf, c, _)
#endif

#define EIGEN_BLAS_TRMV_RM(EIGTYPE, BLASTYPE, EIGPREFIX, BLASPREFIX, BLASPOSTFIX) \
template<typename Index, int Mode, bool ConjLhs, bool ConjRhs> \
struct triangular_matrix_vector_product_trmv<Index,Mode,EIGTYPE,ConjLhs,EIGTYPE,ConjRhs,RowMajor> { \
  enum { \
    IsLower = (Mode&Lower) == Lower, \
    SetDiag = (Mode&(ZeroDiag|UnitDiag)) ? 0 : 1, \
    IsUnitDiag  = (Mode&UnitDiag) ? 1 : 0, \
    IsZeroDiag  = (Mode&ZeroDiag) ? 1 : 0, \
    LowUp = IsLower ? Lower : Upper \
  }; \
 static void run(Index _rows, Index _cols, const EIGTYPE* _lhs, Index lhsStride, \
                 const EIGTYPE* _rhs, Index rhsIncr, EIGTYPE* _res, Index resIncr, EIGTYPE alpha) \
 { \
   if (IsZeroDiag) { \
     triangular_matrix_vector_product<Index,Mode,EIGTYPE,ConjLhs,EIGTYPE,ConjRhs,RowMajor,BuiltIn>::run( \
       _rows, _cols, _lhs, lhsStride, _rhs, rhsIncr, _res, resIncr, alpha); \
     return; \
   }\
   Index size = (std::min)(_rows,_cols); \
   Index rows = IsLower ? _rows : size; \
   Index cols = IsLower ? size : _cols; \
\
   typedef VectorX##EIGPREFIX VectorRhs; \
   EIGTYPE *x, *y;\
\
 \
   Map<const VectorRhs, 0, InnerStride<> > rhs(_rhs,cols,InnerStride<>(rhsIncr)); \
   VectorRhs x_tmp; \
   if (ConjRhs) x_tmp = rhs.conjugate(); else x_tmp = rhs; \
   x = x_tmp.data(); \
\
\
\
   char trans, uplo, diag; \
   BlasIndex m, n, lda, incx, incy; \
   EIGTYPE const *a; \
   EIGTYPE beta(1); \
\
 \
   n = convert_index<BlasIndex>(size); \
   lda = convert_index<BlasIndex>(lhsStride); \
   incx = 1; \
   incy = convert_index<BlasIndex>(resIncr); \
\
 \
   trans = ConjLhs ? 'C' : 'T'; \
   uplo = IsLower ? 'U' : 'L'; \
   diag = IsUnitDiag ? 'U' : 'N'; \
\
 \
   BLASPREFIX##trmv##BLASPOSTFIX(&uplo, &trans, &diag, &n, (const BLASTYPE*)_lhs, &lda, (BLASTYPE*)x, &incx); \
\
 \
   BLASPREFIX##axpy##BLASPOSTFIX(&n, (const BLASTYPE*)&numext::real_ref(alpha),(const BLASTYPE*)x, &incx, (BLASTYPE*)_res, &incy); \
 \
   if (size<(std::max)(rows,cols)) { \
     if (ConjRhs) x_tmp = rhs.conjugate(); else x_tmp = rhs; \
     x = x_tmp.data(); \
     if (size<rows) { \
       y = _res + size*resIncr; \
       a = _lhs + size*lda; \
       m = convert_index<BlasIndex>(rows-size); \
       n = convert_index<BlasIndex>(size); \
     } \
     else { \
       x += size; \
       y = _res; \
       a = _lhs + size; \
       m = convert_index<BlasIndex>(size); \
       n = convert_index<BlasIndex>(cols-size); \
     } \
     BLASPREFIX##gemv##BLASPOSTFIX(&trans, &n, &m, (const BLASTYPE*)&numext::real_ref(alpha), (const BLASTYPE*)a, &lda, (const BLASTYPE*)x, &incx, (const BLASTYPE*)&numext::real_ref(beta), (BLASTYPE*)y, &incy); \
   } \
  } \
};

#ifdef EIGEN_USE_MKL
EIGEN_BLAS_TRMV_RM(double,   double, d,  d,)
EIGEN_BLAS_TRMV_RM(dcomplex, MKL_Complex16, cd, z,)
EIGEN_BLAS_TRMV_RM(float,    float,  f,  s,)
EIGEN_BLAS_TRMV_RM(scomplex, MKL_Complex8,  cf, c,)
#else
EIGEN_BLAS_TRMV_RM(double,   double, d,  d,_)
EIGEN_BLAS_TRMV_RM(dcomplex, double, cd, z,_)
EIGEN_BLAS_TRMV_RM(float,    float,  f,  s,_)
EIGEN_BLAS_TRMV_RM(scomplex, float,  cf, c,_)
#endif

} 

} 

#endif 
