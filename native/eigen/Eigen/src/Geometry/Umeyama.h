#ifndef EIGEN_UMEYAMA_H
#define EIGEN_UMEYAMA_H

namespace Eigen { 

#ifndef EIGEN_PARSED_BY_DOXYGEN

namespace internal {

template<typename MatrixType, typename OtherMatrixType>
struct umeyama_transform_matrix_type
{
  enum {
    MinRowsAtCompileTime = EIGEN_SIZE_MIN_PREFER_DYNAMIC(MatrixType::RowsAtCompileTime, OtherMatrixType::RowsAtCompileTime),

    HomogeneousDimension = int(MinRowsAtCompileTime) == Dynamic ? Dynamic : int(MinRowsAtCompileTime)+1
  };

  typedef Matrix<typename traits<MatrixType>::Scalar,
    HomogeneousDimension,
    HomogeneousDimension,
    AutoAlign | (traits<MatrixType>::Flags & RowMajorBit ? RowMajor : ColMajor),
    HomogeneousDimension,
    HomogeneousDimension
  > type;
};

}

#endif

template <typename Derived, typename OtherDerived>
typename internal::umeyama_transform_matrix_type<Derived, OtherDerived>::type
umeyama(const MatrixBase<Derived>& src, const MatrixBase<OtherDerived>& dst, bool with_scaling = true)
{
  typedef typename internal::umeyama_transform_matrix_type<Derived, OtherDerived>::type TransformationMatrixType;
  typedef typename internal::traits<TransformationMatrixType>::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;

  EIGEN_STATIC_ASSERT(!NumTraits<Scalar>::IsComplex, NUMERIC_TYPE_MUST_BE_REAL)
  EIGEN_STATIC_ASSERT((internal::is_same<Scalar, typename internal::traits<OtherDerived>::Scalar>::value),
    YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)

  enum { Dimension = EIGEN_SIZE_MIN_PREFER_DYNAMIC(Derived::RowsAtCompileTime, OtherDerived::RowsAtCompileTime) };

  typedef Matrix<Scalar, Dimension, 1> VectorType;
  typedef Matrix<Scalar, Dimension, Dimension> MatrixType;
  typedef typename internal::plain_matrix_type_row_major<Derived>::type RowMajorMatrixType;

  const Index m = src.rows(); 
  const Index n = src.cols(); 

  const RealScalar one_over_n = RealScalar(1) / static_cast<RealScalar>(n);

  const VectorType src_mean = src.rowwise().sum() * one_over_n;
  const VectorType dst_mean = dst.rowwise().sum() * one_over_n;

  const RowMajorMatrixType src_demean = src.colwise() - src_mean;
  const RowMajorMatrixType dst_demean = dst.colwise() - dst_mean;

  const Scalar src_var = src_demean.rowwise().squaredNorm().sum() * one_over_n;

  const MatrixType sigma = one_over_n * dst_demean * src_demean.transpose();

  JacobiSVD<MatrixType> svd(sigma, ComputeFullU | ComputeFullV);

  TransformationMatrixType Rt = TransformationMatrixType::Identity(m+1,m+1);

  VectorType S = VectorType::Ones(m);

  if  ( svd.matrixU().determinant() * svd.matrixV().determinant() < 0 )
    S(m-1) = -1;

  Rt.block(0,0,m,m).noalias() = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();

  if (with_scaling)
  {
    
    const Scalar c = Scalar(1)/src_var * svd.singularValues().dot(S);

    Rt.col(m).head(m) = dst_mean;
    Rt.col(m).head(m).noalias() -= c*Rt.topLeftCorner(m,m)*src_mean;
    Rt.block(0,0,m,m) *= c;
  }
  else
  {
    Rt.col(m).head(m) = dst_mean;
    Rt.col(m).head(m).noalias() -= Rt.topLeftCorner(m,m)*src_mean;
  }

  return Rt;
}

} 

#endif 
