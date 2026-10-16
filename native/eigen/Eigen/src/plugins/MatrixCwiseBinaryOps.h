template<typename OtherDerived>
EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const EIGEN_CWISE_BINARY_RETURN_TYPE(Derived,OtherDerived,product)
cwiseProduct(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived> &other) const
{
  return EIGEN_CWISE_BINARY_RETURN_TYPE(Derived,OtherDerived,product)(derived(), other.derived());
}

template<typename OtherDerived>
EIGEN_DEVICE_FUNC
inline const CwiseBinaryOp<numext::equal_to<Scalar>, const Derived, const OtherDerived>
cwiseEqual(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived> &other) const
{
  return CwiseBinaryOp<numext::equal_to<Scalar>, const Derived, const OtherDerived>(derived(), other.derived());
}

template<typename OtherDerived>
EIGEN_DEVICE_FUNC
inline const CwiseBinaryOp<numext::not_equal_to<Scalar>, const Derived, const OtherDerived>
cwiseNotEqual(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived> &other) const
{
  return CwiseBinaryOp<numext::not_equal_to<Scalar>, const Derived, const OtherDerived>(derived(), other.derived());
}

template<typename OtherDerived>
EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_min_op<Scalar,Scalar>, const Derived, const OtherDerived>
cwiseMin(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived> &other) const
{
  return CwiseBinaryOp<internal::scalar_min_op<Scalar,Scalar>, const Derived, const OtherDerived>(derived(), other.derived());
}

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_min_op<Scalar,Scalar>, const Derived, const ConstantReturnType>
cwiseMin(const Scalar &other) const
{
  return cwiseMin(Derived::Constant(rows(), cols(), other));
}

template<typename OtherDerived>
EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_max_op<Scalar,Scalar>, const Derived, const OtherDerived>
cwiseMax(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived> &other) const
{
  return CwiseBinaryOp<internal::scalar_max_op<Scalar,Scalar>, const Derived, const OtherDerived>(derived(), other.derived());
}

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_max_op<Scalar,Scalar>, const Derived, const ConstantReturnType>
cwiseMax(const Scalar &other) const
{
  return cwiseMax(Derived::Constant(rows(), cols(), other));
}

template<typename OtherDerived>
EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_quotient_op<Scalar>, const Derived, const OtherDerived>
cwiseQuotient(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived> &other) const
{
  return CwiseBinaryOp<internal::scalar_quotient_op<Scalar>, const Derived, const OtherDerived>(derived(), other.derived());
}

typedef CwiseBinaryOp<internal::scalar_cmp_op<Scalar,Scalar,internal::cmp_EQ>, const Derived, const ConstantReturnType> CwiseScalarEqualReturnType;

EIGEN_DEVICE_FUNC
inline const CwiseScalarEqualReturnType
cwiseEqual(const Scalar& s) const
{
  return CwiseScalarEqualReturnType(derived(), Derived::Constant(rows(), cols(), s), internal::scalar_cmp_op<Scalar,Scalar,internal::cmp_EQ>());
}
