typedef CwiseUnaryOp<internal::scalar_abs_op<Scalar>, const Derived> CwiseAbsReturnType;
typedef CwiseUnaryOp<internal::scalar_abs2_op<Scalar>, const Derived> CwiseAbs2ReturnType;
typedef CwiseUnaryOp<internal::scalar_arg_op<Scalar>, const Derived> CwiseArgReturnType;
typedef CwiseUnaryOp<internal::scalar_sqrt_op<Scalar>, const Derived> CwiseSqrtReturnType;
typedef CwiseUnaryOp<internal::scalar_sign_op<Scalar>, const Derived> CwiseSignReturnType;
typedef CwiseUnaryOp<internal::scalar_inverse_op<Scalar>, const Derived> CwiseInverseReturnType;

EIGEN_DOC_UNARY_ADDONS(cwiseAbs,absolute value)

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const CwiseAbsReturnType
cwiseAbs() const { return CwiseAbsReturnType(derived()); }

EIGEN_DOC_UNARY_ADDONS(cwiseAbs2,squared absolute value)

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const CwiseAbs2ReturnType
cwiseAbs2() const { return CwiseAbs2ReturnType(derived()); }

EIGEN_DOC_UNARY_ADDONS(cwiseSqrt,square-root)

EIGEN_DEVICE_FUNC
inline const CwiseSqrtReturnType
cwiseSqrt() const { return CwiseSqrtReturnType(derived()); }

EIGEN_DOC_UNARY_ADDONS(cwiseSign,sign function)

EIGEN_DEVICE_FUNC
inline const CwiseSignReturnType
cwiseSign() const { return CwiseSignReturnType(derived()); }

EIGEN_DOC_UNARY_ADDONS(cwiseInverse,inverse)

EIGEN_DEVICE_FUNC
inline const CwiseInverseReturnType
cwiseInverse() const { return CwiseInverseReturnType(derived()); }

EIGEN_DOC_UNARY_ADDONS(cwiseArg,arg)

EIGEN_DEVICE_FUNC
inline const CwiseArgReturnType
cwiseArg() const { return CwiseArgReturnType(derived()); }
