typedef CwiseUnaryOp<internal::scalar_abs_op<Scalar>, const Derived> AbsReturnType;
typedef CwiseUnaryOp<internal::scalar_arg_op<Scalar>, const Derived> ArgReturnType;
typedef CwiseUnaryOp<internal::scalar_abs2_op<Scalar>, const Derived> Abs2ReturnType;
typedef CwiseUnaryOp<internal::scalar_sqrt_op<Scalar>, const Derived> SqrtReturnType;
typedef CwiseUnaryOp<internal::scalar_rsqrt_op<Scalar>, const Derived> RsqrtReturnType;
typedef CwiseUnaryOp<internal::scalar_sign_op<Scalar>, const Derived> SignReturnType;
typedef CwiseUnaryOp<internal::scalar_inverse_op<Scalar>, const Derived> InverseReturnType;
typedef CwiseUnaryOp<internal::scalar_boolean_not_op<Scalar>, const Derived> BooleanNotReturnType;

typedef CwiseUnaryOp<internal::scalar_exp_op<Scalar>, const Derived> ExpReturnType;
typedef CwiseUnaryOp<internal::scalar_expm1_op<Scalar>, const Derived> Expm1ReturnType;
typedef CwiseUnaryOp<internal::scalar_log_op<Scalar>, const Derived> LogReturnType;
typedef CwiseUnaryOp<internal::scalar_log1p_op<Scalar>, const Derived> Log1pReturnType;
typedef CwiseUnaryOp<internal::scalar_log10_op<Scalar>, const Derived> Log10ReturnType;
typedef CwiseUnaryOp<internal::scalar_log2_op<Scalar>, const Derived> Log2ReturnType;
typedef CwiseUnaryOp<internal::scalar_cos_op<Scalar>, const Derived> CosReturnType;
typedef CwiseUnaryOp<internal::scalar_sin_op<Scalar>, const Derived> SinReturnType;
typedef CwiseUnaryOp<internal::scalar_tan_op<Scalar>, const Derived> TanReturnType;
typedef CwiseUnaryOp<internal::scalar_acos_op<Scalar>, const Derived> AcosReturnType;
typedef CwiseUnaryOp<internal::scalar_asin_op<Scalar>, const Derived> AsinReturnType;
typedef CwiseUnaryOp<internal::scalar_atan_op<Scalar>, const Derived> AtanReturnType;
typedef CwiseUnaryOp<internal::scalar_tanh_op<Scalar>, const Derived> TanhReturnType;
typedef CwiseUnaryOp<internal::scalar_logistic_op<Scalar>, const Derived> LogisticReturnType;
typedef CwiseUnaryOp<internal::scalar_sinh_op<Scalar>, const Derived> SinhReturnType;
#if EIGEN_HAS_CXX11_MATH
typedef CwiseUnaryOp<internal::scalar_atanh_op<Scalar>, const Derived> AtanhReturnType;
typedef CwiseUnaryOp<internal::scalar_asinh_op<Scalar>, const Derived> AsinhReturnType;
typedef CwiseUnaryOp<internal::scalar_acosh_op<Scalar>, const Derived> AcoshReturnType;
#endif
typedef CwiseUnaryOp<internal::scalar_cosh_op<Scalar>, const Derived> CoshReturnType;
typedef CwiseUnaryOp<internal::scalar_square_op<Scalar>, const Derived> SquareReturnType;
typedef CwiseUnaryOp<internal::scalar_cube_op<Scalar>, const Derived> CubeReturnType;
typedef CwiseUnaryOp<internal::scalar_round_op<Scalar>, const Derived> RoundReturnType;
typedef CwiseUnaryOp<internal::scalar_rint_op<Scalar>, const Derived> RintReturnType;
typedef CwiseUnaryOp<internal::scalar_floor_op<Scalar>, const Derived> FloorReturnType;
typedef CwiseUnaryOp<internal::scalar_ceil_op<Scalar>, const Derived> CeilReturnType;
typedef CwiseUnaryOp<internal::scalar_isnan_op<Scalar>, const Derived> IsNaNReturnType;
typedef CwiseUnaryOp<internal::scalar_isinf_op<Scalar>, const Derived> IsInfReturnType;
typedef CwiseUnaryOp<internal::scalar_isfinite_op<Scalar>, const Derived> IsFiniteReturnType;

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const AbsReturnType
abs() const
{
  return AbsReturnType(derived());
}

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const ArgReturnType
arg() const
{
  return ArgReturnType(derived());
}

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const Abs2ReturnType
abs2() const
{
  return Abs2ReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const ExpReturnType
exp() const
{
  return ExpReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const Expm1ReturnType
expm1() const
{
  return Expm1ReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const LogReturnType
log() const
{
  return LogReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const Log1pReturnType
log1p() const
{
  return Log1pReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const Log10ReturnType
log10() const
{
  return Log10ReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const Log2ReturnType
log2() const
{
  return Log2ReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const SqrtReturnType
sqrt() const
{
  return SqrtReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const RsqrtReturnType
rsqrt() const
{
  return RsqrtReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const SignReturnType
sign() const
{
  return SignReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const CosReturnType
cos() const
{
  return CosReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const SinReturnType
sin() const
{
  return SinReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const TanReturnType
tan() const
{
  return TanReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const AtanReturnType
atan() const
{
  return AtanReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const AcosReturnType
acos() const
{
  return AcosReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const AsinReturnType
asin() const
{
  return AsinReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const TanhReturnType
tanh() const
{
  return TanhReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const SinhReturnType
sinh() const
{
  return SinhReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const CoshReturnType
cosh() const
{
  return CoshReturnType(derived());
}

#if EIGEN_HAS_CXX11_MATH

EIGEN_DEVICE_FUNC
inline const AtanhReturnType
atanh() const
{
  return AtanhReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const AsinhReturnType
asinh() const
{
  return AsinhReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const AcoshReturnType
acosh() const
{
  return AcoshReturnType(derived());
}
#endif

EIGEN_DEVICE_FUNC
inline const LogisticReturnType
logistic() const
{
  return LogisticReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const InverseReturnType
inverse() const
{
  return InverseReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const SquareReturnType
square() const
{
  return SquareReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const CubeReturnType
cube() const
{
  return CubeReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const RintReturnType
rint() const
{
  return RintReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const RoundReturnType
round() const
{
  return RoundReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const FloorReturnType
floor() const
{
  return FloorReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const CeilReturnType
ceil() const
{
  return CeilReturnType(derived());
}

template<int N> struct ShiftRightXpr {
  typedef CwiseUnaryOp<internal::scalar_shift_right_op<Scalar, N>, const Derived> Type;
};

template<int N>
EIGEN_DEVICE_FUNC
typename ShiftRightXpr<N>::Type
shiftRight() const
{
  return typename ShiftRightXpr<N>::Type(derived());
}

template<int N> struct ShiftLeftXpr {
  typedef CwiseUnaryOp<internal::scalar_shift_left_op<Scalar, N>, const Derived> Type;
};

template<int N>
EIGEN_DEVICE_FUNC
typename ShiftLeftXpr<N>::Type
shiftLeft() const
{
  return typename ShiftLeftXpr<N>::Type(derived());
}

EIGEN_DEVICE_FUNC
inline const IsNaNReturnType
isNaN() const
{
  return IsNaNReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const IsInfReturnType
isInf() const
{
  return IsInfReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const IsFiniteReturnType
isFinite() const
{
  return IsFiniteReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const BooleanNotReturnType
operator!() const
{
  EIGEN_STATIC_ASSERT((internal::is_same<bool,Scalar>::value),
                      THIS_METHOD_IS_ONLY_FOR_EXPRESSIONS_OF_BOOL);
  return BooleanNotReturnType(derived());
}

typedef CwiseUnaryOp<internal::scalar_lgamma_op<Scalar>, const Derived> LgammaReturnType;
typedef CwiseUnaryOp<internal::scalar_digamma_op<Scalar>, const Derived> DigammaReturnType;
typedef CwiseUnaryOp<internal::scalar_erf_op<Scalar>, const Derived> ErfReturnType;
typedef CwiseUnaryOp<internal::scalar_erfc_op<Scalar>, const Derived> ErfcReturnType;
typedef CwiseUnaryOp<internal::scalar_ndtri_op<Scalar>, const Derived> NdtriReturnType;

EIGEN_DEVICE_FUNC
inline const LgammaReturnType
lgamma() const
{
  return LgammaReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const DigammaReturnType
digamma() const
{
  return DigammaReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const ErfReturnType
erf() const
{
  return ErfReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const ErfcReturnType
erfc() const
{
  return ErfcReturnType(derived());
}

EIGEN_DEVICE_FUNC
inline const NdtriReturnType
ndtri() const
{
  return NdtriReturnType(derived());
}
