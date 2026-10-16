#ifndef EIGEN_MATH_FUNCTIONS_AVX_H
#define EIGEN_MATH_FUNCTIONS_AVX_H

namespace Eigen {

namespace internal {

template <>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED Packet8f
psin<Packet8f>(const Packet8f& _x) {
  return psin_float(_x);
}

template <>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED Packet8f
pcos<Packet8f>(const Packet8f& _x) {
  return pcos_float(_x);
}

template <>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED Packet8f
plog<Packet8f>(const Packet8f& _x) {
  return plog_float(_x);
}

template <>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED Packet4d
plog<Packet4d>(const Packet4d& _x) {
  return plog_double(_x);
}

template <>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED Packet8f
plog2<Packet8f>(const Packet8f& _x) {
  return plog2_float(_x);
}

template <>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED Packet4d
plog2<Packet4d>(const Packet4d& _x) {
  return plog2_double(_x);
}

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet8f plog1p<Packet8f>(const Packet8f& _x) {
  return generic_plog1p(_x);
}

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet8f pexpm1<Packet8f>(const Packet8f& _x) {
  return generic_expm1(_x);
}

template <>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED Packet8f
pexp<Packet8f>(const Packet8f& _x) {
  return pexp_float(_x);
}

template <>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED Packet8f
ptanh<Packet8f>(const Packet8f& _x) {
  return internal::generic_fast_tanh_float(_x);
}

template <>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED Packet4d
pexp<Packet4d>(const Packet4d& _x) {
  return pexp_double(_x);
}

#if EIGEN_FAST_MATH
template <>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet8f psqrt<Packet8f>(const Packet8f& _x) {
  Packet8f minus_half_x = pmul(_x, pset1<Packet8f>(-0.5f));
  Packet8f denormal_mask = pandnot(
      pcmp_lt(_x, pset1<Packet8f>((std::numeric_limits<float>::min)())),
      pcmp_lt(_x, pzero(_x)));

  Packet8f x = _mm256_rsqrt_ps(_x);
  
  x = pmul(x, pmadd(minus_half_x, pmul(x,x), pset1<Packet8f>(1.5f)));
  
  return pandnot(pmul(_x,x), denormal_mask);
}

#else

template <> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet8f psqrt<Packet8f>(const Packet8f& _x) {
  return _mm256_sqrt_ps(_x);
}

#endif

template <> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4d psqrt<Packet4d>(const Packet4d& _x) {
  return _mm256_sqrt_pd(_x);
}

#if EIGEN_FAST_MATH
template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet8f prsqrt<Packet8f>(const Packet8f& _x) {
  _EIGEN_DECLARE_CONST_Packet8f_FROM_INT(inf, 0x7f800000);
  _EIGEN_DECLARE_CONST_Packet8f(one_point_five, 1.5f);
  _EIGEN_DECLARE_CONST_Packet8f(minus_half, -0.5f);
  _EIGEN_DECLARE_CONST_Packet8f_FROM_INT(flt_min, 0x00800000);

  Packet8f neg_half = pmul(_x, p8f_minus_half);

  Packet8f lt_min_mask = _mm256_cmp_ps(_x, p8f_flt_min, _CMP_LT_OQ);
  Packet8f inf_mask =  _mm256_cmp_ps(_x, p8f_inf, _CMP_EQ_OQ);
  Packet8f not_normal_finite_mask = _mm256_or_ps(lt_min_mask, inf_mask);

  Packet8f y_approx = _mm256_rsqrt_ps(_x);

  Packet8f y_newton = pmul(y_approx, pmadd(y_approx, pmul(neg_half, y_approx), p8f_one_point_five));

  return pselect<Packet8f>(not_normal_finite_mask, y_approx, y_newton);
}

#else
template <> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet8f prsqrt<Packet8f>(const Packet8f& _x) {
  _EIGEN_DECLARE_CONST_Packet8f(one, 1.0f);
  return _mm256_div_ps(p8f_one, _mm256_sqrt_ps(_x));
}
#endif

template <> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4d prsqrt<Packet4d>(const Packet4d& _x) {
  _EIGEN_DECLARE_CONST_Packet4d(one, 1.0);
  return _mm256_div_pd(p4d_one, _mm256_sqrt_pd(_x));
}

F16_PACKET_FUNCTION(Packet8f, Packet8h, psin)
F16_PACKET_FUNCTION(Packet8f, Packet8h, pcos)
F16_PACKET_FUNCTION(Packet8f, Packet8h, plog)
F16_PACKET_FUNCTION(Packet8f, Packet8h, plog2)
F16_PACKET_FUNCTION(Packet8f, Packet8h, plog1p)
F16_PACKET_FUNCTION(Packet8f, Packet8h, pexpm1)
F16_PACKET_FUNCTION(Packet8f, Packet8h, pexp)
F16_PACKET_FUNCTION(Packet8f, Packet8h, ptanh)
F16_PACKET_FUNCTION(Packet8f, Packet8h, psqrt)
F16_PACKET_FUNCTION(Packet8f, Packet8h, prsqrt)

template <>
EIGEN_STRONG_INLINE Packet8h pfrexp(const Packet8h& a, Packet8h& exponent) {
  Packet8f fexponent;
  const Packet8h out = float2half(pfrexp<Packet8f>(half2float(a), fexponent));
  exponent = float2half(fexponent);
  return out;
}

template <>
EIGEN_STRONG_INLINE Packet8h pldexp(const Packet8h& a, const Packet8h& exponent) {
  return float2half(pldexp<Packet8f>(half2float(a), half2float(exponent)));
}

BF16_PACKET_FUNCTION(Packet8f, Packet8bf, psin)
BF16_PACKET_FUNCTION(Packet8f, Packet8bf, pcos)
BF16_PACKET_FUNCTION(Packet8f, Packet8bf, plog)
BF16_PACKET_FUNCTION(Packet8f, Packet8bf, plog2)
BF16_PACKET_FUNCTION(Packet8f, Packet8bf, plog1p)
BF16_PACKET_FUNCTION(Packet8f, Packet8bf, pexpm1)
BF16_PACKET_FUNCTION(Packet8f, Packet8bf, pexp)
BF16_PACKET_FUNCTION(Packet8f, Packet8bf, ptanh)
BF16_PACKET_FUNCTION(Packet8f, Packet8bf, psqrt)
BF16_PACKET_FUNCTION(Packet8f, Packet8bf, prsqrt)

template <>
EIGEN_STRONG_INLINE Packet8bf pfrexp(const Packet8bf& a, Packet8bf& exponent) {
  Packet8f fexponent;
  const Packet8bf out = F32ToBf16(pfrexp<Packet8f>(Bf16ToF32(a), fexponent));
  exponent = F32ToBf16(fexponent);
  return out;
}

template <>
EIGEN_STRONG_INLINE Packet8bf pldexp(const Packet8bf& a, const Packet8bf& exponent) {
  return F32ToBf16(pldexp<Packet8f>(Bf16ToF32(a), Bf16ToF32(exponent)));
}

}  

}  

#endif  
