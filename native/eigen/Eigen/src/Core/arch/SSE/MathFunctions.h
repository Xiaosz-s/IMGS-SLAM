#ifndef EIGEN_MATH_FUNCTIONS_SSE_H
#define EIGEN_MATH_FUNCTIONS_SSE_H

namespace Eigen {

namespace internal {

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4f plog<Packet4f>(const Packet4f& _x) {
  return plog_float(_x);
}

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet2d plog<Packet2d>(const Packet2d& _x) {
  return plog_double(_x);
}

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4f plog2<Packet4f>(const Packet4f& _x) {
  return plog2_float(_x);
}

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet2d plog2<Packet2d>(const Packet2d& _x) {
  return plog2_double(_x);
}

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4f plog1p<Packet4f>(const Packet4f& _x) {
  return generic_plog1p(_x);
}

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4f pexpm1<Packet4f>(const Packet4f& _x) {
  return generic_expm1(_x);
}

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4f pexp<Packet4f>(const Packet4f& _x)
{
  return pexp_float(_x);
}

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet2d pexp<Packet2d>(const Packet2d& x)
{
  return pexp_double(x);
}

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4f psin<Packet4f>(const Packet4f& _x)
{
  return psin_float(_x);
}

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4f pcos<Packet4f>(const Packet4f& _x)
{
  return pcos_float(_x);
}

#if EIGEN_FAST_MATH

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4f psqrt<Packet4f>(const Packet4f& _x)
{
  Packet4f minus_half_x = pmul(_x, pset1<Packet4f>(-0.5f));
  Packet4f denormal_mask = pandnot(
      pcmp_lt(_x, pset1<Packet4f>((std::numeric_limits<float>::min)())),
      pcmp_lt(_x, pzero(_x)));

  Packet4f x = _mm_rsqrt_ps(_x);
  
  x = pmul(x, pmadd(minus_half_x, pmul(x,x), pset1<Packet4f>(1.5f)));
  
  return pandnot(pmul(_x,x), denormal_mask);
}

#else

template<>EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4f psqrt<Packet4f>(const Packet4f& x) { return _mm_sqrt_ps(x); }

#endif

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet2d psqrt<Packet2d>(const Packet2d& x) { return _mm_sqrt_pd(x); }

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet16b psqrt<Packet16b>(const Packet16b& x) { return x; }

#if EIGEN_FAST_MATH

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4f prsqrt<Packet4f>(const Packet4f& _x) {
  _EIGEN_DECLARE_CONST_Packet4f(one_point_five, 1.5f);
  _EIGEN_DECLARE_CONST_Packet4f(minus_half, -0.5f);
  _EIGEN_DECLARE_CONST_Packet4f_FROM_INT(inf, 0x7f800000u);
  _EIGEN_DECLARE_CONST_Packet4f_FROM_INT(flt_min, 0x00800000u);

  Packet4f neg_half = pmul(_x, p4f_minus_half);

  Packet4f lt_min_mask = _mm_cmplt_ps(_x, p4f_flt_min);
  Packet4f inf_mask = _mm_cmpeq_ps(_x, p4f_inf);
  Packet4f not_normal_finite_mask = _mm_or_ps(lt_min_mask, inf_mask);

  Packet4f y_approx = _mm_rsqrt_ps(_x);

  Packet4f y_newton = pmul(
      y_approx, pmadd(y_approx, pmul(neg_half, y_approx), p4f_one_point_five));

  return pselect<Packet4f>(not_normal_finite_mask, y_approx, y_newton);
}

#else

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet4f prsqrt<Packet4f>(const Packet4f& x) {
  
  return _mm_div_ps(pset1<Packet4f>(1.0f), _mm_sqrt_ps(x));
}

#endif

template<> EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED
Packet2d prsqrt<Packet2d>(const Packet2d& x) {
  return _mm_div_pd(pset1<Packet2d>(1.0), _mm_sqrt_pd(x));
}

template <>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED Packet4f
ptanh<Packet4f>(const Packet4f& x) {
  return internal::generic_fast_tanh_float(x);
}

} 

namespace numext {

template<>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE
float sqrt(const float &x)
{
  return internal::pfirst(internal::Packet4f(_mm_sqrt_ss(_mm_set_ss(x))));
}

template<>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE
double sqrt(const double &x)
{
#if EIGEN_COMP_GNUC_STRICT

  return internal::pfirst(internal::Packet2d(__builtin_ia32_sqrtsd(_mm_set_sd(x))));
#else
  return internal::pfirst(internal::Packet2d(_mm_sqrt_pd(_mm_set_sd(x))));
#endif
}

} 

} 

#endif 
