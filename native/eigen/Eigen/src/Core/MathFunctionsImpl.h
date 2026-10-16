#ifndef EIGEN_MATHFUNCTIONSIMPL_H
#define EIGEN_MATHFUNCTIONSIMPL_H

namespace Eigen {

namespace internal {

template<typename T>
T generic_fast_tanh_float(const T& a_x)
{
  
#ifdef EIGEN_VECTORIZE_FMA
  const T plus_clamp = pset1<T>(7.99881172180175781f);
  const T minus_clamp = pset1<T>(-7.99881172180175781f);
#else
  const T plus_clamp = pset1<T>(7.90531110763549805f);
  const T minus_clamp = pset1<T>(-7.90531110763549805f);
#endif
  const T tiny = pset1<T>(0.0004f);
  const T x = pmax(pmin(a_x, plus_clamp), minus_clamp);
  const T tiny_mask = pcmp_lt(pabs(a_x), tiny);
  
  const T alpha_1 = pset1<T>(4.89352455891786e-03f);
  const T alpha_3 = pset1<T>(6.37261928875436e-04f);
  const T alpha_5 = pset1<T>(1.48572235717979e-05f);
  const T alpha_7 = pset1<T>(5.12229709037114e-08f);
  const T alpha_9 = pset1<T>(-8.60467152213735e-11f);
  const T alpha_11 = pset1<T>(2.00018790482477e-13f);
  const T alpha_13 = pset1<T>(-2.76076847742355e-16f);

  const T beta_0 = pset1<T>(4.89352518554385e-03f);
  const T beta_2 = pset1<T>(2.26843463243900e-03f);
  const T beta_4 = pset1<T>(1.18534705686654e-04f);
  const T beta_6 = pset1<T>(1.19825839466702e-06f);

  const T x2 = pmul(x, x);

  T p = pmadd(x2, alpha_13, alpha_11);
  p = pmadd(x2, p, alpha_9);
  p = pmadd(x2, p, alpha_7);
  p = pmadd(x2, p, alpha_5);
  p = pmadd(x2, p, alpha_3);
  p = pmadd(x2, p, alpha_1);
  p = pmul(x, p);

  T q = pmadd(x2, beta_6, beta_4);
  q = pmadd(x2, q, beta_2);
  q = pmadd(x2, q, beta_0);

  return pselect(tiny_mask, x, pdiv(p, q));
}

template<typename RealScalar>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
RealScalar positive_real_hypot(const RealScalar& x, const RealScalar& y)
{
  
  if ((numext::isinf)(x) || (numext::isinf)(y))
    return NumTraits<RealScalar>::infinity();
  if ((numext::isnan)(x) || (numext::isnan)(y))
    return NumTraits<RealScalar>::quiet_NaN();
    
  EIGEN_USING_STD(sqrt);
  RealScalar p, qp;
  p = numext::maxi(x,y);
  if(p==RealScalar(0)) return RealScalar(0);
  qp = numext::mini(y,x) / p;
  return p * sqrt(RealScalar(1) + qp*qp);
}

template<typename Scalar>
struct hypot_impl
{
  typedef typename NumTraits<Scalar>::Real RealScalar;
  static EIGEN_DEVICE_FUNC
  inline RealScalar run(const Scalar& x, const Scalar& y)
  {
    EIGEN_USING_STD(abs);
    return positive_real_hypot<RealScalar>(abs(x), abs(y));
  }
};

template<typename T>
EIGEN_DEVICE_FUNC std::complex<T> complex_sqrt(const std::complex<T>& z) {

  const T x = numext::real(z);
  const T y = numext::imag(z);
  const T zero = T(0);
  const T w = numext::sqrt(T(0.5) * (numext::abs(x) + numext::hypot(x, y)));

  return
    (numext::isinf)(y) ? std::complex<T>(NumTraits<T>::infinity(), y)
      : x == zero ? std::complex<T>(w, y < zero ? -w : w)
      : x > zero ? std::complex<T>(w, y / (2 * w))
      : std::complex<T>(numext::abs(y) / (2 * w), y < zero ? -w : w );
}

template<typename T>
EIGEN_DEVICE_FUNC std::complex<T> complex_rsqrt(const std::complex<T>& z) {

  const T x = numext::real(z);
  const T y = numext::imag(z);
  const T zero = T(0);

  const T abs_z = numext::hypot(x, y);
  const T w = numext::sqrt(T(0.5) * (numext::abs(x) + abs_z));
  const T woz = w / abs_z;
  
  return
    abs_z == zero ? std::complex<T>(NumTraits<T>::infinity(), NumTraits<T>::quiet_NaN())
      : ((numext::isinf)(x) || (numext::isinf)(y)) ? std::complex<T>(zero, zero)
      : x == zero ? std::complex<T>(woz, y < zero ? woz : -woz)
      : x > zero ? std::complex<T>(woz, -y / (2 * w * abs_z))
      : std::complex<T>(numext::abs(y) / (2 * w * abs_z), y < zero ? woz : -woz );
}

template<typename T>
EIGEN_DEVICE_FUNC std::complex<T> complex_log(const std::complex<T>& z) {
  
  T a = numext::abs(z);
  EIGEN_USING_STD(atan2);
  T b = atan2(z.imag(), z.real());
  return std::complex<T>(numext::log(a), b);
}

} 

} 

#endif 
