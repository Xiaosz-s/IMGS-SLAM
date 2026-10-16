#ifndef EIGEN_COMPLEX_CUDA_H
#define EIGEN_COMPLEX_CUDA_H

#if defined(EIGEN_CUDACC) && defined(EIGEN_GPU_COMPILE_PHASE)

#if !(defined(EIGEN_COMP_ICC) && defined(_USE_COMPLEX_SPECIALIZATION_))

#define EIGEN_USING_STD_COMPLEX_OPERATORS           \
  using Eigen::complex_operator_detail::operator+;  \
  using Eigen::complex_operator_detail::operator-;  \
  using Eigen::complex_operator_detail::operator*;  \
  using Eigen::complex_operator_detail::operator/;  \
  using Eigen::complex_operator_detail::operator+=; \
  using Eigen::complex_operator_detail::operator-=; \
  using Eigen::complex_operator_detail::operator*=; \
  using Eigen::complex_operator_detail::operator/=; \
  using Eigen::complex_operator_detail::operator==; \
  using Eigen::complex_operator_detail::operator!=;

namespace Eigen {

namespace complex_operator_detail {

template<typename T>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
std::complex<T> complex_multiply(const std::complex<T>& a, const std::complex<T>& b) {
  const T a_real = numext::real(a);
  const T a_imag = numext::imag(a);
  const T b_real = numext::real(b);
  const T b_imag = numext::imag(b);
  return std::complex<T>(
      a_real * b_real - a_imag * b_imag,
      a_imag * b_real + a_real * b_imag);
}

template<typename T>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
std::complex<T> complex_divide_fast(const std::complex<T>& a, const std::complex<T>& b) {
  const T a_real = numext::real(a);
  const T a_imag = numext::imag(a);
  const T b_real = numext::real(b);
  const T b_imag = numext::imag(b);
  const T norm = (b_real * b_real + b_imag * b_imag);
  return std::complex<T>((a_real * b_real + a_imag * b_imag) / norm,
                          (a_imag * b_real - a_real * b_imag) / norm);
}

template<typename T>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
std::complex<T> complex_divide_stable(const std::complex<T>& a, const std::complex<T>& b) {
  const T a_real = numext::real(a);
  const T a_imag = numext::imag(a);
  const T b_real = numext::real(b);
  const T b_imag = numext::imag(b);

  const bool scale_imag = numext::abs(b_imag) <= numext::abs(b_real);
  const T rscale = scale_imag ? T(1) : b_real / b_imag;
  const T iscale = scale_imag ? b_imag / b_real : T(1);
  const T denominator = b_real * rscale + b_imag * iscale;
  return std::complex<T>((a_real * rscale + a_imag * iscale) / denominator, 
                         (a_imag * rscale - a_real * iscale) / denominator);
}

template<typename T>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
std::complex<T> complex_divide(const std::complex<T>& a, const std::complex<T>& b) {
#if EIGEN_FAST_MATH
  return complex_divide_fast(a, b);
#else
  return complex_divide_stable(a, b);
#endif
}

#define EIGEN_CREATE_STD_COMPLEX_OPERATOR_SPECIALIZATIONS(T)                                    \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator+(const std::complex<T>& a) { return a; }                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator-(const std::complex<T>& a) {                                           \
  return std::complex<T>(-numext::real(a), -numext::imag(a));                                   \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator+(const std::complex<T>& a, const std::complex<T>& b) {                 \
  return std::complex<T>(numext::real(a) + numext::real(b), numext::imag(a) + numext::imag(b)); \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator+(const std::complex<T>& a, const T& b) {                               \
  return std::complex<T>(numext::real(a) + b, numext::imag(a));                                 \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator+(const T& a, const std::complex<T>& b) {                               \
  return std::complex<T>(a + numext::real(b), numext::imag(b));                                 \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator-(const std::complex<T>& a, const std::complex<T>& b) {                 \
  return std::complex<T>(numext::real(a) - numext::real(b), numext::imag(a) - numext::imag(b)); \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator-(const std::complex<T>& a, const T& b) {                               \
  return std::complex<T>(numext::real(a) - b, numext::imag(a));                                 \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator-(const T& a, const std::complex<T>& b) {                               \
  return std::complex<T>(a - numext::real(b), -numext::imag(b));                                \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator*(const std::complex<T>& a, const std::complex<T>& b) {                 \
  return complex_multiply(a, b);                                                                \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator*(const std::complex<T>& a, const T& b) {                               \
  return std::complex<T>(numext::real(a) * b, numext::imag(a) * b);                             \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator*(const T& a, const std::complex<T>& b) {                               \
  return std::complex<T>(a * numext::real(b), a * numext::imag(b));                             \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator/(const std::complex<T>& a, const std::complex<T>& b) {                 \
  return complex_divide(a, b);                                                                  \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator/(const std::complex<T>& a, const T& b) {                               \
  return std::complex<T>(numext::real(a) / b, numext::imag(a) / b);                             \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T> operator/(const T& a, const std::complex<T>& b) {                               \
  return complex_divide(std::complex<T>(a, 0), b);                                              \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T>& operator+=(std::complex<T>& a, const std::complex<T>& b) {                     \
  numext::real_ref(a) += numext::real(b);                                                       \
  numext::imag_ref(a) += numext::imag(b);                                                       \
  return a;                                                                                     \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T>& operator-=(std::complex<T>& a, const std::complex<T>& b) {                     \
  numext::real_ref(a) -= numext::real(b);                                                       \
  numext::imag_ref(a) -= numext::imag(b);                                                       \
  return a;                                                                                     \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T>& operator*=(std::complex<T>& a, const std::complex<T>& b) {                     \
  a = complex_multiply(a, b);                                                                   \
  return a;                                                                                     \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
std::complex<T>& operator/=(std::complex<T>& a, const std::complex<T>& b) {                     \
  a = complex_divide(a, b);                                                                     \
  return  a;                                                                                    \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
bool operator==(const std::complex<T>& a, const std::complex<T>& b) {                           \
  return numext::real(a) == numext::real(b) && numext::imag(a) == numext::imag(b);              \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
bool operator==(const std::complex<T>& a, const T& b) {                                         \
  return numext::real(a) == b && numext::imag(a) == 0;                                          \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
bool operator==(const T& a, const std::complex<T>& b) {                                         \
  return a  == numext::real(b) && 0 == numext::imag(b);                                         \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
bool operator!=(const std::complex<T>& a, const std::complex<T>& b) {                           \
  return !(a == b);                                                                             \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
bool operator!=(const std::complex<T>& a, const T& b) {                                         \
  return !(a == b);                                                                             \
}                                                                                               \
                                                                                                \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                           \
bool operator!=(const T& a, const std::complex<T>& b) {                                         \
  return !(a == b);                                                                             \
}

EIGEN_CREATE_STD_COMPLEX_OPERATOR_SPECIALIZATIONS(float)
EIGEN_CREATE_STD_COMPLEX_OPERATOR_SPECIALIZATIONS(double)

#undef EIGEN_CREATE_STD_COMPLEX_OPERATOR_SPECIALIZATIONS

}  

EIGEN_USING_STD_COMPLEX_OPERATORS

namespace numext {
EIGEN_USING_STD_COMPLEX_OPERATORS
}  

namespace internal {
EIGEN_USING_STD_COMPLEX_OPERATORS

}  
}  

#endif  

#endif  

#endif  
