#ifndef EIGEN_MATH_FUNCTIONS_SVE_H
#define EIGEN_MATH_FUNCTIONS_SVE_H

namespace Eigen {
namespace internal {

template <>
EIGEN_STRONG_INLINE EIGEN_UNUSED PacketXf pexp<PacketXf>(const PacketXf& x) {
  return pexp_float(x);
}

template <>
EIGEN_STRONG_INLINE EIGEN_UNUSED PacketXf plog<PacketXf>(const PacketXf& x) {
  return plog_float(x);
}

template <>
EIGEN_STRONG_INLINE EIGEN_UNUSED PacketXf psin<PacketXf>(const PacketXf& x) {
  return psin_float(x);
}

template <>
EIGEN_STRONG_INLINE EIGEN_UNUSED PacketXf pcos<PacketXf>(const PacketXf& x) {
  return pcos_float(x);
}

template <>
EIGEN_STRONG_INLINE EIGEN_UNUSED PacketXf ptanh<PacketXf>(const PacketXf& x) {
  return internal::generic_fast_tanh_float(x);
}
}  
}  

#endif  
