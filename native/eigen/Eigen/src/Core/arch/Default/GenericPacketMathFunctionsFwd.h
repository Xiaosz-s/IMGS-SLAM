#ifndef EIGEN_ARCH_GENERIC_PACKET_MATH_FUNCTIONS_FWD_H
#define EIGEN_ARCH_GENERIC_PACKET_MATH_FUNCTIONS_FWD_H

namespace Eigen {
namespace internal {

template<typename Packet> EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC
Packet pfrexp_generic(const Packet& a, Packet& exponent);

template<typename Packet> EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC
Packet pfrexp_generic_get_biased_exponent(const Packet& p);

template<typename Packet> EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC
Packet pldexp_generic(const Packet& a, const Packet& exponent);

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet plog_float(const Packet _x);

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet plog2_float(const Packet _x);

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet plog_double(const Packet _x);

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet plog2_double(const Packet _x);

template<typename Packet>
Packet generic_plog1p(const Packet& x);

template<typename Packet>
Packet generic_expm1(const Packet& x);

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet pexp_float(const Packet _x);

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet pexp_double(const Packet _x);

template<typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet psin_float(const Packet& x);

template<typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet pcos_float(const Packet& x);

template<typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet psqrt_complex(const Packet& a);

template<typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet pdiv_complex(const Packet& x, const Packet& y);

template <typename Packet, int N> struct ppolevl;

} 
} 

#endif 
