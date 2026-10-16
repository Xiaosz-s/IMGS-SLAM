#ifndef EIGEN_HALF_H
#define EIGEN_HALF_H

#include <sstream>

#if defined(EIGEN_HAS_GPU_FP16) || defined(EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC)

  #pragma push_macro("EIGEN_CONSTEXPR")
  #undef EIGEN_CONSTEXPR
  #define EIGEN_CONSTEXPR
#endif

#define F16_PACKET_FUNCTION(PACKET_F, PACKET_F16, METHOD)           \
  template <>                                                       \
  EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC EIGEN_UNUSED                \
  PACKET_F16 METHOD<PACKET_F16>(const PACKET_F16& _x) {             \
    return float2half(METHOD<PACKET_F>(half2float(_x)));            \
  }

namespace Eigen {

struct half;

namespace half_impl {

#if !defined(EIGEN_HAS_GPU_FP16) || !defined(EIGEN_GPU_COMPILE_PHASE)

struct __half_raw {
#if (defined(EIGEN_HAS_GPU_FP16) && !defined(EIGEN_GPU_COMPILE_PHASE))

  EIGEN_DEVICE_FUNC __half_raw() {}
#else
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR __half_raw() : x(0) {}
#endif
#if defined(EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC)
  explicit EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR __half_raw(numext::uint16_t raw) : x(numext::bit_cast<__fp16>(raw)) {
  }
  __fp16 x;
#else
  explicit EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR __half_raw(numext::uint16_t raw) : x(raw) {}
  numext::uint16_t x;
#endif
};

#elif defined(EIGEN_HAS_HIP_FP16)

#elif defined(EIGEN_HAS_CUDA_FP16)
  #if EIGEN_CUDA_SDK_VER < 90000
    
    typedef __half __half_raw;
  #endif 
#elif defined(SYCL_DEVICE_ONLY)
  typedef cl::sycl::half __half_raw;
#endif

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR __half_raw raw_uint16_to_half(numext::uint16_t x);
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC __half_raw float_to_half_rtne(float ff);
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC float half_to_float(__half_raw h);

struct half_base : public __half_raw {
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR half_base() {}
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR half_base(const __half_raw& h) : __half_raw(h) {}

#if defined(EIGEN_HAS_GPU_FP16)
 #if defined(EIGEN_HAS_HIP_FP16)
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR half_base(const __half& h) { x = __half_as_ushort(h); }
 #elif defined(EIGEN_HAS_CUDA_FP16)
  #if EIGEN_CUDA_SDK_VER >= 90000
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR half_base(const __half& h) : __half_raw(*(__half_raw*)&h) {}
  #endif
 #endif
#endif
};

} 

struct half : public half_impl::half_base {

#if !defined(EIGEN_HAS_GPU_FP16) || !defined(EIGEN_GPU_COMPILE_PHASE)

  typedef half_impl::__half_raw __half_raw;
#elif defined(EIGEN_HAS_HIP_FP16)

#elif defined(EIGEN_HAS_CUDA_FP16)

  #if defined(EIGEN_CUDA_SDK_VER) && EIGEN_CUDA_SDK_VER < 90000
    typedef half_impl::__half_raw __half_raw;
  #endif
#endif

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR half() {}

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR half(const __half_raw& h) : half_impl::half_base(h) {}

#if defined(EIGEN_HAS_GPU_FP16)
 #if defined(EIGEN_HAS_HIP_FP16)
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR half(const __half& h) : half_impl::half_base(h) {}
 #elif defined(EIGEN_HAS_CUDA_FP16)
  #if defined(EIGEN_CUDA_SDK_VER) && EIGEN_CUDA_SDK_VER >= 90000
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR half(const __half& h) : half_impl::half_base(h) {}
  #endif
 #endif
#endif

  explicit EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR half(bool b)
      : half_impl::half_base(half_impl::raw_uint16_to_half(b ? 0x3c00 : 0)) {}
  template<class T>
  explicit EIGEN_DEVICE_FUNC half(T val)
      : half_impl::half_base(half_impl::float_to_half_rtne(static_cast<float>(val))) {}
  explicit EIGEN_DEVICE_FUNC half(float f)
      : half_impl::half_base(half_impl::float_to_half_rtne(f)) {}

  template<typename RealScalar>
  explicit EIGEN_DEVICE_FUNC half(std::complex<RealScalar> c)
      : half_impl::half_base(half_impl::float_to_half_rtne(static_cast<float>(c.real()))) {}

   EIGEN_DEVICE_FUNC operator float() const {  
    return half_impl::half_to_float(*this);
  }

#if defined(EIGEN_HAS_GPU_FP16) && !defined(EIGEN_GPU_COMPILE_PHASE)
  EIGEN_DEVICE_FUNC operator __half() const {
    ::__half_raw hr;
    hr.x = x;
    return __half(hr);
  }
#endif
};

} 

namespace std {
template<>
struct numeric_limits<Eigen::half> {
  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool has_infinity = true;
  static const bool has_quiet_NaN = true;
  static const bool has_signaling_NaN = true;
  static const float_denorm_style has_denorm = denorm_present;
  static const bool has_denorm_loss = false;
  static const std::float_round_style round_style = std::round_to_nearest;
  static const bool is_iec559 = false;
  static const bool is_bounded = false;
  static const bool is_modulo = false;
  static const int digits = 11;
  static const int digits10 = 3;      
  static const int max_digits10 = 5;  
  static const int radix = 2;
  static const int min_exponent = -13;
  static const int min_exponent10 = -4;
  static const int max_exponent = 16;
  static const int max_exponent10 = 4;
  static const bool traps = true;
  static const bool tinyness_before = false;

  static Eigen::half (min)() { return Eigen::half_impl::raw_uint16_to_half(0x400); }
  static Eigen::half lowest() { return Eigen::half_impl::raw_uint16_to_half(0xfbff); }
  static Eigen::half (max)() { return Eigen::half_impl::raw_uint16_to_half(0x7bff); }
  static Eigen::half epsilon() { return Eigen::half_impl::raw_uint16_to_half(0x0800); }
  static Eigen::half round_error() { return Eigen::half(0.5); }
  static Eigen::half infinity() { return Eigen::half_impl::raw_uint16_to_half(0x7c00); }
  static Eigen::half quiet_NaN() { return Eigen::half_impl::raw_uint16_to_half(0x7e00); }
  static Eigen::half signaling_NaN() { return Eigen::half_impl::raw_uint16_to_half(0x7d00); }
  static Eigen::half denorm_min() { return Eigen::half_impl::raw_uint16_to_half(0x1); }
};

template<>
struct numeric_limits<const Eigen::half> : numeric_limits<Eigen::half> {};
template<>
struct numeric_limits<volatile Eigen::half> : numeric_limits<Eigen::half> {};
template<>
struct numeric_limits<const volatile Eigen::half> : numeric_limits<Eigen::half> {};
} 

namespace Eigen {

namespace half_impl {

#if (defined(EIGEN_HAS_CUDA_FP16) && defined(EIGEN_CUDA_ARCH) && \
     EIGEN_CUDA_ARCH >= 530) ||                                  \
    (defined(EIGEN_HAS_HIP_FP16) && defined(HIP_DEVICE_COMPILE))

#define EIGEN_HAS_NATIVE_FP16
#endif

#if defined(EIGEN_HAS_NATIVE_FP16)
EIGEN_STRONG_INLINE __device__ half operator + (const half& a, const half& b) {
#if defined(EIGEN_CUDA_SDK_VER) && EIGEN_CUDA_SDK_VER >= 90000
  return __hadd(::__half(a), ::__half(b));
#else
  return __hadd(a, b);
#endif
}
EIGEN_STRONG_INLINE __device__ half operator * (const half& a, const half& b) {
  return __hmul(a, b);
}
EIGEN_STRONG_INLINE __device__ half operator - (const half& a, const half& b) {
  return __hsub(a, b);
}
EIGEN_STRONG_INLINE __device__ half operator / (const half& a, const half& b) {
#if defined(EIGEN_CUDA_SDK_VER) && EIGEN_CUDA_SDK_VER >= 90000
  return __hdiv(a, b);
#else
  float num = __half2float(a);
  float denom = __half2float(b);
  return __float2half(num / denom);
#endif
}
EIGEN_STRONG_INLINE __device__ half operator - (const half& a) {
  return __hneg(a);
}
EIGEN_STRONG_INLINE __device__ half& operator += (half& a, const half& b) {
  a = a + b;
  return a;
}
EIGEN_STRONG_INLINE __device__ half& operator *= (half& a, const half& b) {
  a = a * b;
  return a;
}
EIGEN_STRONG_INLINE __device__ half& operator -= (half& a, const half& b) {
  a = a - b;
  return a;
}
EIGEN_STRONG_INLINE __device__ half& operator /= (half& a, const half& b) {
  a = a / b;
  return a;
}
EIGEN_STRONG_INLINE __device__ bool operator == (const half& a, const half& b) {
  return __heq(a, b);
}
EIGEN_STRONG_INLINE __device__ bool operator != (const half& a, const half& b) {
  return __hne(a, b);
}
EIGEN_STRONG_INLINE __device__ bool operator < (const half& a, const half& b) {
  return __hlt(a, b);
}
EIGEN_STRONG_INLINE __device__ bool operator <= (const half& a, const half& b) {
  return __hle(a, b);
}
EIGEN_STRONG_INLINE __device__ bool operator > (const half& a, const half& b) {
  return __hgt(a, b);
}
EIGEN_STRONG_INLINE __device__ bool operator >= (const half& a, const half& b) {
  return __hge(a, b);
}
#endif

#if defined(EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC)
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator + (const half& a, const half& b) {
  return half(vaddh_f16(a.x, b.x));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator * (const half& a, const half& b) {
  return half(vmulh_f16(a.x, b.x));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator - (const half& a, const half& b) {
  return half(vsubh_f16(a.x, b.x));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator / (const half& a, const half& b) {
  return half(vdivh_f16(a.x, b.x));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator - (const half& a) {
  return half(vnegh_f16(a.x));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half& operator += (half& a, const half& b) {
  a = half(vaddh_f16(a.x, b.x));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half& operator *= (half& a, const half& b) {
  a = half(vmulh_f16(a.x, b.x));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half& operator -= (half& a, const half& b) {
  a = half(vsubh_f16(a.x, b.x));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half& operator /= (half& a, const half& b) {
  a = half(vdivh_f16(a.x, b.x));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator == (const half& a, const half& b) {
  return vceqh_f16(a.x, b.x);
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator != (const half& a, const half& b) {
  return !vceqh_f16(a.x, b.x);
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator < (const half& a, const half& b) {
  return vclth_f16(a.x, b.x);
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator <= (const half& a, const half& b) {
  return vcleh_f16(a.x, b.x);
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator > (const half& a, const half& b) {
  return vcgth_f16(a.x, b.x);
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator >= (const half& a, const half& b) {
  return vcgeh_f16(a.x, b.x);
}

#elif !defined(EIGEN_HAS_NATIVE_FP16) || (EIGEN_COMP_CLANG && !EIGEN_COMP_NVCC) 

#if EIGEN_COMP_CLANG && defined(EIGEN_CUDACC)

#pragma push_macro("EIGEN_DEVICE_FUNC")
#undef EIGEN_DEVICE_FUNC
#if defined(EIGEN_HAS_CUDA_FP16) && defined(EIGEN_HAS_NATIVE_FP16)
#define EIGEN_DEVICE_FUNC __host__
#else 
#define EIGEN_DEVICE_FUNC __host__ __device__
#endif
#endif

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator + (const half& a, const half& b) {
  return half(float(a) + float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator * (const half& a, const half& b) {
  return half(float(a) * float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator - (const half& a, const half& b) {
  return half(float(a) - float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator / (const half& a, const half& b) {
  return half(float(a) / float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator - (const half& a) {
  half result;
  result.x = a.x ^ 0x8000;
  return result;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half& operator += (half& a, const half& b) {
  a = half(float(a) + float(b));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half& operator *= (half& a, const half& b) {
  a = half(float(a) * float(b));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half& operator -= (half& a, const half& b) {
  a = half(float(a) - float(b));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half& operator /= (half& a, const half& b) {
  a = half(float(a) / float(b));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator == (const half& a, const half& b) {
  return numext::equal_strict(float(a),float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator != (const half& a, const half& b) {
  return numext::not_equal_strict(float(a), float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator < (const half& a, const half& b) {
  return float(a) < float(b);
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator <= (const half& a, const half& b) {
  return float(a) <= float(b);
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator > (const half& a, const half& b) {
  return float(a) > float(b);
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator >= (const half& a, const half& b) {
  return float(a) >= float(b);
}

#if defined(__clang__) && defined(__CUDA__)
#pragma pop_macro("EIGEN_DEVICE_FUNC")
#endif
#endif  

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator / (const half& a, Index b) {
  return half(static_cast<float>(a) / static_cast<float>(b));
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator++(half& a) {
  a += half(1);
  return a;
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator--(half& a) {
  a -= half(1);
  return a;
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator++(half& a, int) {
  half original_value = a;
  ++a;
  return original_value;
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half operator--(half& a, int) {
  half original_value = a;
  --a;
  return original_value;
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR __half_raw raw_uint16_to_half(numext::uint16_t x) {

#if defined(EIGEN_HAS_GPU_FP16)
   __half_raw h;
   h.x = x;
  return h;
#else
  return __half_raw(x);
#endif
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC numext::uint16_t raw_half_as_uint16(const __half_raw& h) {

#if defined(EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC)
  return numext::bit_cast<numext::uint16_t>(h.x);
#elif defined(SYCL_DEVICE_ONLY)
  return numext::bit_cast<numext::uint16_t>(h);
#else
  return h.x;
#endif
}

union float32_bits {
  unsigned int u;
  float f;
};

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC __half_raw float_to_half_rtne(float ff) {
#if (defined(EIGEN_HAS_CUDA_FP16) && defined(EIGEN_CUDA_ARCH) && EIGEN_CUDA_ARCH >= 300) || \
  (defined(EIGEN_HAS_HIP_FP16) && defined(EIGEN_HIP_DEVICE_COMPILE))
  __half tmp_ff = __float2half(ff);
  return *(__half_raw*)&tmp_ff;

#elif defined(EIGEN_HAS_FP16_C)
  __half_raw h;
  h.x = _cvtss_sh(ff, 0);
  return h;

#elif defined(EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC)
  __half_raw h;
  h.x = static_cast<__fp16>(ff);
  return h;

#else
  float32_bits f; f.f = ff;

  const float32_bits f32infty = { 255 << 23 };
  const float32_bits f16max = { (127 + 16) << 23 };
  const float32_bits denorm_magic = { ((127 - 15) + (23 - 10) + 1) << 23 };
  unsigned int sign_mask = 0x80000000u;
  __half_raw o;
  o.x = static_cast<numext::uint16_t>(0x0u);

  unsigned int sign = f.u & sign_mask;
  f.u ^= sign;

  if (f.u >= f16max.u) {  
    o.x = (f.u > f32infty.u) ? 0x7e00 : 0x7c00; 
  } else {  
    if (f.u < (113 << 23)) {  

      f.f += denorm_magic.f;

      o.x = static_cast<numext::uint16_t>(f.u - denorm_magic.u);
    } else {
      unsigned int mant_odd = (f.u >> 13) & 1; 

      f.u += 0xc8000fffU;
      
      f.u += mant_odd;
      
      o.x = static_cast<numext::uint16_t>(f.u >> 13);
    }
  }

  o.x |= static_cast<numext::uint16_t>(sign >> 16);
  return o;
#endif
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC float half_to_float(__half_raw h) {
#if (defined(EIGEN_HAS_CUDA_FP16) && defined(EIGEN_CUDA_ARCH) && EIGEN_CUDA_ARCH >= 300) || \
  (defined(EIGEN_HAS_HIP_FP16) && defined(EIGEN_HIP_DEVICE_COMPILE))
  return __half2float(h);
#elif defined(EIGEN_HAS_FP16_C)
  return _cvtsh_ss(h.x);
#elif defined(EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC)
  return static_cast<float>(h.x);
#else
  const float32_bits magic = { 113 << 23 };
  const unsigned int shifted_exp = 0x7c00 << 13; 
  float32_bits o;

  o.u = (h.x & 0x7fff) << 13;             
  unsigned int exp = shifted_exp & o.u;   
  o.u += (127 - 15) << 23;                

  if (exp == shifted_exp) {     
    o.u += (128 - 16) << 23;    
  } else if (exp == 0) {        
    o.u += 1 << 23;             
    o.f -= magic.f;             
  }

  o.u |= (h.x & 0x8000) << 16;    
  return o.f;
#endif
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool (isinf)(const half& a) {
#ifdef EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC
  return (numext::bit_cast<numext::uint16_t>(a.x) & 0x7fff) == 0x7c00;
#else
  return (a.x & 0x7fff) == 0x7c00;
#endif
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool (isnan)(const half& a) {
#if (defined(EIGEN_HAS_CUDA_FP16) && defined(EIGEN_CUDA_ARCH) && EIGEN_CUDA_ARCH >= 530) || \
  (defined(EIGEN_HAS_HIP_FP16) && defined(EIGEN_HIP_DEVICE_COMPILE))
  return __hisnan(a);
#elif defined(EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC)
  return (numext::bit_cast<numext::uint16_t>(a.x) & 0x7fff) > 0x7c00;
#else
  return (a.x & 0x7fff) > 0x7c00;
#endif
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool (isfinite)(const half& a) {
  return !(isinf EIGEN_NOT_A_MACRO (a)) && !(isnan EIGEN_NOT_A_MACRO (a));
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half abs(const half& a) {
#if defined(EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC)
  return half(vabsh_f16(a.x));
#else
  half result;
  result.x = a.x & 0x7FFF;
  return result;
#endif
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half exp(const half& a) {
#if (EIGEN_CUDA_SDK_VER >= 80000 && defined EIGEN_CUDA_ARCH && EIGEN_CUDA_ARCH >= 530) || \
  defined(EIGEN_HIP_DEVICE_COMPILE)
  return half(hexp(a));
#else
   return half(::expf(float(a)));
#endif
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half expm1(const half& a) {
  return half(numext::expm1(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half log(const half& a) {
#if (defined(EIGEN_HAS_CUDA_FP16) && EIGEN_CUDA_SDK_VER >= 80000 && defined(EIGEN_CUDA_ARCH) && EIGEN_CUDA_ARCH >= 530) || \
  (defined(EIGEN_HAS_HIP_FP16) && defined(EIGEN_HIP_DEVICE_COMPILE))
  return half(::hlog(a));
#else
  return half(::logf(float(a)));
#endif
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half log1p(const half& a) {
  return half(numext::log1p(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half log10(const half& a) {
  return half(::log10f(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half log2(const half& a) {
  return half(static_cast<float>(EIGEN_LOG2E) * ::logf(float(a)));
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half sqrt(const half& a) {
#if (EIGEN_CUDA_SDK_VER >= 80000 && defined EIGEN_CUDA_ARCH && EIGEN_CUDA_ARCH >= 530) || \
  defined(EIGEN_HIP_DEVICE_COMPILE)
  return half(hsqrt(a));
#else
    return half(::sqrtf(float(a)));
#endif
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half pow(const half& a, const half& b) {
  return half(::powf(float(a), float(b)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half sin(const half& a) {
  return half(::sinf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half cos(const half& a) {
  return half(::cosf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half tan(const half& a) {
  return half(::tanf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half tanh(const half& a) {
  return half(::tanhf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half asin(const half& a) {
  return half(::asinf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half acos(const half& a) {
  return half(::acosf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half floor(const half& a) {
#if (EIGEN_CUDA_SDK_VER >= 80000 && defined EIGEN_CUDA_ARCH && EIGEN_CUDA_ARCH >= 300) || \
  defined(EIGEN_HIP_DEVICE_COMPILE)
  return half(hfloor(a));
#else
  return half(::floorf(float(a)));
#endif
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half ceil(const half& a) {
#if (EIGEN_CUDA_SDK_VER >= 80000 && defined EIGEN_CUDA_ARCH && EIGEN_CUDA_ARCH >= 300) || \
  defined(EIGEN_HIP_DEVICE_COMPILE)
  return half(hceil(a));
#else
  return half(::ceilf(float(a)));
#endif
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half rint(const half& a) {
  return half(::rintf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half round(const half& a) {
  return half(::roundf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half fmod(const half& a, const half& b) {
  return half(::fmodf(float(a), float(b)));
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half (min)(const half& a, const half& b) {
#if (defined(EIGEN_HAS_CUDA_FP16) && defined(EIGEN_CUDA_ARCH) && EIGEN_CUDA_ARCH >= 530) || \
  (defined(EIGEN_HAS_HIP_FP16) && defined(EIGEN_HIP_DEVICE_COMPILE))
  return __hlt(b, a) ? b : a;
#else
  const float f1 = static_cast<float>(a);
  const float f2 = static_cast<float>(b);
  return f2 < f1 ? b : a;
#endif
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC half (max)(const half& a, const half& b) {
#if (defined(EIGEN_HAS_CUDA_FP16) && defined(EIGEN_CUDA_ARCH) && EIGEN_CUDA_ARCH >= 530) || \
  (defined(EIGEN_HAS_HIP_FP16) && defined(EIGEN_HIP_DEVICE_COMPILE))
  return __hlt(a, b) ? b : a;
#else
  const float f1 = static_cast<float>(a);
  const float f2 = static_cast<float>(b);
  return f1 < f2 ? b : a;
#endif
}

#ifndef EIGEN_NO_IO
EIGEN_ALWAYS_INLINE std::ostream& operator << (std::ostream& os, const half& v) {
  os << static_cast<float>(v);
  return os;
}
#endif

} 

namespace internal {

template<>
struct random_default_impl<half, false, false>
{
  static inline half run(const half& x, const half& y)
  {
    return x + (y-x) * half(float(std::rand()) / float(RAND_MAX));
  }
  static inline half run()
  {
    return run(half(-1.f), half(1.f));
  }
};

template<> struct is_arithmetic<half> { enum { value = true }; };

} 

template<> struct NumTraits<Eigen::half>
    : GenericNumTraits<Eigen::half>
{
  enum {
    IsSigned = true,
    IsInteger = false,
    IsComplex = false,
    RequireInitialization = false
  };

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::half epsilon() {
    return half_impl::raw_uint16_to_half(0x0800);
  }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::half dummy_precision() {
    return half_impl::raw_uint16_to_half(0x211f); 
  }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::half highest() {
    return half_impl::raw_uint16_to_half(0x7bff);
  }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::half lowest() {
    return half_impl::raw_uint16_to_half(0xfbff);
  }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::half infinity() {
    return half_impl::raw_uint16_to_half(0x7c00);
  }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::half quiet_NaN() {
    return half_impl::raw_uint16_to_half(0x7e00);
  }
};

} 

#if defined(EIGEN_HAS_GPU_FP16) || defined(EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC)
  #pragma pop_macro("EIGEN_CONSTEXPR")
#endif

namespace Eigen {
namespace numext {

#if defined(EIGEN_GPU_COMPILE_PHASE)

template <>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE bool(isnan)(const Eigen::half& h) {
  return (half_impl::isnan)(h);
}

template <>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE bool(isinf)(const Eigen::half& h) {
  return (half_impl::isinf)(h);
}

template <>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE bool(isfinite)(const Eigen::half& h) {
  return (half_impl::isfinite)(h);
}

#endif

template <>
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC Eigen::half bit_cast<Eigen::half, uint16_t>(const uint16_t& src) {
  return Eigen::half(Eigen::half_impl::raw_uint16_to_half(src));
}

template <>
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC uint16_t bit_cast<uint16_t, Eigen::half>(const Eigen::half& src) {
  return Eigen::half_impl::raw_half_as_uint16(src);
}

}  
}  

#if (defined(EIGEN_CUDACC) && (!defined(EIGEN_CUDA_ARCH) || EIGEN_CUDA_ARCH >= 300)) \
    || defined(EIGEN_HIPCC)

#if defined(EIGEN_HAS_CUDA_FP16) && EIGEN_CUDA_SDK_VER >= 90000

__device__ EIGEN_STRONG_INLINE Eigen::half __shfl_sync(unsigned mask, Eigen::half var, int srcLane, int width=warpSize) {
  const __half h = var;
  return static_cast<Eigen::half>(__shfl_sync(mask, h, srcLane, width));
}

__device__ EIGEN_STRONG_INLINE Eigen::half __shfl_up_sync(unsigned mask, Eigen::half var, unsigned int delta, int width=warpSize) {
  const __half h = var;
  return static_cast<Eigen::half>(__shfl_up_sync(mask, h, delta, width));
}

__device__ EIGEN_STRONG_INLINE Eigen::half __shfl_down_sync(unsigned mask, Eigen::half var, unsigned int delta, int width=warpSize) {
  const __half h = var;
  return static_cast<Eigen::half>(__shfl_down_sync(mask, h, delta, width));
}

__device__ EIGEN_STRONG_INLINE Eigen::half __shfl_xor_sync(unsigned mask, Eigen::half var, int laneMask, int width=warpSize) {
  const __half h = var;
  return static_cast<Eigen::half>(__shfl_xor_sync(mask, h, laneMask, width));
}

#else 

__device__ EIGEN_STRONG_INLINE Eigen::half __shfl(Eigen::half var, int srcLane, int width=warpSize) {
  const int ivar = static_cast<int>(Eigen::numext::bit_cast<Eigen::numext::uint16_t>(var));
  return Eigen::numext::bit_cast<Eigen::half>(static_cast<Eigen::numext::uint16_t>(__shfl(ivar, srcLane, width)));
}

__device__ EIGEN_STRONG_INLINE Eigen::half __shfl_up(Eigen::half var, unsigned int delta, int width=warpSize) {
  const int ivar = static_cast<int>(Eigen::numext::bit_cast<Eigen::numext::uint16_t>(var));
  return Eigen::numext::bit_cast<Eigen::half>(static_cast<Eigen::numext::uint16_t>(__shfl_up(ivar, delta, width)));
}

__device__ EIGEN_STRONG_INLINE Eigen::half __shfl_down(Eigen::half var, unsigned int delta, int width=warpSize) {
  const int ivar = static_cast<int>(Eigen::numext::bit_cast<Eigen::numext::uint16_t>(var));
  return Eigen::numext::bit_cast<Eigen::half>(static_cast<Eigen::numext::uint16_t>(__shfl_down(ivar, delta, width)));
}

__device__ EIGEN_STRONG_INLINE Eigen::half __shfl_xor(Eigen::half var, int laneMask, int width=warpSize) {
  const int ivar = static_cast<int>(Eigen::numext::bit_cast<Eigen::numext::uint16_t>(var));
  return Eigen::numext::bit_cast<Eigen::half>(static_cast<Eigen::numext::uint16_t>(__shfl_xor(ivar, laneMask, width)));
}

#endif 
#endif 

#if (defined(EIGEN_CUDACC) && (!defined(EIGEN_CUDA_ARCH) || EIGEN_CUDA_ARCH >= 350)) \
    || defined(EIGEN_HIPCC)
EIGEN_STRONG_INLINE __device__ Eigen::half __ldg(const Eigen::half* ptr) {
  return Eigen::half_impl::raw_uint16_to_half(__ldg(reinterpret_cast<const Eigen::numext::uint16_t*>(ptr)));
}
#endif 

#if EIGEN_HAS_STD_HASH
namespace std {
template <>
struct hash<Eigen::half> {
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE std::size_t operator()(const Eigen::half& a) const {
    return static_cast<std::size_t>(Eigen::numext::bit_cast<Eigen::numext::uint16_t>(a));
  }
};
} 
#endif

#endif 
