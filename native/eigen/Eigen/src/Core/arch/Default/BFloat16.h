#ifndef EIGEN_BFLOAT16_H
#define EIGEN_BFLOAT16_H

#define BF16_PACKET_FUNCTION(PACKET_F, PACKET_BF16, METHOD)         \
  template <>                                                       \
  EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS EIGEN_UNUSED  \
  PACKET_BF16 METHOD<PACKET_BF16>(const PACKET_BF16& _x) {          \
    return F32ToBf16(METHOD<PACKET_F>(Bf16ToF32(_x)));              \
  }

namespace Eigen {

struct bfloat16;

namespace bfloat16_impl {

struct __bfloat16_raw {
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR __bfloat16_raw() : value(0) {}
  explicit EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR __bfloat16_raw(unsigned short raw) : value(raw) {}
  unsigned short value;
};

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR __bfloat16_raw raw_uint16_to_bfloat16(unsigned short value);
template <bool AssumeArgumentIsNormalOrInfinityOrZero>
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC __bfloat16_raw float_to_bfloat16_rtne(float ff);

template <>
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC __bfloat16_raw float_to_bfloat16_rtne<false>(float ff);
template <>
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC __bfloat16_raw float_to_bfloat16_rtne<true>(float ff);
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC float bfloat16_to_float(__bfloat16_raw h);

struct bfloat16_base : public __bfloat16_raw {
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR bfloat16_base() {}
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR bfloat16_base(const __bfloat16_raw& h) : __bfloat16_raw(h) {}
};

} 

struct bfloat16 : public bfloat16_impl::bfloat16_base {

  typedef bfloat16_impl::__bfloat16_raw __bfloat16_raw;

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR bfloat16() {}

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR bfloat16(const __bfloat16_raw& h) : bfloat16_impl::bfloat16_base(h) {}

  explicit EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR bfloat16(bool b)
      : bfloat16_impl::bfloat16_base(bfloat16_impl::raw_uint16_to_bfloat16(b ? 0x3f80 : 0)) {}

  template<class T>
  explicit EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR bfloat16(T val)
      : bfloat16_impl::bfloat16_base(bfloat16_impl::float_to_bfloat16_rtne<internal::is_integral<T>::value>(static_cast<float>(val))) {}

  explicit EIGEN_DEVICE_FUNC bfloat16(float f)
      : bfloat16_impl::bfloat16_base(bfloat16_impl::float_to_bfloat16_rtne<false>(f)) {}

  template<typename RealScalar>
  explicit EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR bfloat16(const std::complex<RealScalar>& val)
      : bfloat16_impl::bfloat16_base(bfloat16_impl::float_to_bfloat16_rtne<false>(static_cast<float>(val.real()))) {}

  EIGEN_DEVICE_FUNC operator float() const {  
    return bfloat16_impl::bfloat16_to_float(*this);
  }
};
} 

namespace std {
template<>
struct numeric_limits<Eigen::bfloat16> {
  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool has_infinity = true;
  static const bool has_quiet_NaN = true;
  static const bool has_signaling_NaN = true;
  static const float_denorm_style has_denorm = std::denorm_absent;
  static const bool has_denorm_loss = false;
  static const std::float_round_style round_style = numeric_limits<float>::round_style;
  static const bool is_iec559 = false;
  static const bool is_bounded = true;
  static const bool is_modulo = false;
  static const int digits = 8;
  static const int digits10 = 2;
  static const int max_digits10 = 4;
  static const int radix = 2;
  static const int min_exponent = numeric_limits<float>::min_exponent;
  static const int min_exponent10 = numeric_limits<float>::min_exponent10;
  static const int max_exponent = numeric_limits<float>::max_exponent;
  static const int max_exponent10 = numeric_limits<float>::max_exponent10;
  static const bool traps = numeric_limits<float>::traps;
  static const bool tinyness_before = numeric_limits<float>::tinyness_before;

  static Eigen::bfloat16 (min)() { return Eigen::bfloat16_impl::raw_uint16_to_bfloat16(0x0080); }
  static Eigen::bfloat16 lowest() { return Eigen::bfloat16_impl::raw_uint16_to_bfloat16(0xff7f); }
  static Eigen::bfloat16 (max)() { return Eigen::bfloat16_impl::raw_uint16_to_bfloat16(0x7f7f); }
  static Eigen::bfloat16 epsilon() { return Eigen::bfloat16_impl::raw_uint16_to_bfloat16(0x3c00); }
  static Eigen::bfloat16 round_error() { return Eigen::bfloat16(0x3f00); }
  static Eigen::bfloat16 infinity() { return Eigen::bfloat16_impl::raw_uint16_to_bfloat16(0x7f80); }
  static Eigen::bfloat16 quiet_NaN() { return Eigen::bfloat16_impl::raw_uint16_to_bfloat16(0x7fc0); }
  static Eigen::bfloat16 signaling_NaN() { return Eigen::bfloat16_impl::raw_uint16_to_bfloat16(0x7f81); }
  static Eigen::bfloat16 denorm_min() { return Eigen::bfloat16_impl::raw_uint16_to_bfloat16(0x0001); }
};

template<>
struct numeric_limits<const Eigen::bfloat16> : numeric_limits<Eigen::bfloat16> {};
template<>
struct numeric_limits<volatile Eigen::bfloat16> : numeric_limits<Eigen::bfloat16> {};
template<>
struct numeric_limits<const volatile Eigen::bfloat16> : numeric_limits<Eigen::bfloat16> {};
} 

namespace Eigen {

namespace bfloat16_impl {

#if !defined(EIGEN_HAS_NATIVE_BF16) || (EIGEN_COMP_CLANG && !EIGEN_COMP_NVCC) 

#if EIGEN_COMP_CLANG && defined(EIGEN_CUDACC)

#pragma push_macro("EIGEN_DEVICE_FUNC")
#undef EIGEN_DEVICE_FUNC
#if defined(EIGEN_HAS_CUDA_BF16) && defined(EIGEN_HAS_NATIVE_BF16)
#define EIGEN_DEVICE_FUNC __host__
#else 
#define EIGEN_DEVICE_FUNC __host__ __device__
#endif
#endif

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator + (const bfloat16& a, const bfloat16& b) {
  return bfloat16(float(a) + float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator + (const bfloat16& a, const int& b) {
  return bfloat16(float(a) + static_cast<float>(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator + (const int& a, const bfloat16& b) {
  return bfloat16(static_cast<float>(a) + float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator * (const bfloat16& a, const bfloat16& b) {
  return bfloat16(float(a) * float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator - (const bfloat16& a, const bfloat16& b) {
  return bfloat16(float(a) - float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator / (const bfloat16& a, const bfloat16& b) {
  return bfloat16(float(a) / float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator - (const bfloat16& a) {
  bfloat16 result;
  result.value = a.value ^ 0x8000;
  return result;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16& operator += (bfloat16& a, const bfloat16& b) {
  a = bfloat16(float(a) + float(b));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16& operator *= (bfloat16& a, const bfloat16& b) {
  a = bfloat16(float(a) * float(b));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16& operator -= (bfloat16& a, const bfloat16& b) {
  a = bfloat16(float(a) - float(b));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16& operator /= (bfloat16& a, const bfloat16& b) {
  a = bfloat16(float(a) / float(b));
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator++(bfloat16& a) {
  a += bfloat16(1);
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator--(bfloat16& a) {
  a -= bfloat16(1);
  return a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator++(bfloat16& a, int) {
  bfloat16 original_value = a;
  ++a;
  return original_value;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator--(bfloat16& a, int) {
  bfloat16 original_value = a;
  --a;
  return original_value;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator == (const bfloat16& a, const bfloat16& b) {
  return numext::equal_strict(float(a),float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator != (const bfloat16& a, const bfloat16& b) {
  return numext::not_equal_strict(float(a), float(b));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator < (const bfloat16& a, const bfloat16& b) {
  return float(a) < float(b);
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator <= (const bfloat16& a, const bfloat16& b) {
  return float(a) <= float(b);
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator > (const bfloat16& a, const bfloat16& b) {
  return float(a) > float(b);
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool operator >= (const bfloat16& a, const bfloat16& b) {
  return float(a) >= float(b);
}

#if EIGEN_COMP_CLANG && defined(EIGEN_CUDACC)
#pragma pop_macro("EIGEN_DEVICE_FUNC")
#endif
#endif  

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 operator / (const bfloat16& a, Index b) {
  return bfloat16(static_cast<float>(a) / static_cast<float>(b));
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC __bfloat16_raw truncate_to_bfloat16(const float v) {
  __bfloat16_raw output;
  if (Eigen::numext::isnan EIGEN_NOT_A_MACRO(v)) {
    output.value = std::signbit(v) ? 0xFFC0: 0x7FC0;
    return output;
  }
  const uint16_t* p = reinterpret_cast<const uint16_t*>(&v);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  output.value = p[0];
#else
  output.value = p[1];
#endif
  return output;
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR __bfloat16_raw raw_uint16_to_bfloat16(numext::uint16_t value) {
  return __bfloat16_raw(value);
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR numext::uint16_t raw_bfloat16_as_uint16(const __bfloat16_raw& bf) {
  return bf.value;
}

template <>
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC __bfloat16_raw float_to_bfloat16_rtne<false>(float ff) {
#if (defined(EIGEN_HAS_CUDA_BF16) && defined(EIGEN_HAS_HIP_BF16))
  
#else
  __bfloat16_raw output;

  if (Eigen::numext::isnan EIGEN_NOT_A_MACRO(ff)) {

    output.value = std::signbit(ff) ? 0xFFC0: 0x7FC0;
  } else {

    output = float_to_bfloat16_rtne<true>(ff);
  }
  return output;
#endif
}

template <>
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC __bfloat16_raw float_to_bfloat16_rtne<true>(float ff) {
#if (defined(EIGEN_HAS_CUDA_BF16) && defined(EIGEN_HAS_HIP_BF16))
    
#else
    numext::uint32_t input = numext::bit_cast<numext::uint32_t>(ff);
    __bfloat16_raw output;

    numext::uint32_t lsb = (input >> 16) & 1;
    numext::uint32_t rounding_bias = 0x7fff + lsb;
    input += rounding_bias;
    output.value = static_cast<numext::uint16_t>(input >> 16);
    return output;
#endif
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC float bfloat16_to_float(__bfloat16_raw h) {
    float result = 0;
    unsigned short* q = reinterpret_cast<unsigned short*>(&result);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    q[0] = h.value;
#else
    q[1] = h.value;
#endif
    return result;
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool (isinf)(const bfloat16& a) {
  EIGEN_USING_STD(isinf);
  return (isinf)(float(a));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool (isnan)(const bfloat16& a) {
  EIGEN_USING_STD(isnan);
  return (isnan)(float(a));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bool (isfinite)(const bfloat16& a) {
  return !(isinf EIGEN_NOT_A_MACRO (a)) && !(isnan EIGEN_NOT_A_MACRO (a));
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 abs(const bfloat16& a) {
  bfloat16 result;
  result.value = a.value & 0x7FFF;
  return result;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 exp(const bfloat16& a) {
   return bfloat16(::expf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 expm1(const bfloat16& a) {
  return bfloat16(numext::expm1(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 log(const bfloat16& a) {
  return bfloat16(::logf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 log1p(const bfloat16& a) {
  return bfloat16(numext::log1p(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 log10(const bfloat16& a) {
  return bfloat16(::log10f(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 log2(const bfloat16& a) {
  return bfloat16(static_cast<float>(EIGEN_LOG2E) * ::logf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 sqrt(const bfloat16& a) {
    return bfloat16(::sqrtf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 pow(const bfloat16& a, const bfloat16& b) {
  return bfloat16(::powf(float(a), float(b)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 sin(const bfloat16& a) {
  return bfloat16(::sinf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 cos(const bfloat16& a) {
  return bfloat16(::cosf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 tan(const bfloat16& a) {
  return bfloat16(::tanf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 asin(const bfloat16& a) {
  return bfloat16(::asinf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 acos(const bfloat16& a) {
  return bfloat16(::acosf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 atan(const bfloat16& a) {
  return bfloat16(::atanf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 sinh(const bfloat16& a) {
  return bfloat16(::sinhf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 cosh(const bfloat16& a) {
  return bfloat16(::coshf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 tanh(const bfloat16& a) {
  return bfloat16(::tanhf(float(a)));
}
#if EIGEN_HAS_CXX11_MATH
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 asinh(const bfloat16& a) {
  return bfloat16(::asinhf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 acosh(const bfloat16& a) {
  return bfloat16(::acoshf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 atanh(const bfloat16& a) {
  return bfloat16(::atanhf(float(a)));
}
#endif
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 floor(const bfloat16& a) {
  return bfloat16(::floorf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 ceil(const bfloat16& a) {
  return bfloat16(::ceilf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 rint(const bfloat16& a) {
  return bfloat16(::rintf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 round(const bfloat16& a) {
  return bfloat16(::roundf(float(a)));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 fmod(const bfloat16& a, const bfloat16& b) {
  return bfloat16(::fmodf(float(a), float(b)));
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 (min)(const bfloat16& a, const bfloat16& b) {
  const float f1 = static_cast<float>(a);
  const float f2 = static_cast<float>(b);
  return f2 < f1 ? b : a;
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 (max)(const bfloat16& a, const bfloat16& b) {
  const float f1 = static_cast<float>(a);
  const float f2 = static_cast<float>(b);
  return f1 < f2 ? b : a;
}

EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 fmin(const bfloat16& a, const bfloat16& b) {
  const float f1 = static_cast<float>(a);
  const float f2 = static_cast<float>(b);
  return bfloat16(::fminf(f1, f2));
}
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC bfloat16 fmax(const bfloat16& a, const bfloat16& b) {
  const float f1 = static_cast<float>(a);
  const float f2 = static_cast<float>(b);
  return bfloat16(::fmaxf(f1, f2));
}

#ifndef EIGEN_NO_IO
EIGEN_ALWAYS_INLINE std::ostream& operator << (std::ostream& os, const bfloat16& v) {
  os << static_cast<float>(v);
  return os;
}
#endif

} 

namespace internal {

template<>
struct random_default_impl<bfloat16, false, false>
{
  static inline bfloat16 run(const bfloat16& x, const bfloat16& y)
  {
    return x + (y-x) * bfloat16(float(std::rand()) / float(RAND_MAX));
  }
  static inline bfloat16 run()
  {
    return run(bfloat16(-1.f), bfloat16(1.f));
  }
};

template<> struct is_arithmetic<bfloat16> { enum { value = true }; };

} 

template<> struct NumTraits<Eigen::bfloat16>
    : GenericNumTraits<Eigen::bfloat16>
{
  enum {
    IsSigned = true,
    IsInteger = false,
    IsComplex = false,
    RequireInitialization = false
  };

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::bfloat16 epsilon() {
    return bfloat16_impl::raw_uint16_to_bfloat16(0x3c00);
  }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::bfloat16 dummy_precision() {
    return bfloat16_impl::raw_uint16_to_bfloat16(0x3D4D);  

  }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::bfloat16 highest() {
    return bfloat16_impl::raw_uint16_to_bfloat16(0x7F7F);
  }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::bfloat16 lowest() {
    return bfloat16_impl::raw_uint16_to_bfloat16(0xFF7F);
  }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::bfloat16 infinity() {
    return bfloat16_impl::raw_uint16_to_bfloat16(0x7f80);
  }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static EIGEN_STRONG_INLINE Eigen::bfloat16 quiet_NaN() {
    return bfloat16_impl::raw_uint16_to_bfloat16(0x7fc0);
  }
};

} 

namespace Eigen {
namespace numext {

template<>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE
bool (isnan)(const Eigen::bfloat16& h) {
  return (bfloat16_impl::isnan)(h);
}

template<>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE
bool (isinf)(const Eigen::bfloat16& h) {
  return (bfloat16_impl::isinf)(h);
}

template<>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE
bool (isfinite)(const Eigen::bfloat16& h) {
  return (bfloat16_impl::isfinite)(h);
}

template <>
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC Eigen::bfloat16 bit_cast<Eigen::bfloat16, uint16_t>(const uint16_t& src) {
  return Eigen::bfloat16(Eigen::bfloat16_impl::raw_uint16_to_bfloat16(src));
}

template <>
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC uint16_t bit_cast<uint16_t, Eigen::bfloat16>(const Eigen::bfloat16& src) {
  return Eigen::bfloat16_impl::raw_bfloat16_as_uint16(src);
}

}  
}  

#if EIGEN_HAS_STD_HASH
namespace std {
template <>
struct hash<Eigen::bfloat16> {
  EIGEN_STRONG_INLINE std::size_t operator()(const Eigen::bfloat16& a) const {
    return static_cast<std::size_t>(Eigen::numext::bit_cast<Eigen::numext::uint16_t>(a));
  }
};
} 
#endif

#endif 
