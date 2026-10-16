#ifndef EIGEN_GENERIC_PACKET_MATH_H
#define EIGEN_GENERIC_PACKET_MATH_H

namespace Eigen {

namespace internal {

#ifndef EIGEN_DEBUG_ALIGNED_LOAD
#define EIGEN_DEBUG_ALIGNED_LOAD
#endif

#ifndef EIGEN_DEBUG_UNALIGNED_LOAD
#define EIGEN_DEBUG_UNALIGNED_LOAD
#endif

#ifndef EIGEN_DEBUG_ALIGNED_STORE
#define EIGEN_DEBUG_ALIGNED_STORE
#endif

#ifndef EIGEN_DEBUG_UNALIGNED_STORE
#define EIGEN_DEBUG_UNALIGNED_STORE
#endif

struct default_packet_traits
{
  enum {
    HasHalfPacket = 0,

    HasAdd       = 1,
    HasSub       = 1,
    HasShift     = 1,
    HasMul       = 1,
    HasNegate    = 1,
    HasAbs       = 1,
    HasArg       = 0,
    HasAbs2      = 1,
    HasAbsDiff   = 0,
    HasMin       = 1,
    HasMax       = 1,
    HasConj      = 1,
    HasSetLinear = 1,
    HasBlend     = 0,

    HasCmp       = 0,

    HasDiv    = 0,
    HasSqrt   = 0,
    HasRsqrt  = 0,
    HasExp    = 0,
    HasExpm1  = 0,
    HasLog    = 0,
    HasLog1p  = 0,
    HasLog10  = 0,
    HasPow    = 0,

    HasSin    = 0,
    HasCos    = 0,
    HasTan    = 0,
    HasASin   = 0,
    HasACos   = 0,
    HasATan   = 0,
    HasSinh   = 0,
    HasCosh   = 0,
    HasTanh   = 0,
    HasLGamma = 0,
    HasDiGamma = 0,
    HasZeta = 0,
    HasPolygamma = 0,
    HasErf = 0,
    HasErfc = 0,
    HasNdtri = 0,
    HasBessel = 0,
    HasIGamma = 0,
    HasIGammaDerA = 0,
    HasGammaSampleDerAlpha = 0,
    HasIGammac = 0,
    HasBetaInc = 0,

    HasRound  = 0,
    HasRint   = 0,
    HasFloor  = 0,
    HasCeil   = 0,
    HasSign   = 0
  };
};

template<typename T> struct packet_traits : default_packet_traits
{
  typedef T type;
  typedef T half;
  enum {
    Vectorizable = 0,
    size = 1,
    AlignedOnScalar = 0,
    HasHalfPacket = 0
  };
  enum {
    HasAdd    = 0,
    HasSub    = 0,
    HasMul    = 0,
    HasNegate = 0,
    HasAbs    = 0,
    HasAbs2   = 0,
    HasMin    = 0,
    HasMax    = 0,
    HasConj   = 0,
    HasSetLinear = 0
  };
};

template<typename T> struct packet_traits<const T> : packet_traits<T> { };

template<typename T> struct unpacket_traits
{
  typedef T type;
  typedef T half;
  enum
  {
    size = 1,
    alignment = 1,
    vectorizable = false,
    masked_load_available=false,
    masked_store_available=false
  };
};

template<typename T> struct unpacket_traits<const T> : unpacket_traits<T> { };

template <typename Src, typename Tgt> struct type_casting_traits {
  enum {
    VectorizedCast = 0,
    SrcCoeffRatio = 1,
    TgtCoeffRatio = 1
  };
};

template<typename T, int unique_id = 0>
struct eigen_packet_wrapper
{
  EIGEN_ALWAYS_INLINE operator T&() { return m_val; }
  EIGEN_ALWAYS_INLINE operator const T&() const { return m_val; }
  EIGEN_ALWAYS_INLINE eigen_packet_wrapper() {}
  EIGEN_ALWAYS_INLINE eigen_packet_wrapper(const T &v) : m_val(v) {}
  EIGEN_ALWAYS_INLINE eigen_packet_wrapper& operator=(const T &v) {
    m_val = v;
    return *this;
  }

  T m_val;
};

template<typename Packet>
struct is_scalar {
  typedef typename unpacket_traits<Packet>::type Scalar;
  enum {
    value = internal::is_same<Packet, Scalar>::value
  };
};

template <typename SrcPacket, typename TgtPacket>
EIGEN_DEVICE_FUNC inline TgtPacket
pcast(const SrcPacket& a) {
  return static_cast<TgtPacket>(a);
}
template <typename SrcPacket, typename TgtPacket>
EIGEN_DEVICE_FUNC inline TgtPacket
pcast(const SrcPacket& a, const SrcPacket& ) {
  return static_cast<TgtPacket>(a);
}
template <typename SrcPacket, typename TgtPacket>
EIGEN_DEVICE_FUNC inline TgtPacket
pcast(const SrcPacket& a, const SrcPacket& , const SrcPacket& , const SrcPacket& ) {
  return static_cast<TgtPacket>(a);
}
template <typename SrcPacket, typename TgtPacket>
EIGEN_DEVICE_FUNC inline TgtPacket
pcast(const SrcPacket& a, const SrcPacket& , const SrcPacket& , const SrcPacket& ,
      const SrcPacket& , const SrcPacket& , const SrcPacket& , const SrcPacket& ) {
  return static_cast<TgtPacket>(a);
}

template <typename Target, typename Packet>
EIGEN_DEVICE_FUNC inline Target
preinterpret(const Packet& a); 

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
padd(const Packet& a, const Packet& b) { return a+b; }

template<> EIGEN_DEVICE_FUNC inline bool
padd(const bool& a, const bool& b) { return a || b; }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
psub(const Packet& a, const Packet& b) { return a-b; }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pnegate(const Packet& a) { return -a; }

template<> EIGEN_DEVICE_FUNC inline bool
pnegate(const bool& a) { return !a; }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pconj(const Packet& a) { return numext::conj(a); }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pmul(const Packet& a, const Packet& b) { return a*b; }

template<> EIGEN_DEVICE_FUNC inline bool
pmul(const bool& a, const bool& b) { return a && b; }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pdiv(const Packet& a, const Packet& b) { return a/b; }

template<typename Packet, typename EnableIf = void>
struct ptrue_impl {
  static EIGEN_DEVICE_FUNC inline Packet run(const Packet& ){
    Packet b;
    memset(static_cast<void*>(&b), 0xff, sizeof(Packet));
    return b;
  }
};

template<typename T>
struct ptrue_impl<T, 
    typename internal::enable_if<is_scalar<T>::value && NumTraits<T>::RequireInitialization>::type > {
  static EIGEN_DEVICE_FUNC inline T run(const T& ){
    return T(1);
  }
};

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
ptrue(const Packet& a) {
  return ptrue_impl<Packet>::run(a);
}

template<typename Packet, typename EnableIf = void>
struct pzero_impl {
  static EIGEN_DEVICE_FUNC inline Packet run(const Packet& ) {
    Packet b;
    memset(static_cast<void*>(&b), 0x00, sizeof(Packet));
    return b;
  }
};

template<typename T>
struct pzero_impl<T,
    typename internal::enable_if<is_scalar<T>::value>::type> {
  static EIGEN_DEVICE_FUNC inline T run(const T& ) {
    return T(0);
  }
};

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pzero(const Packet& a) {
  return pzero_impl<Packet>::run(a);
}

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pcmp_le(const Packet& a, const Packet& b)  { return a<=b ? ptrue(a) : pzero(a); }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pcmp_lt(const Packet& a, const Packet& b)  { return a<b ? ptrue(a) : pzero(a); }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pcmp_eq(const Packet& a, const Packet& b) { return a==b ? ptrue(a) : pzero(a); }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pcmp_lt_or_nan(const Packet& a, const Packet& b) { return a>=b ? pzero(a) : ptrue(a); }

template<typename T>
struct bit_and {
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR EIGEN_ALWAYS_INLINE T operator()(const T& a, const T& b) const {
    return a & b;
  }
};

template<typename T>
struct bit_or {
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR EIGEN_ALWAYS_INLINE T operator()(const T& a, const T& b) const {
    return a | b;
  }
};

template<typename T>
struct bit_xor {
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR EIGEN_ALWAYS_INLINE T operator()(const T& a, const T& b) const {
    return a ^ b;
  }
};

template<typename T>
struct bit_not {
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR EIGEN_ALWAYS_INLINE T operator()(const T& a) const {
    return ~a;
  }
};

template<typename T>
struct operator_bitwise_helper {
  EIGEN_DEVICE_FUNC static inline T bitwise_and(const T& a, const T& b) { return bit_and<T>()(a, b); }
  EIGEN_DEVICE_FUNC static inline T bitwise_or(const T& a, const T& b) { return bit_or<T>()(a, b); }
  EIGEN_DEVICE_FUNC static inline T bitwise_xor(const T& a, const T& b) { return bit_xor<T>()(a, b); }
  EIGEN_DEVICE_FUNC static inline T bitwise_not(const T& a) { return bit_not<T>()(a); }
};

template<typename T>
struct bytewise_bitwise_helper {
  EIGEN_DEVICE_FUNC static inline T bitwise_and(const T& a, const T& b) {
    return binary(a, b, bit_and<unsigned char>());
  }
  EIGEN_DEVICE_FUNC static inline T bitwise_or(const T& a, const T& b) { 
    return binary(a, b, bit_or<unsigned char>());
   }
  EIGEN_DEVICE_FUNC static inline T bitwise_xor(const T& a, const T& b) {
    return binary(a, b, bit_xor<unsigned char>());
  }
  EIGEN_DEVICE_FUNC static inline T bitwise_not(const T& a) { 
    return unary(a,bit_not<unsigned char>());
   }
  
 private:
  template<typename Op>
  EIGEN_DEVICE_FUNC static inline T unary(const T& a, Op op) {
    const unsigned char* a_ptr = reinterpret_cast<const unsigned char*>(&a);
    T c;
    unsigned char* c_ptr = reinterpret_cast<unsigned char*>(&c);
    for (size_t i = 0; i < sizeof(T); ++i) {
      *c_ptr++ = op(*a_ptr++);
    }
    return c;
  }

  template<typename Op>
  EIGEN_DEVICE_FUNC static inline T binary(const T& a, const T& b, Op op) {
    const unsigned char* a_ptr = reinterpret_cast<const unsigned char*>(&a);
    const unsigned char* b_ptr = reinterpret_cast<const unsigned char*>(&b);
    T c;
    unsigned char* c_ptr = reinterpret_cast<unsigned char*>(&c);
    for (size_t i = 0; i < sizeof(T); ++i) {
      *c_ptr++ = op(*a_ptr++, *b_ptr++);
    }
    return c;
  }
};

template<typename T, typename EnableIf = void>
struct bitwise_helper : public bytewise_bitwise_helper<T> {};

template<typename T>
struct bitwise_helper<T,
  typename internal::enable_if<
    is_scalar<T>::value && (NumTraits<T>::IsInteger || NumTraits<T>::RequireInitialization)>::type
  > : public operator_bitwise_helper<T> {};

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pand(const Packet& a, const Packet& b) {
  return bitwise_helper<Packet>::bitwise_and(a, b);
}

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
por(const Packet& a, const Packet& b) {
  return bitwise_helper<Packet>::bitwise_or(a, b);
}

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pxor(const Packet& a, const Packet& b) {
  return bitwise_helper<Packet>::bitwise_xor(a, b);
}

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pnot(const Packet& a) {
  return bitwise_helper<Packet>::bitwise_not(a);
}

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pandnot(const Packet& a, const Packet& b) { return pand(a, pnot(b)); }

template<typename Packet, typename EnableIf = void>
struct pselect_impl {
  static EIGEN_DEVICE_FUNC inline Packet run(const Packet& mask, const Packet& a, const Packet& b) {
    return por(pand(a,mask),pandnot(b,mask));
  }
};

template<typename Packet>
struct pselect_impl<Packet, 
    typename internal::enable_if<is_scalar<Packet>::value>::type > {
  static EIGEN_DEVICE_FUNC inline Packet run(const Packet& mask, const Packet& a, const Packet& b) {
    return numext::equal_strict(mask, Packet(0)) ? b : a;
  }
};

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pselect(const Packet& mask, const Packet& a, const Packet& b) {
  return pselect_impl<Packet>::run(mask, a, b);
}

template<> EIGEN_DEVICE_FUNC inline bool pselect<bool>(
    const bool& cond, const bool& a, const bool& b) {
  return cond ? a : b;
}

template<int NaNPropagation>
struct pminmax_impl {
  template <typename Packet, typename Op>
  static EIGEN_DEVICE_FUNC inline Packet run(const Packet& a, const Packet& b, Op op) {
    return op(a,b);
  }
};

template<>
struct pminmax_impl<PropagateNaN> {
  template <typename Packet, typename Op>
  static EIGEN_DEVICE_FUNC inline Packet run(const Packet& a, const Packet& b, Op op) {
  Packet not_nan_mask_a = pcmp_eq(a, a);
  Packet not_nan_mask_b = pcmp_eq(b, b);
  return pselect(not_nan_mask_a,
                 pselect(not_nan_mask_b, op(a, b), b),
                 a);
  }
};

template<>
struct pminmax_impl<PropagateNumbers> {
  template <typename Packet, typename Op>
  static EIGEN_DEVICE_FUNC inline Packet run(const Packet& a, const Packet& b, Op op) {
  Packet not_nan_mask_a = pcmp_eq(a, a);
  Packet not_nan_mask_b = pcmp_eq(b, b);
  return pselect(not_nan_mask_a,
                 pselect(not_nan_mask_b, op(a, b), a),
                 b);
  }
};

#ifndef SYCL_DEVICE_ONLY
#define EIGEN_BINARY_OP_NAN_PROPAGATION(Type, Func) Func
#else
#define EIGEN_BINARY_OP_NAN_PROPAGATION(Type, Func) \
[](const Type& a, const Type& b) { \
        return Func(a, b);}
#endif

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pmin(const Packet& a, const Packet& b) { return numext::mini(a,b); }

template <int NaNPropagation, typename Packet>
EIGEN_DEVICE_FUNC inline Packet pmin(const Packet& a, const Packet& b) {
  return pminmax_impl<NaNPropagation>::run(a, b, EIGEN_BINARY_OP_NAN_PROPAGATION(Packet, (pmin<Packet>)));
}

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pmax(const Packet& a, const Packet& b) { return numext::maxi(a, b); }

template <int NaNPropagation, typename Packet>
EIGEN_DEVICE_FUNC inline Packet pmax(const Packet& a, const Packet& b) {
  return pminmax_impl<NaNPropagation>::run(a, b, EIGEN_BINARY_OP_NAN_PROPAGATION(Packet,(pmax<Packet>)));
}

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pabs(const Packet& a) { return numext::abs(a); }
template<> EIGEN_DEVICE_FUNC inline unsigned int
pabs(const unsigned int& a) { return a; }
template<> EIGEN_DEVICE_FUNC inline unsigned long
pabs(const unsigned long& a) { return a; }
template<> EIGEN_DEVICE_FUNC inline unsigned long long
pabs(const unsigned long long& a) { return a; }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
paddsub(const Packet& a, const Packet& b) {
  return pselect(peven_mask(a), padd(a, b), psub(a, b));
 }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
parg(const Packet& a) { using numext::arg; return arg(a); }

template<int N> EIGEN_DEVICE_FUNC inline int
parithmetic_shift_right(const int& a) { return a >> N; }
template<int N> EIGEN_DEVICE_FUNC inline long int
parithmetic_shift_right(const long int& a) { return a >> N; }

template<int N> EIGEN_DEVICE_FUNC inline int
plogical_shift_right(const int& a) { return static_cast<int>(static_cast<unsigned int>(a) >> N); }
template<int N> EIGEN_DEVICE_FUNC inline long int
plogical_shift_right(const long int& a) { return static_cast<long>(static_cast<unsigned long>(a) >> N); }

template<int N> EIGEN_DEVICE_FUNC inline int
plogical_shift_left(const int& a) { return a << N; }
template<int N> EIGEN_DEVICE_FUNC inline long int
plogical_shift_left(const long int& a) { return a << N; }

template <typename Packet>
EIGEN_DEVICE_FUNC inline Packet pfrexp(const Packet& a, Packet& exponent) {
  int exp;
  EIGEN_USING_STD(frexp);
  Packet result = static_cast<Packet>(frexp(a, &exp));
  exponent = static_cast<Packet>(exp);
  return result;
}

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pldexp(const Packet &a, const Packet &exponent) {
  EIGEN_USING_STD(ldexp)
  return static_cast<Packet>(ldexp(a, static_cast<int>(exponent)));
}

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pabsdiff(const Packet& a, const Packet& b) { return pselect(pcmp_lt(a, b), psub(b, a), psub(a, b)); }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pload(const typename unpacket_traits<Packet>::type* from) { return *from; }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
ploadu(const typename unpacket_traits<Packet>::type* from) { return *from; }

template<typename Packet> EIGEN_DEVICE_FUNC inline
typename enable_if<unpacket_traits<Packet>::masked_load_available, Packet>::type
ploadu(const typename unpacket_traits<Packet>::type* from, typename unpacket_traits<Packet>::mask_t umask);

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pset1(const typename unpacket_traits<Packet>::type& a) { return a; }

template<typename Packet,typename BitsType> EIGEN_DEVICE_FUNC inline Packet
pset1frombits(BitsType a);

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pload1(const typename unpacket_traits<Packet>::type  *a) { return pset1<Packet>(*a); }

template<typename Packet> EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet
ploaddup(const typename unpacket_traits<Packet>::type* from) { return *from; }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
ploadquad(const typename unpacket_traits<Packet>::type* from)
{ return pload1<Packet>(from); }

template<typename Packet> EIGEN_DEVICE_FUNC
inline void pbroadcast4(const typename unpacket_traits<Packet>::type *a,
                        Packet& a0, Packet& a1, Packet& a2, Packet& a3)
{
  a0 = pload1<Packet>(a+0);
  a1 = pload1<Packet>(a+1);
  a2 = pload1<Packet>(a+2);
  a3 = pload1<Packet>(a+3);
}

template<typename Packet> EIGEN_DEVICE_FUNC
inline void pbroadcast2(const typename unpacket_traits<Packet>::type *a,
                        Packet& a0, Packet& a1)
{
  a0 = pload1<Packet>(a+0);
  a1 = pload1<Packet>(a+1);
}

template<typename Packet> EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet
plset(const typename unpacket_traits<Packet>::type& a) { return a; }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
peven_mask(const Packet& ) {
  typedef typename unpacket_traits<Packet>::type Scalar;
  const size_t n = unpacket_traits<Packet>::size;
  EIGEN_ALIGN_TO_BOUNDARY(sizeof(Packet)) Scalar elements[n];
  for(size_t i = 0; i < n; ++i) {
    memset(elements+i, ((i & 1) == 0 ? 0xff : 0), sizeof(Scalar));
  }
  return ploadu<Packet>(elements);
}

template<typename Scalar, typename Packet> EIGEN_DEVICE_FUNC inline void pstore(Scalar* to, const Packet& from)
{ (*to) = from; }

template<typename Scalar, typename Packet> EIGEN_DEVICE_FUNC inline void pstoreu(Scalar* to, const Packet& from)
{  (*to) = from; }

template<typename Scalar, typename Packet>
EIGEN_DEVICE_FUNC inline
typename enable_if<unpacket_traits<Packet>::masked_store_available, void>::type
pstoreu(Scalar* to, const Packet& from, typename unpacket_traits<Packet>::mask_t umask);

 template<typename Scalar, typename Packet> EIGEN_DEVICE_FUNC inline Packet pgather(const Scalar* from, Index )
 { return ploadu<Packet>(from); }

 template<typename Scalar, typename Packet> EIGEN_DEVICE_FUNC inline void pscatter(Scalar* to, const Packet& from, Index )
 { pstore(to, from); }

template<typename Scalar> EIGEN_DEVICE_FUNC inline void prefetch(const Scalar* addr)
{
#if defined(EIGEN_HIP_DEVICE_COMPILE)
  
#elif defined(EIGEN_CUDA_ARCH)
#if defined(__LP64__) || EIGEN_OS_WIN64
  
  asm(" prefetch.L1 [ %1 ];" : "=l"(addr) : "l"(addr));
#else
  
  asm(" prefetch.L1 [ %1 ];" : "=r"(addr) : "r"(addr));
#endif
#elif (!EIGEN_COMP_MSVC) && (EIGEN_COMP_GNUC || EIGEN_COMP_CLANG || EIGEN_COMP_ICC)
  __builtin_prefetch(addr);
#endif
}

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet preverse(const Packet& a)
{ return a; }

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet pcplxflip(const Packet& a)
{
  return Packet(numext::imag(a),numext::real(a));
}

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet psin(const Packet& a) { EIGEN_USING_STD(sin); return sin(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet pcos(const Packet& a) { EIGEN_USING_STD(cos); return cos(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet ptan(const Packet& a) { EIGEN_USING_STD(tan); return tan(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet pasin(const Packet& a) { EIGEN_USING_STD(asin); return asin(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet pacos(const Packet& a) { EIGEN_USING_STD(acos); return acos(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet patan(const Packet& a) { EIGEN_USING_STD(atan); return atan(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet psinh(const Packet& a) { EIGEN_USING_STD(sinh); return sinh(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet pcosh(const Packet& a) { EIGEN_USING_STD(cosh); return cosh(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet ptanh(const Packet& a) { EIGEN_USING_STD(tanh); return tanh(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet pexp(const Packet& a) { EIGEN_USING_STD(exp); return exp(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet pexpm1(const Packet& a) { return numext::expm1(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet plog(const Packet& a) { EIGEN_USING_STD(log); return log(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet plog1p(const Packet& a) { return numext::log1p(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet plog10(const Packet& a) { EIGEN_USING_STD(log10); return log10(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet plog2(const Packet& a) {
  typedef typename internal::unpacket_traits<Packet>::type Scalar;
  return pmul(pset1<Packet>(Scalar(EIGEN_LOG2E)), plog(a)); 
}

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet psqrt(const Packet& a) { return numext::sqrt(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet prsqrt(const Packet& a) {
  typedef typename internal::unpacket_traits<Packet>::type Scalar;
  return pdiv(pset1<Packet>(Scalar(1)), psqrt(a));
}

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet pround(const Packet& a) { using numext::round; return round(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet pfloor(const Packet& a) { using numext::floor; return floor(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet print(const Packet& a) { using numext::rint; return rint(a); }

template<typename Packet> EIGEN_DECLARE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
Packet pceil(const Packet& a) { using numext::ceil; return ceil(a); }

template<typename Packet>
EIGEN_DEVICE_FUNC inline typename unpacket_traits<Packet>::type
pfirst(const Packet& a)
{ return a; }

template<typename Packet>
EIGEN_DEVICE_FUNC inline typename conditional<(unpacket_traits<Packet>::size%8)==0,typename unpacket_traits<Packet>::half,Packet>::type
predux_half_dowto4(const Packet& a)
{ return a; }

template <typename Packet, typename Op>
EIGEN_DEVICE_FUNC inline typename unpacket_traits<Packet>::type
predux_helper(const Packet& a, Op op) {
  typedef typename unpacket_traits<Packet>::type Scalar;
  const size_t n = unpacket_traits<Packet>::size;
  EIGEN_ALIGN_TO_BOUNDARY(sizeof(Packet)) Scalar elements[n];
  pstoreu<Scalar>(elements, a);
  for(size_t k = n / 2; k > 0; k /= 2)  {
    for(size_t i = 0; i < k; ++i) {
      elements[i] = op(elements[i], elements[i + k]);
    }
  }
  return elements[0];
}

template<typename Packet>
EIGEN_DEVICE_FUNC inline typename unpacket_traits<Packet>::type
predux(const Packet& a)
{
  return a;
}

template <typename Packet>
EIGEN_DEVICE_FUNC inline typename unpacket_traits<Packet>::type predux_mul(
    const Packet& a) {
  typedef typename unpacket_traits<Packet>::type Scalar; 
  return predux_helper(a, EIGEN_BINARY_OP_NAN_PROPAGATION(Scalar, (pmul<Scalar>)));
}

template <typename Packet>
EIGEN_DEVICE_FUNC inline typename unpacket_traits<Packet>::type predux_min(
    const Packet &a) {
  typedef typename unpacket_traits<Packet>::type Scalar; 
  return predux_helper(a, EIGEN_BINARY_OP_NAN_PROPAGATION(Scalar, (pmin<PropagateFast, Scalar>)));
}

template <int NaNPropagation, typename Packet>
EIGEN_DEVICE_FUNC inline typename unpacket_traits<Packet>::type predux_min(
    const Packet& a) {
  typedef typename unpacket_traits<Packet>::type Scalar; 
  return predux_helper(a, EIGEN_BINARY_OP_NAN_PROPAGATION(Scalar, (pmin<NaNPropagation, Scalar>)));
}

template <typename Packet>
EIGEN_DEVICE_FUNC inline typename unpacket_traits<Packet>::type predux_max(
    const Packet &a) {
  typedef typename unpacket_traits<Packet>::type Scalar; 
  return predux_helper(a, EIGEN_BINARY_OP_NAN_PROPAGATION(Scalar, (pmax<PropagateFast, Scalar>)));
}

template <int NaNPropagation, typename Packet>
EIGEN_DEVICE_FUNC inline typename unpacket_traits<Packet>::type predux_max(
    const Packet& a) {
  typedef typename unpacket_traits<Packet>::type Scalar; 
  return predux_helper(a, EIGEN_BINARY_OP_NAN_PROPAGATION(Scalar, (pmax<NaNPropagation, Scalar>)));
}

#undef EIGEN_BINARY_OP_NAN_PROPAGATION

template<typename Packet> EIGEN_DEVICE_FUNC inline bool predux_any(const Packet& a)
{

  typedef typename unpacket_traits<Packet>::type Scalar;
  return numext::not_equal_strict(predux(a), Scalar(0));
}

template<typename Packet>
inline void pstore1(typename unpacket_traits<Packet>::type* to, const typename unpacket_traits<Packet>::type& a)
{
  pstore(to, pset1<Packet>(a));
}

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pmadd(const Packet&  a,
         const Packet&  b,
         const Packet&  c)
{ return padd(pmul(a, b),c); }

template<typename Packet, int Alignment>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE Packet ploadt(const typename unpacket_traits<Packet>::type* from)
{
  if(Alignment >= unpacket_traits<Packet>::alignment)
    return pload<Packet>(from);
  else
    return ploadu<Packet>(from);
}

template<typename Scalar, typename Packet, int Alignment>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE void pstoret(Scalar* to, const Packet& from)
{
  if(Alignment >= unpacket_traits<Packet>::alignment)
    pstore(to, from);
  else
    pstoreu(to, from);
}

template<typename Packet, int LoadMode>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE Packet ploadt_ro(const typename unpacket_traits<Packet>::type* from)
{
  return ploadt<Packet, LoadMode>(from);
}

#if !defined(EIGEN_GPUCC)

template<> inline std::complex<float> pmul(const std::complex<float>& a, const std::complex<float>& b)
{ return std::complex<float>(a.real()*b.real() - a.imag()*b.imag(), a.imag()*b.real() + a.real()*b.imag()); }

template<> inline std::complex<double> pmul(const std::complex<double>& a, const std::complex<double>& b)
{ return std::complex<double>(a.real()*b.real() - a.imag()*b.imag(), a.imag()*b.real() + a.real()*b.imag()); }

#endif

template <typename Packet,int N=unpacket_traits<Packet>::size> struct PacketBlock {
  Packet packet[N];
};

template<typename Packet> EIGEN_DEVICE_FUNC inline void
ptranspose(PacketBlock<Packet,1>& ) {
  
}

template <size_t N> struct Selector {
  bool select[N];
};

template<typename Packet> EIGEN_DEVICE_FUNC inline Packet
pblend(const Selector<unpacket_traits<Packet>::size>& ifPacket, const Packet& thenPacket, const Packet& elsePacket) {
  return ifPacket.select[0] ? thenPacket : elsePacket;
}

} 

} 

#endif 
