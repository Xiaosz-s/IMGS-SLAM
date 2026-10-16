#ifndef EIGEN_ARCH_GENERIC_PACKET_MATH_FUNCTIONS_H
#define EIGEN_ARCH_GENERIC_PACKET_MATH_FUNCTIONS_H

namespace Eigen {
namespace internal {

template<typename T> struct make_integer;
template<> struct make_integer<float>    { typedef numext::int32_t type; };
template<> struct make_integer<double>   { typedef numext::int64_t type; };
template<> struct make_integer<half>     { typedef numext::int16_t type; };
template<> struct make_integer<bfloat16> { typedef numext::int16_t type; };

template<typename Packet> EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC  
Packet pfrexp_generic_get_biased_exponent(const Packet& a) {
  typedef typename unpacket_traits<Packet>::type Scalar;
  typedef typename unpacket_traits<Packet>::integer_packet PacketI;
  enum { mantissa_bits = numext::numeric_limits<Scalar>::digits - 1};
  return pcast<PacketI, Packet>(plogical_shift_right<mantissa_bits>(preinterpret<PacketI>(pabs(a))));
}

template<typename Packet> EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC
Packet pfrexp_generic(const Packet& a, Packet& exponent) {
  typedef typename unpacket_traits<Packet>::type Scalar;
  typedef typename make_unsigned<typename make_integer<Scalar>::type>::type ScalarUI;
  enum {
    TotalBits = sizeof(Scalar) * CHAR_BIT,
    MantissaBits = numext::numeric_limits<Scalar>::digits - 1,
    ExponentBits = int(TotalBits) - int(MantissaBits) - 1
  };

  EIGEN_CONSTEXPR ScalarUI scalar_sign_mantissa_mask = 
      ~(((ScalarUI(1) << int(ExponentBits)) - ScalarUI(1)) << int(MantissaBits)); 
  const Packet sign_mantissa_mask = pset1frombits<Packet>(static_cast<ScalarUI>(scalar_sign_mantissa_mask)); 
  const Packet half = pset1<Packet>(Scalar(0.5));
  const Packet zero = pzero(a);
  const Packet normal_min = pset1<Packet>((numext::numeric_limits<Scalar>::min)()); 

  const Packet is_denormal = pcmp_lt(pabs(a), normal_min);
  EIGEN_CONSTEXPR ScalarUI scalar_normalization_offset = ScalarUI(int(MantissaBits) + 1); 
  
  const Scalar scalar_normalization_factor = Scalar(ScalarUI(1) << int(scalar_normalization_offset)); 
  const Packet normalization_factor = pset1<Packet>(scalar_normalization_factor);  
  const Packet normalized_a = pselect(is_denormal, pmul(a, normalization_factor), a);

  const Scalar scalar_exponent_offset = -Scalar((ScalarUI(1)<<(int(ExponentBits)-1)) - ScalarUI(2)); 
  Packet exponent_offset = pset1<Packet>(scalar_exponent_offset);
  const Packet normalization_offset = pset1<Packet>(-Scalar(scalar_normalization_offset)); 
  exponent_offset = pselect(is_denormal, padd(exponent_offset, normalization_offset), exponent_offset);

  exponent = pfrexp_generic_get_biased_exponent(normalized_a);

  const Scalar scalar_non_finite_exponent = Scalar((ScalarUI(1) << int(ExponentBits)) - ScalarUI(1));  
  const Packet non_finite_exponent = pset1<Packet>(scalar_non_finite_exponent);
  const Packet is_zero_or_not_finite = por(pcmp_eq(a, zero), pcmp_eq(exponent, non_finite_exponent));
  const Packet m = pselect(is_zero_or_not_finite, a, por(pand(normalized_a, sign_mantissa_mask), half));
  exponent = pselect(is_zero_or_not_finite, zero, padd(exponent, exponent_offset));  
  return m;
}

template<typename Packet> EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC
Packet pldexp_generic(const Packet& a, const Packet& exponent) {

  typedef typename unpacket_traits<Packet>::integer_packet PacketI;
  typedef typename unpacket_traits<Packet>::type Scalar;
  typedef typename unpacket_traits<PacketI>::type ScalarI;
  enum {
    TotalBits = sizeof(Scalar) * CHAR_BIT,
    MantissaBits = numext::numeric_limits<Scalar>::digits - 1,
    ExponentBits = int(TotalBits) - int(MantissaBits) - 1
  };

  const Packet max_exponent = pset1<Packet>(Scalar((ScalarI(1)<<int(ExponentBits)) + ScalarI(int(MantissaBits) - 1)));  
  const PacketI bias = pset1<PacketI>((ScalarI(1)<<(int(ExponentBits)-1)) - ScalarI(1));  
  const PacketI e = pcast<Packet, PacketI>(pmin(pmax(exponent, pnegate(max_exponent)), max_exponent));
  PacketI b = parithmetic_shift_right<2>(e); 
  Packet c = preinterpret<Packet>(plogical_shift_left<int(MantissaBits)>(padd(b, bias)));  
  Packet out = pmul(pmul(pmul(a, c), c), c);  
  b = psub(psub(psub(e, b), b), b); 
  c = preinterpret<Packet>(plogical_shift_left<int(MantissaBits)>(padd(b, bias)));  
  out = pmul(out, c);
  return out;
}

template<typename Packet>
struct pldexp_fast_impl {
  typedef typename unpacket_traits<Packet>::integer_packet PacketI;
  typedef typename unpacket_traits<Packet>::type Scalar;
  typedef typename unpacket_traits<PacketI>::type ScalarI;
  enum {
    TotalBits = sizeof(Scalar) * CHAR_BIT,
    MantissaBits = numext::numeric_limits<Scalar>::digits - 1,
    ExponentBits = int(TotalBits) - int(MantissaBits) - 1
  };
  
  static EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC
  Packet run(const Packet& a, const Packet& exponent) {
    const Packet bias = pset1<Packet>(Scalar((ScalarI(1)<<(int(ExponentBits)-1)) - ScalarI(1)));  
    const Packet limit = pset1<Packet>(Scalar((ScalarI(1)<<int(ExponentBits)) - ScalarI(1)));     
    
    const PacketI e = pcast<Packet, PacketI>(pmin(pmax(padd(exponent, bias), pzero(limit)), limit)); 
    
    return pmul(a, preinterpret<Packet>(plogical_shift_left<int(MantissaBits)>(e)));
  }
};

template <typename Packet, bool base2>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet plog_impl_float(const Packet _x)
{
  Packet x = _x;

  const Packet cst_1              = pset1<Packet>(1.0f);
  const Packet cst_neg_half       = pset1<Packet>(-0.5f);
  
  const Packet cst_min_norm_pos   = pset1frombits<Packet>( 0x00800000u);
  const Packet cst_minus_inf      = pset1frombits<Packet>( 0xff800000u);
  const Packet cst_pos_inf        = pset1frombits<Packet>( 0x7f800000u);

  const Packet cst_cephes_SQRTHF = pset1<Packet>(0.707106781186547524f);
  const Packet cst_cephes_log_p0 = pset1<Packet>(7.0376836292E-2f);
  const Packet cst_cephes_log_p1 = pset1<Packet>(-1.1514610310E-1f);
  const Packet cst_cephes_log_p2 = pset1<Packet>(1.1676998740E-1f);
  const Packet cst_cephes_log_p3 = pset1<Packet>(-1.2420140846E-1f);
  const Packet cst_cephes_log_p4 = pset1<Packet>(+1.4249322787E-1f);
  const Packet cst_cephes_log_p5 = pset1<Packet>(-1.6668057665E-1f);
  const Packet cst_cephes_log_p6 = pset1<Packet>(+2.0000714765E-1f);
  const Packet cst_cephes_log_p7 = pset1<Packet>(-2.4999993993E-1f);
  const Packet cst_cephes_log_p8 = pset1<Packet>(+3.3333331174E-1f);

  x = pmax(x, cst_min_norm_pos);

  Packet e;
  
  x = pfrexp(x,e);

  Packet mask = pcmp_lt(x, cst_cephes_SQRTHF);
  Packet tmp = pand(x, mask);
  x = psub(x, cst_1);
  e = psub(e, pand(cst_1, mask));
  x = padd(x, tmp);

  Packet x2 = pmul(x, x);
  Packet x3 = pmul(x2, x);

  Packet y, y1, y2;
  y  = pmadd(cst_cephes_log_p0, x, cst_cephes_log_p1);
  y1 = pmadd(cst_cephes_log_p3, x, cst_cephes_log_p4);
  y2 = pmadd(cst_cephes_log_p6, x, cst_cephes_log_p7);
  y  = pmadd(y, x, cst_cephes_log_p2);
  y1 = pmadd(y1, x, cst_cephes_log_p5);
  y2 = pmadd(y2, x, cst_cephes_log_p8);
  y  = pmadd(y, x3, y1);
  y  = pmadd(y, x3, y2);
  y  = pmul(y, x3);

  y = pmadd(cst_neg_half, x2, y);
  x = padd(x, y);

  if (base2) {
    const Packet cst_log2e = pset1<Packet>(static_cast<float>(EIGEN_LOG2E));
    x = pmadd(x, cst_log2e, e);
  } else {
    const Packet cst_ln2 = pset1<Packet>(static_cast<float>(EIGEN_LN2));
    x = pmadd(e, cst_ln2, x);
  }

  Packet invalid_mask = pcmp_lt_or_nan(_x, pzero(_x));
  Packet iszero_mask  = pcmp_eq(_x,pzero(_x));
  Packet pos_inf_mask = pcmp_eq(_x,cst_pos_inf);

  return pselect(iszero_mask, cst_minus_inf,
                              por(pselect(pos_inf_mask,cst_pos_inf,x), invalid_mask));
}

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet plog_float(const Packet _x)
{
  return plog_impl_float<Packet,  false>(_x);
}

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet plog2_float(const Packet _x)
{
  return plog_impl_float<Packet,  true>(_x);
}

template <typename Packet, bool base2>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet plog_impl_double(const Packet _x)
{
  Packet x = _x;

  const Packet cst_1              = pset1<Packet>(1.0);
  const Packet cst_neg_half       = pset1<Packet>(-0.5);
  
  const Packet cst_min_norm_pos   = pset1frombits<Packet>( static_cast<uint64_t>(0x0010000000000000ull));
  const Packet cst_minus_inf      = pset1frombits<Packet>( static_cast<uint64_t>(0xfff0000000000000ull));
  const Packet cst_pos_inf        = pset1frombits<Packet>( static_cast<uint64_t>(0x7ff0000000000000ull));

  const Packet cst_cephes_SQRTHF = pset1<Packet>(0.70710678118654752440E0);
  const Packet cst_cephes_log_p0 = pset1<Packet>(1.01875663804580931796E-4);
  const Packet cst_cephes_log_p1 = pset1<Packet>(4.97494994976747001425E-1);
  const Packet cst_cephes_log_p2 = pset1<Packet>(4.70579119878881725854E0);
  const Packet cst_cephes_log_p3 = pset1<Packet>(1.44989225341610930846E1);
  const Packet cst_cephes_log_p4 = pset1<Packet>(1.79368678507819816313E1);
  const Packet cst_cephes_log_p5 = pset1<Packet>(7.70838733755885391666E0);

  const Packet cst_cephes_log_q0 = pset1<Packet>(1.0);
  const Packet cst_cephes_log_q1 = pset1<Packet>(1.12873587189167450590E1);
  const Packet cst_cephes_log_q2 = pset1<Packet>(4.52279145837532221105E1);
  const Packet cst_cephes_log_q3 = pset1<Packet>(8.29875266912776603211E1);
  const Packet cst_cephes_log_q4 = pset1<Packet>(7.11544750618563894466E1);
  const Packet cst_cephes_log_q5 = pset1<Packet>(2.31251620126765340583E1);

  x = pmax(x, cst_min_norm_pos);

  Packet e;
  
  x = pfrexp(x,e);

  Packet mask = pcmp_lt(x, cst_cephes_SQRTHF);
  Packet tmp = pand(x, mask);
  x = psub(x, cst_1);
  e = psub(e, pand(cst_1, mask));
  x = padd(x, tmp);

  Packet x2 = pmul(x, x);
  Packet x3 = pmul(x2, x);

  Packet y, y1, y_;
  y  = pmadd(cst_cephes_log_p0, x, cst_cephes_log_p1);
  y1 = pmadd(cst_cephes_log_p3, x, cst_cephes_log_p4);
  y  = pmadd(y, x, cst_cephes_log_p2);
  y1 = pmadd(y1, x, cst_cephes_log_p5);
  y_ = pmadd(y, x3, y1);

  y  = pmadd(cst_cephes_log_q0, x, cst_cephes_log_q1);
  y1 = pmadd(cst_cephes_log_q3, x, cst_cephes_log_q4);
  y  = pmadd(y, x, cst_cephes_log_q2);
  y1 = pmadd(y1, x, cst_cephes_log_q5);
  y  = pmadd(y, x3, y1);

  y_ = pmul(y_, x3);
  y  = pdiv(y_, y);

  y = pmadd(cst_neg_half, x2, y);
  x = padd(x, y);

  if (base2) {
    const Packet cst_log2e = pset1<Packet>(static_cast<double>(EIGEN_LOG2E));
    x = pmadd(x, cst_log2e, e);
  } else {
    const Packet cst_ln2 = pset1<Packet>(static_cast<double>(EIGEN_LN2));
    x = pmadd(e, cst_ln2, x);
  }

  Packet invalid_mask = pcmp_lt_or_nan(_x, pzero(_x));
  Packet iszero_mask  = pcmp_eq(_x,pzero(_x));
  Packet pos_inf_mask = pcmp_eq(_x,cst_pos_inf);

  return pselect(iszero_mask, cst_minus_inf,
                              por(pselect(pos_inf_mask,cst_pos_inf,x), invalid_mask));
}

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet plog_double(const Packet _x)
{
  return plog_impl_double<Packet,  false>(_x);
}

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet plog2_double(const Packet _x)
{
  return plog_impl_double<Packet,  true>(_x);
}

template<typename Packet>
Packet generic_plog1p(const Packet& x)
{
  typedef typename unpacket_traits<Packet>::type ScalarType;
  const Packet one = pset1<Packet>(ScalarType(1));
  Packet xp1 = padd(x, one);
  Packet small_mask = pcmp_eq(xp1, one);
  Packet log1 = plog(xp1);
  Packet inf_mask = pcmp_eq(xp1, log1);
  Packet log_large = pmul(x, pdiv(log1, psub(xp1, one)));
  return pselect(por(small_mask, inf_mask), x, log_large);
}

template<typename Packet>
Packet generic_expm1(const Packet& x)
{
  typedef typename unpacket_traits<Packet>::type ScalarType;
  const Packet one = pset1<Packet>(ScalarType(1));
  const Packet neg_one = pset1<Packet>(ScalarType(-1));
  Packet u = pexp(x);
  Packet one_mask = pcmp_eq(u, one);
  Packet u_minus_one = psub(u, one);
  Packet neg_one_mask = pcmp_eq(u_minus_one, neg_one);
  Packet logu = plog(u);

  Packet pos_inf_mask = pcmp_eq(logu, u);
  Packet expm1 = pmul(u_minus_one, pdiv(x, logu));
  expm1 = pselect(pos_inf_mask, u, expm1);
  return pselect(one_mask,
                 x,
                 pselect(neg_one_mask,
                         neg_one,
                         expm1));
}

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet pexp_float(const Packet _x)
{
  const Packet cst_1      = pset1<Packet>(1.0f);
  const Packet cst_half   = pset1<Packet>(0.5f);
  const Packet cst_exp_hi = pset1<Packet>( 88.723f);
  const Packet cst_exp_lo = pset1<Packet>(-88.723f);

  const Packet cst_cephes_LOG2EF = pset1<Packet>(1.44269504088896341f);
  const Packet cst_cephes_exp_p0 = pset1<Packet>(1.9875691500E-4f);
  const Packet cst_cephes_exp_p1 = pset1<Packet>(1.3981999507E-3f);
  const Packet cst_cephes_exp_p2 = pset1<Packet>(8.3334519073E-3f);
  const Packet cst_cephes_exp_p3 = pset1<Packet>(4.1665795894E-2f);
  const Packet cst_cephes_exp_p4 = pset1<Packet>(1.6666665459E-1f);
  const Packet cst_cephes_exp_p5 = pset1<Packet>(5.0000001201E-1f);

  Packet x = pmax(pmin(_x, cst_exp_hi), cst_exp_lo);

  Packet m = pfloor(pmadd(x, cst_cephes_LOG2EF, cst_half));

  const Packet cst_cephes_exp_C1 = pset1<Packet>(-0.693359375f);
  const Packet cst_cephes_exp_C2 = pset1<Packet>(2.12194440e-4f);
  Packet r = pmadd(m, cst_cephes_exp_C1, x);
  r = pmadd(m, cst_cephes_exp_C2, r);

  Packet r2 = pmul(r, r);
  Packet r3 = pmul(r2, r);

  Packet y, y1, y2;
  y  = pmadd(cst_cephes_exp_p0, r, cst_cephes_exp_p1);
  y1 = pmadd(cst_cephes_exp_p3, r, cst_cephes_exp_p4);
  y2 = padd(r, cst_1);
  y  = pmadd(y, r, cst_cephes_exp_p2);
  y1 = pmadd(y1, r, cst_cephes_exp_p5);
  y  = pmadd(y, r3, y1);
  y  = pmadd(y, r2, y2);

  return pmax(pldexp(y,m), _x);
}

template <typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet pexp_double(const Packet _x)
{
  Packet x = _x;

  const Packet cst_1 = pset1<Packet>(1.0);
  const Packet cst_2 = pset1<Packet>(2.0);
  const Packet cst_half = pset1<Packet>(0.5);

  const Packet cst_exp_hi = pset1<Packet>(709.784);
  const Packet cst_exp_lo = pset1<Packet>(-709.784);

  const Packet cst_cephes_LOG2EF = pset1<Packet>(1.4426950408889634073599);
  const Packet cst_cephes_exp_p0 = pset1<Packet>(1.26177193074810590878e-4);
  const Packet cst_cephes_exp_p1 = pset1<Packet>(3.02994407707441961300e-2);
  const Packet cst_cephes_exp_p2 = pset1<Packet>(9.99999999999999999910e-1);
  const Packet cst_cephes_exp_q0 = pset1<Packet>(3.00198505138664455042e-6);
  const Packet cst_cephes_exp_q1 = pset1<Packet>(2.52448340349684104192e-3);
  const Packet cst_cephes_exp_q2 = pset1<Packet>(2.27265548208155028766e-1);
  const Packet cst_cephes_exp_q3 = pset1<Packet>(2.00000000000000000009e0);
  const Packet cst_cephes_exp_C1 = pset1<Packet>(0.693145751953125);
  const Packet cst_cephes_exp_C2 = pset1<Packet>(1.42860682030941723212e-6);

  Packet tmp, fx;

  x = pmax(pmin(x, cst_exp_hi), cst_exp_lo);
  
  fx = pmadd(cst_cephes_LOG2EF, x, cst_half);

  fx = pfloor(fx);

  tmp = pmul(fx, cst_cephes_exp_C1);
  Packet z = pmul(fx, cst_cephes_exp_C2);
  x = psub(x, tmp);
  x = psub(x, z);

  Packet x2 = pmul(x, x);

  Packet px = cst_cephes_exp_p0;
  px = pmadd(px, x2, cst_cephes_exp_p1);
  px = pmadd(px, x2, cst_cephes_exp_p2);
  px = pmul(px, x);

  Packet qx = cst_cephes_exp_q0;
  qx = pmadd(qx, x2, cst_cephes_exp_q1);
  qx = pmadd(qx, x2, cst_cephes_exp_q2);
  qx = pmadd(qx, x2, cst_cephes_exp_q3);

  x = pdiv(px, psub(qx, px));
  x = pmadd(cst_2, x, cst_1);

  return pmax(pldexp(x,fx), _x);
}

inline float trig_reduce_huge (float xf, int *quadrant)
{
  using Eigen::numext::int32_t;
  using Eigen::numext::uint32_t;
  using Eigen::numext::int64_t;
  using Eigen::numext::uint64_t;

  const double pio2_62 = 3.4061215800865545e-19;    
  const uint64_t zero_dot_five = uint64_t(1) << 61; 

  static const uint32_t two_over_pi [] = 
  {
    0x00000028, 0x000028be, 0x0028be60, 0x28be60db,
    0xbe60db93, 0x60db9391, 0xdb939105, 0x9391054a,
    0x91054a7f, 0x054a7f09, 0x4a7f09d5, 0x7f09d5f4,
    0x09d5f47d, 0xd5f47d4d, 0xf47d4d37, 0x7d4d3770,
    0x4d377036, 0x377036d8, 0x7036d8a5, 0x36d8a566,
    0xd8a5664f, 0xa5664f10, 0x664f10e4, 0x4f10e410,
    0x10e41000, 0xe4100000
  };
  
  uint32_t xi = numext::bit_cast<uint32_t>(xf);

  uint32_t e = (xi >> 23) - 118;
  
  xi = ((xi & 0x007fffffu)| 0x00800000u) << (e & 0x7);

  uint32_t i = e >> 3;
  uint32_t twoopi_1  = two_over_pi[i-1];
  uint32_t twoopi_2  = two_over_pi[i+3];
  uint32_t twoopi_3  = two_over_pi[i+7];

  uint64_t p;
  p = uint64_t(xi) * twoopi_3;
  p = uint64_t(xi) * twoopi_2 + (p >> 32);
  p = (uint64_t(xi * twoopi_1) << 32) + p;

  uint64_t q = (p + zero_dot_five) >> 62;
  *quadrant = int(q);

  p -= q<<62;
  return float(double(int64_t(p)) * pio2_62);
}

template<bool ComputeSine,typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
#if EIGEN_GNUC_AT_LEAST(4,4) && EIGEN_COMP_GNUC_STRICT
__attribute__((optimize("-fno-unsafe-math-optimizations")))
#endif
Packet psincos_float(const Packet& _x)
{
  typedef typename unpacket_traits<Packet>::integer_packet PacketI;

  const Packet  cst_2oPI            = pset1<Packet>(0.636619746685028076171875f); 
  const Packet  cst_rounding_magic  = pset1<Packet>(12582912); 
  const PacketI csti_1              = pset1<PacketI>(1);
  const Packet  cst_sign_mask       = pset1frombits<Packet>(0x80000000u);

  Packet x = pabs(_x);

  Packet y = pmul(x, cst_2oPI);

  Packet y_round = padd(y, cst_rounding_magic);
  EIGEN_OPTIMIZATION_BARRIER(y_round)
  PacketI y_int = preinterpret<PacketI>(y_round); 
  y = psub(y_round, cst_rounding_magic); 

  #if defined(EIGEN_HAS_SINGLE_INSTRUCTION_MADD)

  const float huge_th = ComputeSine ? 117435.992f : 71476.0625f;
  x = pmadd(y, pset1<Packet>(-1.57079601287841796875f), x);
  x = pmadd(y, pset1<Packet>(-3.1391647326017846353352069854736328125e-07f), x);
  x = pmadd(y, pset1<Packet>(-5.390302529957764765544681040410068817436695098876953125e-15f), x);
  #else

  const float huge_th = ComputeSine ? 25966.f : 18838.f;
  x = pmadd(y, pset1<Packet>(-1.5703125), x); 
  EIGEN_OPTIMIZATION_BARRIER(x)
  x = pmadd(y, pset1<Packet>(-0.000483989715576171875), x); 
  EIGEN_OPTIMIZATION_BARRIER(x)
  x = pmadd(y, pset1<Packet>(1.62865035235881805419921875e-07), x); 
  x = pmadd(y, pset1<Packet>(5.5644315544167710640977020375430583953857421875e-11), x); 

  #endif

  if(predux_any(pcmp_le(pset1<Packet>(huge_th),pabs(_x))))
  {
    const int PacketSize = unpacket_traits<Packet>::size;
    EIGEN_ALIGN_TO_BOUNDARY(sizeof(Packet)) float vals[PacketSize];
    EIGEN_ALIGN_TO_BOUNDARY(sizeof(Packet)) float x_cpy[PacketSize];
    EIGEN_ALIGN_TO_BOUNDARY(sizeof(Packet)) int y_int2[PacketSize];
    pstoreu(vals, pabs(_x));
    pstoreu(x_cpy, x);
    pstoreu(y_int2, y_int);
    for(int k=0; k<PacketSize;++k)
    {
      float val = vals[k];
      if(val>=huge_th && (numext::isfinite)(val))
        x_cpy[k] = trig_reduce_huge(val,&y_int2[k]);
    }
    x = ploadu<Packet>(x_cpy);
    y_int = ploadu<PacketI>(y_int2);
  }

  Packet sign_bit = ComputeSine ? pxor(_x, preinterpret<Packet>(plogical_shift_left<30>(y_int)))
                                : preinterpret<Packet>(plogical_shift_left<30>(padd(y_int,csti_1)));
  sign_bit = pand(sign_bit, cst_sign_mask); 

  Packet poly_mask = preinterpret<Packet>(pcmp_eq(pand(y_int, csti_1), pzero(y_int)));

  Packet x2 = pmul(x,x);

  Packet y1 =        pset1<Packet>(2.4372266125283204019069671630859375e-05f);
  y1 = pmadd(y1, x2, pset1<Packet>(-0.00138865201734006404876708984375f     ));
  y1 = pmadd(y1, x2, pset1<Packet>(0.041666619479656219482421875f           ));
  y1 = pmadd(y1, x2, pset1<Packet>(-0.5f));
  y1 = pmadd(y1, x2, pset1<Packet>(1.f));

  Packet y2 =        pset1<Packet>(-0.0001959234114083702898469196984621021329076029360294342041015625f);
  y2 = pmadd(y2, x2, pset1<Packet>( 0.0083326873655616851693794799871284340042620897293090820312500000f));
  y2 = pmadd(y2, x2, pset1<Packet>(-0.1666666203982298255503735617821803316473960876464843750000000000f));
  y2 = pmul(y2, x2);
  y2 = pmadd(y2, x, x);

  y = ComputeSine ? pselect(poly_mask,y2,y1)
                  : pselect(poly_mask,y1,y2);

  return pxor(y, sign_bit);
}

template<typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet psin_float(const Packet& x)
{
  return psincos_float<true>(x);
}

template<typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet pcos_float(const Packet& x)
{
  return psincos_float<false>(x);
}

template<typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED Packet pdiv_complex(const Packet& x, const Packet& y) {
  typedef typename unpacket_traits<Packet>::as_real RealPacket;

  const RealPacket y_abs = pabs(y.v);  
  const RealPacket y_abs_flip = pcplxflip(Packet(y_abs)).v; 
  const RealPacket y_max = pmax(y_abs, y_abs_flip); 
  const RealPacket y_scaled = pdiv(y.v, y_max);  
  
  const RealPacket y_scaled_sq = pmul(y_scaled, y_scaled); 
  const RealPacket denom = padd(y_scaled_sq, pcplxflip(Packet(y_scaled_sq)).v);
  Packet result_scaled = pmul(x, pconj(Packet(y_scaled)));  
  
  result_scaled = Packet(pdiv(result_scaled.v, denom));
  
  return Packet(pdiv(result_scaled.v, y_max));
}

template<typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet psqrt_complex(const Packet& a) {
  typedef typename unpacket_traits<Packet>::type Scalar;
  typedef typename Scalar::value_type RealScalar;
  typedef typename unpacket_traits<Packet>::as_real RealPacket;

  RealPacket a_abs = pabs(a.v);           
  RealPacket a_abs_flip = pcplxflip(Packet(a_abs)).v; 
  RealPacket a_max = pmax(a_abs, a_abs_flip);
  RealPacket a_min = pmin(a_abs, a_abs_flip);
  RealPacket a_min_zero_mask = pcmp_eq(a_min, pzero(a_min));
  RealPacket a_max_zero_mask = pcmp_eq(a_max, pzero(a_max));
  RealPacket r = pdiv(a_min, a_max);
  const RealPacket cst_one  = pset1<RealPacket>(RealScalar(1));
  RealPacket l = pmul(a_max, psqrt(padd(cst_one, pmul(r, r))));  
  
  l = pselect(a_min_zero_mask, a_max, l);

  const RealPacket cst_half = pset1<RealPacket>(RealScalar(0.5));
  Packet rho;
  rho.v = psqrt(pmul(cst_half, padd(a_abs, l)));

  RealPacket eta = pandnot(pmul(cst_half, pdiv(a.v, pcplxflip(rho).v)), a_max_zero_mask);
  RealPacket real_mask = peven_mask(a.v);
  Packet positive_real_result;
  
  positive_real_result.v = pselect(real_mask, rho.v, eta);

  const RealScalar neg_zero = RealScalar(numext::bit_cast<float>(0x80000000u));
  const RealPacket cst_imag_sign_mask = pset1<Packet>(Scalar(RealScalar(0.0), neg_zero)).v;
  RealPacket imag_signs = pand(a.v, cst_imag_sign_mask);
  Packet negative_real_result;
  
  negative_real_result.v = por(pabs(pcplxflip(positive_real_result).v), imag_signs);

  Packet negative_real_mask;
  negative_real_mask.v = pcmp_lt(pand(real_mask, a.v), pzero(a.v));
  negative_real_mask.v = por(negative_real_mask.v, pcplxflip(negative_real_mask).v);
  Packet result = pselect(negative_real_mask, negative_real_result, positive_real_result);

  const RealPacket cst_pos_inf = pset1<RealPacket>(NumTraits<RealScalar>::infinity());
  Packet is_inf;
  is_inf.v = pcmp_eq(a_abs, cst_pos_inf);
  Packet is_real_inf;
  is_real_inf.v = pand(is_inf.v, real_mask);
  is_real_inf = por(is_real_inf, pcplxflip(is_real_inf));
  
  Packet real_inf_result;
  real_inf_result.v = pmul(a_abs, pset1<Packet>(Scalar(RealScalar(1.0), RealScalar(0.0))).v);
  real_inf_result.v = pselect(negative_real_mask.v, pcplxflip(real_inf_result).v, real_inf_result.v);
  
  Packet is_imag_inf;
  is_imag_inf.v = pandnot(is_inf.v, real_mask);
  is_imag_inf = por(is_imag_inf, pcplxflip(is_imag_inf));
  Packet imag_inf_result;
  imag_inf_result.v = por(pand(cst_pos_inf, real_mask), pandnot(a.v, real_mask));

  return  pselect(is_imag_inf, imag_inf_result,
                  pselect(is_real_inf, real_inf_result,result));
}

template<typename Packet>
EIGEN_STRONG_INLINE
void absolute_split(const Packet& x, Packet& n, Packet& r) {
  n = pround(x);
  r = psub(x, n);
}

template<typename Packet>
EIGEN_STRONG_INLINE
void fast_twosum(const Packet& x, const Packet& y, Packet& s_hi, Packet& s_lo) {
  s_hi = padd(x, y);
  const Packet t = psub(s_hi, x);
  s_lo = psub(y, t);
}

#ifdef EIGEN_HAS_SINGLE_INSTRUCTION_MADD

template<typename Packet>
EIGEN_STRONG_INLINE
void twoprod(const Packet& x, const Packet& y,
             Packet& p_hi, Packet& p_lo) {
  p_hi = pmul(x, y);
  p_lo = pmadd(x, y, pnegate(p_hi));
}

#else

template<typename Packet>
EIGEN_STRONG_INLINE
void veltkamp_splitting(const Packet& x, Packet& x_hi, Packet& x_lo) {
  typedef typename unpacket_traits<Packet>::type Scalar;
  EIGEN_CONSTEXPR int shift = (NumTraits<Scalar>::digits() + 1) / 2;
  const Scalar shift_scale = Scalar(uint64_t(1) << shift);  
  const Packet gamma = pmul(pset1<Packet>(shift_scale + Scalar(1)), x);
  Packet rho = psub(x, gamma);
  x_hi = padd(rho, gamma);
  x_lo = psub(x, x_hi);
}

template<typename Packet>
EIGEN_STRONG_INLINE
void twoprod(const Packet& x, const Packet& y,
             Packet& p_hi, Packet& p_lo) {
  Packet x_hi, x_lo, y_hi, y_lo;
  veltkamp_splitting(x, x_hi, x_lo);
  veltkamp_splitting(y, y_hi, y_lo);

  p_hi = pmul(x, y);
  p_lo = pmadd(x_hi, y_hi, pnegate(p_hi));
  p_lo = pmadd(x_hi, y_lo, p_lo);
  p_lo = pmadd(x_lo, y_hi, p_lo);
  p_lo = pmadd(x_lo, y_lo, p_lo);
}

#endif  

template<typename Packet>
EIGEN_STRONG_INLINE
  void twosum(const Packet& x_hi, const Packet& x_lo,
              const Packet& y_hi, const Packet& y_lo,
              Packet& s_hi, Packet& s_lo) {
  const Packet x_greater_mask = pcmp_lt(pabs(y_hi), pabs(x_hi));
  Packet r_hi_1, r_lo_1;
  fast_twosum(x_hi, y_hi,r_hi_1, r_lo_1);
  Packet r_hi_2, r_lo_2;
  fast_twosum(y_hi, x_hi,r_hi_2, r_lo_2);
  const Packet r_hi = pselect(x_greater_mask, r_hi_1, r_hi_2);

  const Packet s1 = padd(padd(y_lo, r_lo_1), x_lo);
  const Packet s2 = padd(padd(x_lo, r_lo_2), y_lo);
  const Packet s = pselect(x_greater_mask, s1, s2);

  fast_twosum(r_hi, s, s_hi, s_lo);
}

template<typename Packet>
EIGEN_STRONG_INLINE
  void fast_twosum(const Packet& x_hi, const Packet& x_lo,
              const Packet& y_hi, const Packet& y_lo,
              Packet& s_hi, Packet& s_lo) {
  Packet r_hi, r_lo;
  fast_twosum(x_hi, y_hi, r_hi, r_lo);
  const Packet s = padd(padd(y_lo, r_lo), x_lo);
  fast_twosum(r_hi, s, s_hi, s_lo);
}

template<typename Packet>
EIGEN_STRONG_INLINE
void fast_twosum(const Packet& x,
                 const Packet& y_hi, const Packet& y_lo,
                 Packet& s_hi, Packet& s_lo) {
  Packet r_hi, r_lo;
  fast_twosum(x, y_hi, r_hi, r_lo);
  const Packet s = padd(y_lo, r_lo);
  fast_twosum(r_hi, s, s_hi, s_lo);
}

template<typename Packet>
EIGEN_STRONG_INLINE
void twoprod(const Packet& x_hi, const Packet& x_lo, const Packet& y,
             Packet& p_hi, Packet& p_lo) {
  Packet c_hi, c_lo1;
  twoprod(x_hi, y, c_hi, c_lo1);
  const Packet c_lo2 = pmul(x_lo, y);
  Packet t_hi, t_lo1;
  fast_twosum(c_hi, c_lo2, t_hi, t_lo1);
  const Packet t_lo2 = padd(t_lo1, c_lo1);
  fast_twosum(t_hi, t_lo2, p_hi, p_lo);
}

template<typename Packet>
EIGEN_STRONG_INLINE
void twoprod(const Packet& x_hi, const Packet& x_lo,
             const Packet& y_hi, const Packet& y_lo,
             Packet& p_hi, Packet& p_lo) {
  Packet p_hi_hi, p_hi_lo;
  twoprod(x_hi, x_lo, y_hi, p_hi_hi, p_hi_lo);
  Packet p_lo_hi, p_lo_lo;
  twoprod(x_hi, x_lo, y_lo, p_lo_hi, p_lo_lo);
  fast_twosum(p_hi_hi, p_hi_lo, p_lo_hi, p_lo_lo, p_hi, p_lo);
}

template <typename Packet>
void doubleword_reciprocal(const Packet& x, Packet& recip_hi, Packet& recip_lo) {
  typedef typename unpacket_traits<Packet>::type Scalar;
  
  Packet approx_recip = prsqrt(x);
  approx_recip = pmul(approx_recip, approx_recip);

  Packet t1_hi, t1_lo;
  twoprod(pnegate(x), approx_recip, t1_hi, t1_lo);
  
  Packet t2_hi, t2_lo;
  fast_twosum(pset1<Packet>(Scalar(2)), t1_hi, t2_hi, t2_lo);
  Packet t3_hi, t3_lo;
  fast_twosum(t2_hi, padd(t2_lo, t1_lo), t3_hi, t3_lo);
  
  twoprod(t3_hi, t3_lo, approx_recip, recip_hi, recip_lo);
}

template <typename Scalar>
struct accurate_log2 {
  template <typename Packet>
  EIGEN_STRONG_INLINE
  void operator()(const Packet& x, Packet& log2_x_hi, Packet& log2_x_lo) {
    log2_x_hi = plog2(x);
    log2_x_lo = pzero(x);
  }
};

template <>
struct accurate_log2<float> {
  template <typename Packet>
  EIGEN_STRONG_INLINE
  void operator()(const Packet& z, Packet& log2_x_hi, Packet& log2_x_lo) {

    const Packet p6 = pset1<Packet>( 9.703654795885e-2f);
    const Packet p5 = pset1<Packet>(-0.1690667718648f);
    const Packet p4 = pset1<Packet>( 0.1720575392246f);
    const Packet p3 = pset1<Packet>(-0.1789081543684f);
    const Packet p2 = pset1<Packet>( 0.2050433009862f);
    const Packet p1 = pset1<Packet>(-0.2404672354459f);
    const Packet p0 = pset1<Packet>( 0.2885761857032f);

    const Packet C3_hi = pset1<Packet>(-0.360674142838f);
    const Packet C3_lo = pset1<Packet>(-6.13283912543e-09f);
    const Packet C2_hi = pset1<Packet>(0.480897903442f);
    const Packet C2_lo = pset1<Packet>(-1.44861207474e-08f);
    const Packet C1_hi = pset1<Packet>(-0.721347510815f);
    const Packet C1_lo = pset1<Packet>(-4.84483164698e-09f);
    const Packet C0_hi = pset1<Packet>(1.44269502163f);
    const Packet C0_lo = pset1<Packet>(2.01711713999e-08f);
    const Packet one = pset1<Packet>(1.0f);

    const Packet x = psub(z, one);

    Packet x2 = pmul(x,x);
    Packet p_even = pmadd(p6, x2, p4);
    p_even = pmadd(p_even, x2, p2);
    p_even = pmadd(p_even, x2, p0);
    Packet p_odd = pmadd(p5, x2, p3);
    p_odd = pmadd(p_odd, x2, p1);
    Packet p = pmadd(p_odd, x, p_even);

    Packet q_hi, q_lo;
    Packet t_hi, t_lo;
    
    twoprod(p, x, t_hi, t_lo);
    fast_twosum(C3_hi, C3_lo, t_hi, t_lo, q_hi, q_lo);
    
    twoprod(q_hi, q_lo, x, t_hi, t_lo);
    fast_twosum(C2_hi, C2_lo, t_hi, t_lo, q_hi, q_lo);
    
    twoprod(q_hi, q_lo, x, t_hi, t_lo);
    fast_twosum(C1_hi, C1_lo, t_hi, t_lo, q_hi, q_lo);
    
    twoprod(q_hi, q_lo, x, t_hi, t_lo);
    fast_twosum(C0_hi, C0_lo, t_hi, t_lo, q_hi, q_lo);

    twoprod(q_hi, q_lo, x, log2_x_hi, log2_x_lo);
  }
};

template <>
struct accurate_log2<double> {
  template <typename Packet>
  EIGEN_STRONG_INLINE
  void operator()(const Packet& x, Packet& log2_x_hi, Packet& log2_x_lo) {

    const Packet q12 = pset1<Packet>(2.87074255468000586e-9);
    const Packet q10 = pset1<Packet>(2.38957980901884082e-8);
    const Packet q8 = pset1<Packet>(2.31032094540014656e-7);
    const Packet q6 = pset1<Packet>(2.27279857398537278e-6);
    const Packet q4 = pset1<Packet>(2.31271023278625638e-5);
    const Packet q2 = pset1<Packet>(2.47556738444535513e-4);
    const Packet q0 = pset1<Packet>(2.88543873228900172e-3);
    const Packet C_hi = pset1<Packet>(0.0400377511598501157);
    const Packet C_lo = pset1<Packet>(-4.77726582251425391e-19);
    const Packet one = pset1<Packet>(1.0);

    const Packet cst_2_log2e_hi = pset1<Packet>(2.88539008177792677);
    const Packet cst_2_log2e_lo = pset1<Packet>(4.07660016854549667e-17);
    
    Packet num_hi, num_lo;
    twoprod(cst_2_log2e_hi, cst_2_log2e_lo, psub(x, one), num_hi, num_lo);

    Packet denom_hi, denom_lo;
    doubleword_reciprocal(padd(x, one), denom_hi, denom_lo);
    
    Packet r_hi, r_lo;
    twoprod(num_hi, num_lo, denom_hi, denom_lo, r_hi, r_lo);
    
    Packet r2_hi, r2_lo;
    twoprod(r_hi, r_lo, r_hi, r_lo, r2_hi, r2_lo);
    
    Packet r4_hi, r4_lo;
    twoprod(r2_hi, r2_lo, r2_hi, r2_lo, r4_hi, r4_lo);

    Packet q_even = pmadd(q12, r4_hi, q8);
    Packet q_odd = pmadd(q10, r4_hi, q6);
    q_even = pmadd(q_even, r4_hi, q4);
    q_odd = pmadd(q_odd, r4_hi, q2);
    q_even = pmadd(q_even, r4_hi, q0);
    Packet q = pmadd(q_odd, r2_hi, q_even);

    Packet p_hi, p_lo;
    twoprod(r2_hi, r2_lo, q, p_hi, p_lo);
    
    Packet p1_hi, p1_lo;
    fast_twosum(C_hi, C_lo, p_hi, p_lo, p1_hi, p1_lo);
    
    Packet p2_hi, p2_lo;
    twoprod(r2_hi, r2_lo, p1_hi, p1_lo, p2_hi, p2_lo);
    
    Packet p3_hi, p3_lo;
    fast_twosum(one, p2_hi, p2_lo, p3_hi, p3_lo);

    twoprod(p3_hi, p3_lo, r_hi, r_lo, log2_x_hi, log2_x_lo);
  }
};

template <typename Scalar>
struct fast_accurate_exp2 {
  template <typename Packet>
  EIGEN_STRONG_INLINE
  Packet operator()(const Packet& x) {
    
    return pexp(pmul(pset1<Packet>(Scalar(EIGEN_LN2)), x));
  }
};

template <>
struct fast_accurate_exp2<float> {
  template <typename Packet>
  EIGEN_STRONG_INLINE
  Packet operator()(const Packet& x) {

    const Packet p4 = pset1<Packet>(1.539513905e-4f);
    const Packet p3 = pset1<Packet>(1.340007293e-3f);
    const Packet p2 = pset1<Packet>(9.618283249e-3f);
    const Packet p1 = pset1<Packet>(5.550328270e-2f);
    const Packet p0 = pset1<Packet>(0.2402264923f);

    const Packet C_hi = pset1<Packet>(0.6931471825f);
    const Packet C_lo = pset1<Packet>(2.36836577e-08f);
    const Packet one = pset1<Packet>(1.0f);

    Packet x2 = pmul(x,x);
    Packet p_even = pmadd(p4, x2, p2);
    Packet p_odd = pmadd(p3, x2, p1);
    p_even = pmadd(p_even, x2, p0);
    Packet p = pmadd(p_odd, x, p_even);

    Packet p_hi, p_lo;
    
    twoprod(p, x, p_hi, p_lo);
    
    Packet q1_hi, q1_lo;
    twosum(p_hi, p_lo, C_hi, C_lo, q1_hi, q1_lo);
    
    Packet q2_hi, q2_lo;
    twoprod(q1_hi, q1_lo, x, q2_hi, q2_lo);
    
    Packet q3_hi, q3_lo;

    fast_twosum(one, q2_hi, q3_hi, q3_lo);
    return padd(q3_hi, padd(q2_lo, q3_lo));
  }
};

template <>
struct fast_accurate_exp2<double> {
  template <typename Packet>
  EIGEN_STRONG_INLINE
  Packet operator()(const Packet& x) {

    const Packet p9 = pset1<Packet>(4.431642109085495276e-10);
    const Packet p8 = pset1<Packet>(7.073829923303358410e-9);
    const Packet p7 = pset1<Packet>(1.017822306737031311e-7);
    const Packet p6 = pset1<Packet>(1.321543498017646657e-6);
    const Packet p5 = pset1<Packet>(1.525273342728892877e-5);
    const Packet p4 = pset1<Packet>(1.540353045780084423e-4);
    const Packet p3 = pset1<Packet>(1.333355814685869807e-3);
    const Packet p2 = pset1<Packet>(9.618129107593478832e-3);
    const Packet p1 = pset1<Packet>(5.550410866481961247e-2);
    const Packet p0 = pset1<Packet>(0.240226506959101332);
    const Packet C_hi = pset1<Packet>(0.693147180559945286); 
    const Packet C_lo = pset1<Packet>(4.81927865669806721e-17);
    const Packet one = pset1<Packet>(1.0);

    Packet x2 = pmul(x,x);
    Packet p_even = pmadd(p8, x2, p6);
    Packet p_odd = pmadd(p9, x2, p7);
    p_even = pmadd(p_even, x2, p4);
    p_odd = pmadd(p_odd, x2, p5);
    p_even = pmadd(p_even, x2, p2);
    p_odd = pmadd(p_odd, x2, p3);
    p_even = pmadd(p_even, x2, p0);
    p_odd = pmadd(p_odd, x2, p1);
    Packet p = pmadd(p_odd, x, p_even);

    Packet p_hi, p_lo;
    
    twoprod(p, x, p_hi, p_lo);
    
    Packet q1_hi, q1_lo;
    twosum(p_hi, p_lo, C_hi, C_lo, q1_hi, q1_lo);
    
    Packet q2_hi, q2_lo;
    twoprod(q1_hi, q1_lo, x, q2_hi, q2_lo);
    
    Packet q3_hi, q3_lo;

    fast_twosum(one, q2_hi, q3_hi, q3_lo);
    return padd(q3_hi, padd(q2_lo, q3_lo));
  }
};

template <typename Packet>
EIGEN_STRONG_INLINE Packet generic_pow_impl(const Packet& x, const Packet& y) {
  typedef typename unpacket_traits<Packet>::type Scalar;
  
  Packet e_x;
  Packet m_x = pfrexp(x, e_x);

  EIGEN_CONSTEXPR Scalar sqrt_half = Scalar(0.70710678118654752440);
  const Packet m_x_scale_mask = pcmp_lt(m_x, pset1<Packet>(sqrt_half));
  m_x = pselect(m_x_scale_mask, pmul(pset1<Packet>(Scalar(2)), m_x), m_x);
  e_x = pselect(m_x_scale_mask, psub(e_x, pset1<Packet>(Scalar(1))), e_x);

  Packet rx_hi, rx_lo;
  accurate_log2<Scalar>()(m_x, rx_hi, rx_lo);

  Packet f1_hi, f1_lo, f2_hi, f2_lo;
  twoprod(e_x, y, f1_hi, f1_lo);
  twoprod(rx_hi, rx_lo, y, f2_hi, f2_lo);

  Packet f_hi, f_lo;
  fast_twosum(f1_hi, f1_lo, f2_hi, f2_lo, f_hi, f_lo);

  Packet n_z, r_z;
  absolute_split(f_hi, n_z, r_z);
  r_z = padd(r_z, f_lo);
  Packet n_r;
  absolute_split(r_z, n_r, r_z);
  n_z = padd(n_z, n_r);

  const Packet e_r = fast_accurate_exp2<Scalar>()(r_z);
  return pldexp(e_r, n_z);
}

template<typename Packet>
EIGEN_DEFINE_FUNCTION_ALLOWING_MULTIPLE_DEFINITIONS
EIGEN_UNUSED
Packet generic_pow(const Packet& x, const Packet& y) {
  typedef typename unpacket_traits<Packet>::type Scalar;

  const Packet cst_pos_inf = pset1<Packet>(NumTraits<Scalar>::infinity());
  const Packet cst_zero = pset1<Packet>(Scalar(0));
  const Packet cst_one = pset1<Packet>(Scalar(1));
  const Packet cst_nan = pset1<Packet>(NumTraits<Scalar>::quiet_NaN());

  const Packet abs_x = pabs(x);
  
  const Packet x_is_zero = pcmp_eq(x, cst_zero);
  const Packet x_is_neg = pcmp_lt(x, cst_zero);
  const Packet abs_x_is_inf = pcmp_eq(abs_x, cst_pos_inf);
  const Packet abs_x_is_one =  pcmp_eq(abs_x, cst_one);
  const Packet abs_x_is_gt_one = pcmp_lt(cst_one, abs_x);
  const Packet abs_x_is_lt_one = pcmp_lt(abs_x, cst_one);
  const Packet x_is_one =  pandnot(abs_x_is_one, x_is_neg);
  const Packet x_is_neg_one =  pand(abs_x_is_one, x_is_neg);
  const Packet x_is_nan = pandnot(ptrue(x), pcmp_eq(x, x));

  const Packet y_is_one = pcmp_eq(y, cst_one);
  const Packet y_is_zero = pcmp_eq(y, cst_zero);
  const Packet y_is_neg = pcmp_lt(y, cst_zero);
  const Packet y_is_pos = pandnot(ptrue(y), por(y_is_zero, y_is_neg));
  const Packet y_is_nan = pandnot(ptrue(y), pcmp_eq(y, y));
  const Packet abs_y_is_inf = pcmp_eq(pabs(y), cst_pos_inf);
  EIGEN_CONSTEXPR Scalar huge_exponent =
      (NumTraits<Scalar>::max_exponent() * Scalar(EIGEN_LN2)) /
       NumTraits<Scalar>::epsilon();
  const Packet abs_y_is_huge = pcmp_le(pset1<Packet>(huge_exponent), pabs(y));

  const Packet y_is_int = pcmp_eq(pfloor(y), y);
  const Packet y_div_2 = pmul(y, pset1<Packet>(Scalar(0.5)));
  const Packet y_is_even = pcmp_eq(pround(y_div_2), y_div_2);

  const Packet invalid_negative_x = pandnot(pandnot(pandnot(x_is_neg, abs_x_is_inf),
                                                    y_is_int),
                                            abs_y_is_inf);
  const Packet pow_is_one = por(por(x_is_one, y_is_zero),
                                pand(x_is_neg_one,
                                     por(abs_y_is_inf, pandnot(y_is_even, invalid_negative_x))));
  const Packet pow_is_nan = por(invalid_negative_x, por(x_is_nan, y_is_nan));
  const Packet pow_is_zero = por(por(por(pand(x_is_zero, y_is_pos),
                                         pand(abs_x_is_inf, y_is_neg)),
                                     pand(pand(abs_x_is_lt_one, abs_y_is_huge),
                                          y_is_pos)),
                                 pand(pand(abs_x_is_gt_one, abs_y_is_huge),
                                      y_is_neg));
  const Packet pow_is_inf = por(por(por(pand(x_is_zero, y_is_neg),
                                        pand(abs_x_is_inf, y_is_pos)),
                                    pand(pand(abs_x_is_lt_one, abs_y_is_huge),
                                         y_is_neg)),
                                pand(pand(abs_x_is_gt_one, abs_y_is_huge),
                                     y_is_pos));

  const Packet negate_pow_abs = pandnot(x_is_neg, y_is_even);
  const Packet pow_abs = generic_pow_impl(abs_x, y);
  return pselect(y_is_one, x,
                 pselect(pow_is_one, cst_one,
                         pselect(pow_is_nan, cst_nan,
                                 pselect(pow_is_inf, cst_pos_inf,
                                         pselect(pow_is_zero, cst_zero,
                                                 pselect(negate_pow_abs, pnegate(pow_abs), pow_abs))))));
}

template <typename Packet, int N>
struct ppolevl {
  static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet run(const Packet& x, const typename unpacket_traits<Packet>::type coeff[]) {
    EIGEN_STATIC_ASSERT((N > 0), YOU_MADE_A_PROGRAMMING_MISTAKE);
    return pmadd(ppolevl<Packet, N-1>::run(x, coeff), x, pset1<Packet>(coeff[N]));
  }
};

template <typename Packet>
struct ppolevl<Packet, 0> {
  static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet run(const Packet& x, const typename unpacket_traits<Packet>::type coeff[]) {
    EIGEN_UNUSED_VARIABLE(x);
    return pset1<Packet>(coeff[0]);
  }
};

template <typename Packet, int N>
struct pchebevl {
  EIGEN_DEVICE_FUNC
  static EIGEN_STRONG_INLINE Packet run(Packet x, const typename unpacket_traits<Packet>::type coef[]) {
    typedef typename unpacket_traits<Packet>::type Scalar;
    Packet b0 = pset1<Packet>(coef[0]);
    Packet b1 = pset1<Packet>(static_cast<Scalar>(0.f));
    Packet b2;

    for (int i = 1; i < N; i++) {
      b2 = b1;
      b1 = b0;
      b0 = psub(pmadd(x, b1, pset1<Packet>(coef[i])), b2);
    }

    return pmul(pset1<Packet>(static_cast<Scalar>(0.5f)), psub(b0, b2));
  }
};

} 
} 

#endif 
