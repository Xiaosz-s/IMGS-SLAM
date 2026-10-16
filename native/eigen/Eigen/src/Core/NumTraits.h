#ifndef EIGEN_NUMTRAITS_H
#define EIGEN_NUMTRAITS_H

namespace Eigen {

namespace internal {

template< typename T,
          bool use_numeric_limits = std::numeric_limits<T>::is_specialized,
          bool is_integer = NumTraits<T>::IsInteger>
struct default_digits10_impl
{
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static int run() { return std::numeric_limits<T>::digits10; }
};

template<typename T>
struct default_digits10_impl<T,false,false> 
{
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static int run() {
    using std::log10;
    using std::ceil;
    typedef typename NumTraits<T>::Real Real;
    return int(ceil(-log10(NumTraits<Real>::epsilon())));
  }
};

template<typename T>
struct default_digits10_impl<T,false,true> 
{
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static int run() { return 0; }
};

template< typename T,
          bool use_numeric_limits = std::numeric_limits<T>::is_specialized,
          bool is_integer = NumTraits<T>::IsInteger>
struct default_digits_impl
{
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static int run() { return std::numeric_limits<T>::digits; }
};

template<typename T>
struct default_digits_impl<T,false,false> 
{
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static int run() {
    using std::log;
    using std::ceil;
    typedef typename NumTraits<T>::Real Real;
    return int(ceil(-log(NumTraits<Real>::epsilon())/log(static_cast<Real>(2))));
  }
};

template<typename T>
struct default_digits_impl<T,false,true> 
{
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static int run() { return 0; }
};

} 

namespace numext {

template <typename Tgt, typename Src>
EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC Tgt bit_cast(const Src& src) {
#if EIGEN_HAS_TYPE_TRAITS
  
  EIGEN_STATIC_ASSERT(std::is_trivially_copyable<Src>::value, THIS_TYPE_IS_NOT_SUPPORTED);
  EIGEN_STATIC_ASSERT(std::is_trivially_copyable<Tgt>::value && std::is_default_constructible<Tgt>::value,
                      THIS_TYPE_IS_NOT_SUPPORTED);
#endif

  EIGEN_STATIC_ASSERT(sizeof(Src) == sizeof(Tgt), THIS_TYPE_IS_NOT_SUPPORTED);
  Tgt tgt;
  EIGEN_USING_STD(memcpy)
  memcpy(&tgt, &src, sizeof(Tgt));
  return tgt;
}
}  

template<typename T> struct GenericNumTraits
{
  enum {
    IsInteger = std::numeric_limits<T>::is_integer,
    IsSigned = std::numeric_limits<T>::is_signed,
    IsComplex = 0,
    RequireInitialization = internal::is_arithmetic<T>::value ? 0 : 1,
    ReadCost = 1,
    AddCost = 1,
    MulCost = 1
  };

  typedef T Real;
  typedef typename internal::conditional<
                     IsInteger,
                     typename internal::conditional<sizeof(T)<=2, float, double>::type,
                     T
                   >::type NonInteger;
  typedef T Nested;
  typedef T Literal;

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline Real epsilon()
  {
    return numext::numeric_limits<T>::epsilon();
  }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline int digits10()
  {
    return internal::default_digits10_impl<T>::run();
  }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline int digits()
  {
    return internal::default_digits_impl<T>::run();
  }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline int min_exponent()
  {
    return numext::numeric_limits<T>::min_exponent;
  }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline int max_exponent()
  {
    return numext::numeric_limits<T>::max_exponent;
  }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline Real dummy_precision()
  {
    
    return Real(0);
  }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline T highest() {
    return (numext::numeric_limits<T>::max)();
  }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline T lowest()  {
    return IsInteger ? (numext::numeric_limits<T>::min)()
                     : static_cast<T>(-(numext::numeric_limits<T>::max)());
  }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline T infinity() {
    return numext::numeric_limits<T>::infinity();
  }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline T quiet_NaN() {
    return numext::numeric_limits<T>::quiet_NaN();
  }
};

template<typename T> struct NumTraits : GenericNumTraits<T>
{};

template<> struct NumTraits<float>
  : GenericNumTraits<float>
{
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline float dummy_precision() { return 1e-5f; }
};

template<> struct NumTraits<double> : GenericNumTraits<double>
{
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline double dummy_precision() { return 1e-12; }
};

template<> struct NumTraits<long double>
  : GenericNumTraits<long double>
{
  EIGEN_CONSTEXPR
  static inline long double dummy_precision() { return 1e-15l; }
};

template<typename Real_> struct NumTraits<std::complex<Real_> >
  : GenericNumTraits<std::complex<Real_> >
{
  typedef Real_ Real;
  typedef typename NumTraits<Real_>::Literal Literal;
  enum {
    IsComplex = 1,
    RequireInitialization = NumTraits<Real_>::RequireInitialization,
    ReadCost = 2 * NumTraits<Real_>::ReadCost,
    AddCost = 2 * NumTraits<Real>::AddCost,
    MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
  };

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline Real epsilon() { return NumTraits<Real>::epsilon(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline Real dummy_precision() { return NumTraits<Real>::dummy_precision(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline int digits10() { return NumTraits<Real>::digits10(); }
};

template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct NumTraits<Array<Scalar, Rows, Cols, Options, MaxRows, MaxCols> >
{
  typedef Array<Scalar, Rows, Cols, Options, MaxRows, MaxCols> ArrayType;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Array<RealScalar, Rows, Cols, Options, MaxRows, MaxCols> Real;
  typedef typename NumTraits<Scalar>::NonInteger NonIntegerScalar;
  typedef Array<NonIntegerScalar, Rows, Cols, Options, MaxRows, MaxCols> NonInteger;
  typedef ArrayType & Nested;
  typedef typename NumTraits<Scalar>::Literal Literal;

  enum {
    IsComplex = NumTraits<Scalar>::IsComplex,
    IsInteger = NumTraits<Scalar>::IsInteger,
    IsSigned  = NumTraits<Scalar>::IsSigned,
    RequireInitialization = 1,
    ReadCost = ArrayType::SizeAtCompileTime==Dynamic ? HugeCost : ArrayType::SizeAtCompileTime * int(NumTraits<Scalar>::ReadCost),
    AddCost  = ArrayType::SizeAtCompileTime==Dynamic ? HugeCost : ArrayType::SizeAtCompileTime * int(NumTraits<Scalar>::AddCost),
    MulCost  = ArrayType::SizeAtCompileTime==Dynamic ? HugeCost : ArrayType::SizeAtCompileTime * int(NumTraits<Scalar>::MulCost)
  };

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline RealScalar epsilon() { return NumTraits<RealScalar>::epsilon(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline RealScalar dummy_precision() { return NumTraits<RealScalar>::dummy_precision(); }

  EIGEN_CONSTEXPR
  static inline int digits10() { return NumTraits<Scalar>::digits10(); }
};

template<> struct NumTraits<std::string>
  : GenericNumTraits<std::string>
{
  enum {
    RequireInitialization = 1,
    ReadCost = HugeCost,
    AddCost  = HugeCost,
    MulCost  = HugeCost
  };

  EIGEN_CONSTEXPR
  static inline int digits10() { return 0; }

private:
  static inline std::string epsilon();
  static inline std::string dummy_precision();
  static inline std::string lowest();
  static inline std::string highest();
  static inline std::string infinity();
  static inline std::string quiet_NaN();
};

template<> struct NumTraits<void> {};

template<> struct NumTraits<bool> : GenericNumTraits<bool> {};

} 

#endif 
