#ifndef EIGEN_INTEGRAL_CONSTANT_H
#define EIGEN_INTEGRAL_CONSTANT_H

namespace Eigen {

namespace internal {

template<int N> class FixedInt;
template<int N> class VariableAndFixedInt;

template<int N> class FixedInt
{
public:
  static const int value = N;
  EIGEN_CONSTEXPR operator int() const { return value; }
  FixedInt() {}
  FixedInt( VariableAndFixedInt<N> other) {
    #ifndef EIGEN_INTERNAL_DEBUGGING
    EIGEN_UNUSED_VARIABLE(other);
    #endif
    eigen_internal_assert(int(other)==N);
  }

  FixedInt<-N> operator-() const { return FixedInt<-N>(); }
  template<int M>
  FixedInt<N+M> operator+( FixedInt<M>) const { return FixedInt<N+M>(); }
  template<int M>
  FixedInt<N-M> operator-( FixedInt<M>) const { return FixedInt<N-M>(); }
  template<int M>
  FixedInt<N*M> operator*( FixedInt<M>) const { return FixedInt<N*M>(); }
  template<int M>
  FixedInt<N/M> operator/( FixedInt<M>) const { return FixedInt<N/M>(); }
  template<int M>
  FixedInt<N%M> operator%( FixedInt<M>) const { return FixedInt<N%M>(); }
  template<int M>
  FixedInt<N|M> operator|( FixedInt<M>) const { return FixedInt<N|M>(); }
  template<int M>
  FixedInt<N&M> operator&( FixedInt<M>) const { return FixedInt<N&M>(); }

#if EIGEN_HAS_CXX14_VARIABLE_TEMPLATES
  
  FixedInt operator() () const { return *this; }

  VariableAndFixedInt<N> operator() (int val) const { return VariableAndFixedInt<N>(val); }
#else
  FixedInt ( FixedInt<N> (*)() ) {}
#endif

#if EIGEN_HAS_CXX11
  FixedInt(std::integral_constant<int,N>) {}
#endif
};

template<int N> class VariableAndFixedInt
{
public:
  static const int value = N;
  operator int() const { return m_value; }
  VariableAndFixedInt(int val) { m_value = val; }
protected:
  int m_value;
};

template<typename T, int Default=Dynamic> struct get_fixed_value {
  static const int value = Default;
};

template<int N,int Default> struct get_fixed_value<FixedInt<N>,Default> {
  static const int value = N;
};

#if !EIGEN_HAS_CXX14_VARIABLE_TEMPLATES
template<int N,int Default> struct get_fixed_value<FixedInt<N> (*)(),Default> {
  static const int value = N;
};
#endif

template<int N,int Default> struct get_fixed_value<VariableAndFixedInt<N>,Default> {
  static const int value = N ;
};

template<typename T, int N, int Default>
struct get_fixed_value<variable_if_dynamic<T,N>,Default> {
  static const int value = N;
};

template<typename T> EIGEN_DEVICE_FUNC Index get_runtime_value(const T &x) { return x; }
#if !EIGEN_HAS_CXX14_VARIABLE_TEMPLATES
template<int N> EIGEN_DEVICE_FUNC Index get_runtime_value(FixedInt<N> (*)()) { return N; }
#endif

template<typename T, int DynamicKey=Dynamic, typename EnableIf=void> struct cleanup_index_type { typedef T type; };

template<typename T, int DynamicKey> struct cleanup_index_type<T,DynamicKey,typename internal::enable_if<internal::is_integral<T>::value>::type> { typedef Index type; };

#if !EIGEN_HAS_CXX14_VARIABLE_TEMPLATES

template<int N, int DynamicKey> struct cleanup_index_type<FixedInt<N> (*)(), DynamicKey> { typedef FixedInt<N> type; };
#endif

template<int N, int DynamicKey> struct cleanup_index_type<VariableAndFixedInt<N>, DynamicKey> { typedef FixedInt<N> type; };

template<int DynamicKey> struct cleanup_index_type<VariableAndFixedInt<DynamicKey>, DynamicKey> { typedef Index type; };

#if EIGEN_HAS_CXX11
template<int N, int DynamicKey> struct cleanup_index_type<std::integral_constant<int,N>, DynamicKey> { typedef FixedInt<N> type; };
#endif

} 

#ifndef EIGEN_PARSED_BY_DOXYGEN

#if EIGEN_HAS_CXX14_VARIABLE_TEMPLATES
template<int N>
static const internal::FixedInt<N> fix{};
#else
template<int N>
inline internal::FixedInt<N> fix() { return internal::FixedInt<N>(); }

template<int N,typename T>
inline internal::VariableAndFixedInt<N> fix(T val) { return internal::VariableAndFixedInt<N>(internal::convert_index<int>(val)); }
#endif

#else 

template<int N>
static const auto fix();

template<int N>
static const auto fix(int val);

#endif 

} 

#endif 
