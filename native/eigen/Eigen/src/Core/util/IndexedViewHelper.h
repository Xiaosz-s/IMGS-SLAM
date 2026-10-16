#ifndef EIGEN_INDEXED_VIEW_HELPER_H
#define EIGEN_INDEXED_VIEW_HELPER_H

namespace Eigen {

namespace internal {
struct symbolic_last_tag {};
}

static const symbolic::SymbolExpr<internal::symbolic_last_tag> last; 

#ifdef EIGEN_PARSED_BY_DOXYGEN
static const auto lastp1 = last+fix<1>;
#else

static const symbolic::AddExpr<symbolic::SymbolExpr<internal::symbolic_last_tag>,symbolic::ValueExpr<Eigen::internal::FixedInt<1> > > lastp1(last+fix<1>());
#endif

namespace internal {

inline Index eval_expr_given_size(Index x, Index )   { return x; }

template<int N>
FixedInt<N> eval_expr_given_size(FixedInt<N> x, Index )   { return x; }

template<typename Derived>
Index eval_expr_given_size(const symbolic::BaseExpr<Derived> &x, Index size)
{
  return x.derived().eval(last=size-1);
}

template<typename T, typename EnableIf = void> struct get_compile_time_incr {
  enum { value = UndefinedIncr };
};

template<typename T>
EIGEN_CONSTEXPR Index first(const T& x) EIGEN_NOEXCEPT { return x.first(); }

template<typename T,int XprSize,typename EnableIf=void>
struct IndexedViewCompatibleType {
  typedef T type;
};

template<typename T,typename Q>
const T& makeIndexedViewCompatible(const T& x, Index , Q) { return x; }

struct SingleRange {
  enum {
    SizeAtCompileTime = 1
  };
  SingleRange(Index val) : m_value(val) {}
  Index operator[](Index) const { return m_value; }
  static EIGEN_CONSTEXPR Index size() EIGEN_NOEXCEPT { return 1; }
  Index first() const EIGEN_NOEXCEPT { return m_value; }
  Index m_value;
};

template<> struct get_compile_time_incr<SingleRange> {
  enum { value = 1 }; 
};

template<typename T, int XprSize>
struct IndexedViewCompatibleType<T,XprSize,typename internal::enable_if<internal::is_integral<T>::value>::type> {

  typedef SingleRange type;
};

template<typename T, int XprSize>
struct IndexedViewCompatibleType<T, XprSize, typename enable_if<symbolic::is_symbolic<T>::value>::type> {
  typedef SingleRange type;
};

template<typename T>
typename enable_if<symbolic::is_symbolic<T>::value,SingleRange>::type
makeIndexedViewCompatible(const T& id, Index size, SpecializedType) {
  return eval_expr_given_size(id,size);
}

struct all_t { all_t() {} };

template<int XprSize>
struct AllRange {
  enum { SizeAtCompileTime = XprSize };
  AllRange(Index size = XprSize) : m_size(size) {}
  EIGEN_CONSTEXPR Index operator[](Index i) const EIGEN_NOEXCEPT { return i; }
  EIGEN_CONSTEXPR Index size() const EIGEN_NOEXCEPT { return m_size.value(); }
  EIGEN_CONSTEXPR Index first() const EIGEN_NOEXCEPT { return 0; }
  variable_if_dynamic<Index,XprSize> m_size;
};

template<int XprSize>
struct IndexedViewCompatibleType<all_t,XprSize> {
  typedef AllRange<XprSize> type;
};

template<typename XprSizeType>
inline AllRange<get_fixed_value<XprSizeType>::value> makeIndexedViewCompatible(all_t , XprSizeType size, SpecializedType) {
  return AllRange<get_fixed_value<XprSizeType>::value>(size);
}

template<int Size> struct get_compile_time_incr<AllRange<Size> > {
  enum { value = 1 };
};

} 

static const Eigen::internal::all_t all; 

namespace placeholders {
  typedef symbolic::SymbolExpr<internal::symbolic_last_tag> last_t;
  typedef symbolic::AddExpr<symbolic::SymbolExpr<internal::symbolic_last_tag>,symbolic::ValueExpr<Eigen::internal::FixedInt<1> > > end_t;
  typedef Eigen::internal::all_t all_t;

  EIGEN_DEPRECATED static const all_t  all  = Eigen::all;    
  EIGEN_DEPRECATED static const last_t last = Eigen::last;   
  EIGEN_DEPRECATED static const end_t  end  = Eigen::lastp1; 
}

} 

#endif 
