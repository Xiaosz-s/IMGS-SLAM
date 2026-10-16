#ifndef EIGEN_SYMBOLIC_INDEX_H
#define EIGEN_SYMBOLIC_INDEX_H

namespace Eigen {

namespace symbolic {

template<typename Tag> class Symbol;
template<typename Arg0> class NegateExpr;
template<typename Arg1,typename Arg2> class AddExpr;
template<typename Arg1,typename Arg2> class ProductExpr;
template<typename Arg1,typename Arg2> class QuotientExpr;

template<typename IndexType=Index>
class ValueExpr {
public:
  ValueExpr(IndexType val) : m_value(val) {}
  template<typename T>
  IndexType eval_impl(const T&) const { return m_value; }
protected:
  IndexType m_value;
};

template<int N>
class ValueExpr<internal::FixedInt<N> > {
public:
  ValueExpr() {}
  template<typename T>
  EIGEN_CONSTEXPR Index eval_impl(const T&) const { return N; }
};

template<typename Derived>
class BaseExpr
{
public:
  const Derived& derived() const { return *static_cast<const Derived*>(this); }

  template<typename T>
  Index eval(const T& values) const { return derived().eval_impl(values); }

#if EIGEN_HAS_CXX14
  template<typename... Types>
  Index eval(Types&&... values) const { return derived().eval_impl(std::make_tuple(values...)); }
#endif

  NegateExpr<Derived> operator-() const { return NegateExpr<Derived>(derived()); }

  AddExpr<Derived,ValueExpr<> > operator+(Index b) const
  { return AddExpr<Derived,ValueExpr<> >(derived(),  b); }
  AddExpr<Derived,ValueExpr<> > operator-(Index a) const
  { return AddExpr<Derived,ValueExpr<> >(derived(), -a); }
  ProductExpr<Derived,ValueExpr<> > operator*(Index a) const
  { return ProductExpr<Derived,ValueExpr<> >(derived(),a); }
  QuotientExpr<Derived,ValueExpr<> > operator/(Index a) const
  { return QuotientExpr<Derived,ValueExpr<> >(derived(),a); }

  friend AddExpr<Derived,ValueExpr<> > operator+(Index a, const BaseExpr& b)
  { return AddExpr<Derived,ValueExpr<> >(b.derived(), a); }
  friend AddExpr<NegateExpr<Derived>,ValueExpr<> > operator-(Index a, const BaseExpr& b)
  { return AddExpr<NegateExpr<Derived>,ValueExpr<> >(-b.derived(), a); }
  friend ProductExpr<ValueExpr<>,Derived> operator*(Index a, const BaseExpr& b)
  { return ProductExpr<ValueExpr<>,Derived>(a,b.derived()); }
  friend QuotientExpr<ValueExpr<>,Derived> operator/(Index a, const BaseExpr& b)
  { return QuotientExpr<ValueExpr<>,Derived>(a,b.derived()); }

  template<int N>
  AddExpr<Derived,ValueExpr<internal::FixedInt<N> > > operator+(internal::FixedInt<N>) const
  { return AddExpr<Derived,ValueExpr<internal::FixedInt<N> > >(derived(), ValueExpr<internal::FixedInt<N> >()); }
  template<int N>
  AddExpr<Derived,ValueExpr<internal::FixedInt<-N> > > operator-(internal::FixedInt<N>) const
  { return AddExpr<Derived,ValueExpr<internal::FixedInt<-N> > >(derived(), ValueExpr<internal::FixedInt<-N> >()); }
  template<int N>
  ProductExpr<Derived,ValueExpr<internal::FixedInt<N> > > operator*(internal::FixedInt<N>) const
  { return ProductExpr<Derived,ValueExpr<internal::FixedInt<N> > >(derived(),ValueExpr<internal::FixedInt<N> >()); }
  template<int N>
  QuotientExpr<Derived,ValueExpr<internal::FixedInt<N> > > operator/(internal::FixedInt<N>) const
  { return QuotientExpr<Derived,ValueExpr<internal::FixedInt<N> > >(derived(),ValueExpr<internal::FixedInt<N> >()); }

  template<int N>
  friend AddExpr<Derived,ValueExpr<internal::FixedInt<N> > > operator+(internal::FixedInt<N>, const BaseExpr& b)
  { return AddExpr<Derived,ValueExpr<internal::FixedInt<N> > >(b.derived(), ValueExpr<internal::FixedInt<N> >()); }
  template<int N>
  friend AddExpr<NegateExpr<Derived>,ValueExpr<internal::FixedInt<N> > > operator-(internal::FixedInt<N>, const BaseExpr& b)
  { return AddExpr<NegateExpr<Derived>,ValueExpr<internal::FixedInt<N> > >(-b.derived(), ValueExpr<internal::FixedInt<N> >()); }
  template<int N>
  friend ProductExpr<ValueExpr<internal::FixedInt<N> >,Derived> operator*(internal::FixedInt<N>, const BaseExpr& b)
  { return ProductExpr<ValueExpr<internal::FixedInt<N> >,Derived>(ValueExpr<internal::FixedInt<N> >(),b.derived()); }
  template<int N>
  friend QuotientExpr<ValueExpr<internal::FixedInt<N> >,Derived> operator/(internal::FixedInt<N>, const BaseExpr& b)
  { return QuotientExpr<ValueExpr<internal::FixedInt<N> > ,Derived>(ValueExpr<internal::FixedInt<N> >(),b.derived()); }

#if (!EIGEN_HAS_CXX14)
  template<int N>
  AddExpr<Derived,ValueExpr<internal::FixedInt<N> > > operator+(internal::FixedInt<N> (*)()) const
  { return AddExpr<Derived,ValueExpr<internal::FixedInt<N> > >(derived(), ValueExpr<internal::FixedInt<N> >()); }
  template<int N>
  AddExpr<Derived,ValueExpr<internal::FixedInt<-N> > > operator-(internal::FixedInt<N> (*)()) const
  { return AddExpr<Derived,ValueExpr<internal::FixedInt<-N> > >(derived(), ValueExpr<internal::FixedInt<-N> >()); }
  template<int N>
  ProductExpr<Derived,ValueExpr<internal::FixedInt<N> > > operator*(internal::FixedInt<N> (*)()) const
  { return ProductExpr<Derived,ValueExpr<internal::FixedInt<N> > >(derived(),ValueExpr<internal::FixedInt<N> >()); }
  template<int N>
  QuotientExpr<Derived,ValueExpr<internal::FixedInt<N> > > operator/(internal::FixedInt<N> (*)()) const
  { return QuotientExpr<Derived,ValueExpr<internal::FixedInt<N> > >(derived(),ValueExpr<internal::FixedInt<N> >()); }

  template<int N>
  friend AddExpr<Derived,ValueExpr<internal::FixedInt<N> > > operator+(internal::FixedInt<N> (*)(), const BaseExpr& b)
  { return AddExpr<Derived,ValueExpr<internal::FixedInt<N> > >(b.derived(), ValueExpr<internal::FixedInt<N> >()); }
  template<int N>
  friend AddExpr<NegateExpr<Derived>,ValueExpr<internal::FixedInt<N> > > operator-(internal::FixedInt<N> (*)(), const BaseExpr& b)
  { return AddExpr<NegateExpr<Derived>,ValueExpr<internal::FixedInt<N> > >(-b.derived(), ValueExpr<internal::FixedInt<N> >()); }
  template<int N>
  friend ProductExpr<ValueExpr<internal::FixedInt<N> >,Derived> operator*(internal::FixedInt<N> (*)(), const BaseExpr& b)
  { return ProductExpr<ValueExpr<internal::FixedInt<N> >,Derived>(ValueExpr<internal::FixedInt<N> >(),b.derived()); }
  template<int N>
  friend QuotientExpr<ValueExpr<internal::FixedInt<N> >,Derived> operator/(internal::FixedInt<N> (*)(), const BaseExpr& b)
  { return QuotientExpr<ValueExpr<internal::FixedInt<N> > ,Derived>(ValueExpr<internal::FixedInt<N> >(),b.derived()); }
#endif

  template<typename OtherDerived>
  AddExpr<Derived,OtherDerived> operator+(const BaseExpr<OtherDerived> &b) const
  { return AddExpr<Derived,OtherDerived>(derived(),  b.derived()); }

  template<typename OtherDerived>
  AddExpr<Derived,NegateExpr<OtherDerived> > operator-(const BaseExpr<OtherDerived> &b) const
  { return AddExpr<Derived,NegateExpr<OtherDerived> >(derived(), -b.derived()); }

  template<typename OtherDerived>
  ProductExpr<Derived,OtherDerived> operator*(const BaseExpr<OtherDerived> &b) const
  { return ProductExpr<Derived,OtherDerived>(derived(), b.derived()); }

  template<typename OtherDerived>
  QuotientExpr<Derived,OtherDerived> operator/(const BaseExpr<OtherDerived> &b) const
  { return QuotientExpr<Derived,OtherDerived>(derived(), b.derived()); }
};

template<typename T>
struct is_symbolic {
  
  enum { value = internal::is_convertible<T,BaseExpr<T> >::value };
};

template<typename Tag>
class SymbolValue
{
public:
  
  SymbolValue(Index val) : m_value(val) {}

  Index value() const { return m_value; }
protected:
  Index m_value;
};

template<typename tag>
class SymbolExpr : public BaseExpr<SymbolExpr<tag> >
{
public:
  
  typedef tag Tag;

  SymbolExpr() {}

  SymbolValue<Tag> operator=(Index val) const {
    return SymbolValue<Tag>(val);
  }

  Index eval_impl(const SymbolValue<Tag> &values) const { return values.value(); }

#if EIGEN_HAS_CXX14
  
  template<typename... Types>
  Index eval_impl(const std::tuple<Types...>& values) const { return std::get<SymbolValue<Tag> >(values).value(); }
#endif
};

template<typename Arg0>
class NegateExpr : public BaseExpr<NegateExpr<Arg0> >
{
public:
  NegateExpr(const Arg0& arg0) : m_arg0(arg0) {}

  template<typename T>
  Index eval_impl(const T& values) const { return -m_arg0.eval_impl(values); }
protected:
  Arg0 m_arg0;
};

template<typename Arg0, typename Arg1>
class AddExpr : public BaseExpr<AddExpr<Arg0,Arg1> >
{
public:
  AddExpr(const Arg0& arg0, const Arg1& arg1) : m_arg0(arg0), m_arg1(arg1) {}

  template<typename T>
  Index eval_impl(const T& values) const { return m_arg0.eval_impl(values) + m_arg1.eval_impl(values); }
protected:
  Arg0 m_arg0;
  Arg1 m_arg1;
};

template<typename Arg0, typename Arg1>
class ProductExpr : public BaseExpr<ProductExpr<Arg0,Arg1> >
{
public:
  ProductExpr(const Arg0& arg0, const Arg1& arg1) : m_arg0(arg0), m_arg1(arg1) {}

  template<typename T>
  Index eval_impl(const T& values) const { return m_arg0.eval_impl(values) * m_arg1.eval_impl(values); }
protected:
  Arg0 m_arg0;
  Arg1 m_arg1;
};

template<typename Arg0, typename Arg1>
class QuotientExpr : public BaseExpr<QuotientExpr<Arg0,Arg1> >
{
public:
  QuotientExpr(const Arg0& arg0, const Arg1& arg1) : m_arg0(arg0), m_arg1(arg1) {}

  template<typename T>
  Index eval_impl(const T& values) const { return m_arg0.eval_impl(values) / m_arg1.eval_impl(values); }
protected:
  Arg0 m_arg0;
  Arg1 m_arg1;
};

} 

} 

#endif 
