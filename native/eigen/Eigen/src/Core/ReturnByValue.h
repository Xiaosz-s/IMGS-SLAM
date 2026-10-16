#ifndef EIGEN_RETURNBYVALUE_H
#define EIGEN_RETURNBYVALUE_H

namespace Eigen {

namespace internal {

template<typename Derived>
struct traits<ReturnByValue<Derived> >
  : public traits<typename traits<Derived>::ReturnType>
{
  enum {

    Flags = (traits<typename traits<Derived>::ReturnType>::Flags
             | EvalBeforeNestingBit) & ~DirectAccessBit
  };
};

template<typename Derived,int n,typename PlainObject>
struct nested_eval<ReturnByValue<Derived>, n, PlainObject>
{
  typedef typename traits<Derived>::ReturnType type;
};

} 

template<typename Derived> class ReturnByValue
  : public internal::dense_xpr_base< ReturnByValue<Derived> >::type, internal::no_assignment_operator
{
  public:
    typedef typename internal::traits<Derived>::ReturnType ReturnType;

    typedef typename internal::dense_xpr_base<ReturnByValue>::type Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(ReturnByValue)

    template<typename Dest>
    EIGEN_DEVICE_FUNC
    inline void evalTo(Dest& dst) const
    { static_cast<const Derived*>(this)->evalTo(dst); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index rows() const EIGEN_NOEXCEPT { return static_cast<const Derived*>(this)->rows(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index cols() const EIGEN_NOEXCEPT { return static_cast<const Derived*>(this)->cols(); }

#ifndef EIGEN_PARSED_BY_DOXYGEN
#define Unusable YOU_ARE_TRYING_TO_ACCESS_A_SINGLE_COEFFICIENT_IN_A_SPECIAL_EXPRESSION_WHERE_THAT_IS_NOT_ALLOWED_BECAUSE_THAT_WOULD_BE_INEFFICIENT
    class Unusable{
      Unusable(const Unusable&) {}
      Unusable& operator=(const Unusable&) {return *this;}
    };
    const Unusable& coeff(Index) const { return *reinterpret_cast<const Unusable*>(this); }
    const Unusable& coeff(Index,Index) const { return *reinterpret_cast<const Unusable*>(this); }
    Unusable& coeffRef(Index) { return *reinterpret_cast<Unusable*>(this); }
    Unusable& coeffRef(Index,Index) { return *reinterpret_cast<Unusable*>(this); }
#undef Unusable
#endif
};

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC Derived& DenseBase<Derived>::operator=(const ReturnByValue<OtherDerived>& other)
{
  other.evalTo(derived());
  return derived();
}

namespace internal {

template<typename Derived>
struct evaluator<ReturnByValue<Derived> >
  : public evaluator<typename internal::traits<Derived>::ReturnType>
{
  typedef ReturnByValue<Derived> XprType;
  typedef typename internal::traits<Derived>::ReturnType PlainObject;
  typedef evaluator<PlainObject> Base;

  EIGEN_DEVICE_FUNC explicit evaluator(const XprType& xpr)
    : m_result(xpr.rows(), xpr.cols())
  {
    ::new (static_cast<Base*>(this)) Base(m_result);
    xpr.evalTo(m_result);
  }

protected:
  PlainObject m_result;
};

} 

} 

#endif 
