#ifndef EIGEN_INVERSE_H
#define EIGEN_INVERSE_H

namespace Eigen {

template<typename XprType,typename StorageKind> class InverseImpl;

namespace internal {

template<typename XprType>
struct traits<Inverse<XprType> >
  : traits<typename XprType::PlainObject>
{
  typedef typename XprType::PlainObject PlainObject;
  typedef traits<PlainObject> BaseTraits;
  enum {
    Flags = BaseTraits::Flags & RowMajorBit
  };
};

} 

template<typename XprType>
class Inverse : public InverseImpl<XprType,typename internal::traits<XprType>::StorageKind>
{
public:
  typedef typename XprType::StorageIndex StorageIndex;
  typedef typename XprType::Scalar                            Scalar;
  typedef typename internal::ref_selector<XprType>::type      XprTypeNested;
  typedef typename internal::remove_all<XprTypeNested>::type  XprTypeNestedCleaned;
  typedef typename internal::ref_selector<Inverse>::type Nested;
  typedef typename internal::remove_all<XprType>::type NestedExpression;

  explicit EIGEN_DEVICE_FUNC Inverse(const XprType &xpr)
    : m_xpr(xpr)
  {}

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR  Index rows() const EIGEN_NOEXCEPT { return m_xpr.cols(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR  Index cols() const EIGEN_NOEXCEPT { return m_xpr.rows(); }

  EIGEN_DEVICE_FUNC const XprTypeNestedCleaned& nestedExpression() const { return m_xpr; }

protected:
  XprTypeNested m_xpr;
};

template<typename XprType, typename StorageKind>
class InverseImpl
  : public internal::generic_xpr_base<Inverse<XprType> >::type
{
public:
  typedef typename internal::generic_xpr_base<Inverse<XprType> >::type Base;
  typedef typename XprType::Scalar Scalar;
private:

  Scalar coeff(Index row, Index col) const;
  Scalar coeff(Index i) const;
};

namespace internal {

template<typename ArgType>
struct unary_evaluator<Inverse<ArgType> >
  : public evaluator<typename Inverse<ArgType>::PlainObject>
{
  typedef Inverse<ArgType> InverseType;
  typedef typename InverseType::PlainObject PlainObject;
  typedef evaluator<PlainObject> Base;

  enum { Flags = Base::Flags | EvalBeforeNestingBit };

  unary_evaluator(const InverseType& inv_xpr)
    : m_result(inv_xpr.rows(), inv_xpr.cols())
  {
    ::new (static_cast<Base*>(this)) Base(m_result);
    internal::call_assignment_no_alias(m_result, inv_xpr);
  }

protected:
  PlainObject m_result;
};

} 

} 

#endif 
