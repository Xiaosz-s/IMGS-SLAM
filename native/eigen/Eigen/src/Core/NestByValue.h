#ifndef EIGEN_NESTBYVALUE_H
#define EIGEN_NESTBYVALUE_H

namespace Eigen {

namespace internal {
template<typename ExpressionType>
struct traits<NestByValue<ExpressionType> > : public traits<ExpressionType>
{
  enum {
    Flags = traits<ExpressionType>::Flags & ~NestByRefBit
  };
};
}

template<typename ExpressionType> class NestByValue
  : public internal::dense_xpr_base< NestByValue<ExpressionType> >::type
{
  public:

    typedef typename internal::dense_xpr_base<NestByValue>::type Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(NestByValue)

    EIGEN_DEVICE_FUNC explicit inline NestByValue(const ExpressionType& matrix) : m_expression(matrix) {}

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR inline Index rows() const EIGEN_NOEXCEPT { return m_expression.rows(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR inline Index cols() const EIGEN_NOEXCEPT { return m_expression.cols(); }

    EIGEN_DEVICE_FUNC operator const ExpressionType&() const { return m_expression; }

    EIGEN_DEVICE_FUNC const ExpressionType& nestedExpression() const { return m_expression; }

  protected:
    const ExpressionType m_expression;
};

template<typename Derived>
EIGEN_DEVICE_FUNC inline const NestByValue<Derived>
DenseBase<Derived>::nestByValue() const
{
  return NestByValue<Derived>(derived());
}

namespace internal {

template<typename ArgType>
struct evaluator<NestByValue<ArgType> >
  : public evaluator<ArgType>
{
  typedef evaluator<ArgType> Base;

  EIGEN_DEVICE_FUNC explicit evaluator(const NestByValue<ArgType>& xpr)
    : Base(xpr.nestedExpression())
  {}
};
}

} 

#endif 
