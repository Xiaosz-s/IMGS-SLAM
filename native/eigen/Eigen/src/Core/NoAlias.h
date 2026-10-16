#ifndef EIGEN_NOALIAS_H
#define EIGEN_NOALIAS_H

namespace Eigen {

template<typename ExpressionType, template <typename> class StorageBase>
class NoAlias
{
  public:
    typedef typename ExpressionType::Scalar Scalar;
    
    EIGEN_DEVICE_FUNC
    explicit NoAlias(ExpressionType& expression) : m_expression(expression) {}
    
    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE ExpressionType& operator=(const StorageBase<OtherDerived>& other)
    {
      call_assignment_no_alias(m_expression, other.derived(), internal::assign_op<Scalar,typename OtherDerived::Scalar>());
      return m_expression;
    }
    
    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE ExpressionType& operator+=(const StorageBase<OtherDerived>& other)
    {
      call_assignment_no_alias(m_expression, other.derived(), internal::add_assign_op<Scalar,typename OtherDerived::Scalar>());
      return m_expression;
    }
    
    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE ExpressionType& operator-=(const StorageBase<OtherDerived>& other)
    {
      call_assignment_no_alias(m_expression, other.derived(), internal::sub_assign_op<Scalar,typename OtherDerived::Scalar>());
      return m_expression;
    }

    EIGEN_DEVICE_FUNC
    ExpressionType& expression() const
    {
      return m_expression;
    }

  protected:
    ExpressionType& m_expression;
};

template<typename Derived>
NoAlias<Derived,MatrixBase> EIGEN_DEVICE_FUNC MatrixBase<Derived>::noalias()
{
  return NoAlias<Derived, Eigen::MatrixBase >(derived());
}

} 

#endif 
