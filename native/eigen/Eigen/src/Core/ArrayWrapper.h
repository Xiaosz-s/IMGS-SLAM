#ifndef EIGEN_ARRAYWRAPPER_H
#define EIGEN_ARRAYWRAPPER_H

namespace Eigen {

namespace internal {
template<typename ExpressionType>
struct traits<ArrayWrapper<ExpressionType> >
  : public traits<typename remove_all<typename ExpressionType::Nested>::type >
{
  typedef ArrayXpr XprKind;
  
  enum {
    Flags0 = traits<typename remove_all<typename ExpressionType::Nested>::type >::Flags,
    LvalueBitFlag = is_lvalue<ExpressionType>::value ? LvalueBit : 0,
    Flags = (Flags0 & ~(NestByRefBit | LvalueBit)) | LvalueBitFlag
  };
};
}

template<typename ExpressionType>
class ArrayWrapper : public ArrayBase<ArrayWrapper<ExpressionType> >
{
  public:
    typedef ArrayBase<ArrayWrapper> Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(ArrayWrapper)
    EIGEN_INHERIT_ASSIGNMENT_OPERATORS(ArrayWrapper)
    typedef typename internal::remove_all<ExpressionType>::type NestedExpression;

    typedef typename internal::conditional<
                       internal::is_lvalue<ExpressionType>::value,
                       Scalar,
                       const Scalar
                     >::type ScalarWithConstIfNotLvalue;

    typedef typename internal::ref_selector<ExpressionType>::non_const_type NestedExpressionType;

    using Base::coeffRef;

    EIGEN_DEVICE_FUNC
    explicit EIGEN_STRONG_INLINE ArrayWrapper(ExpressionType& matrix) : m_expression(matrix) {}

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index rows() const EIGEN_NOEXCEPT { return m_expression.rows(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index cols() const EIGEN_NOEXCEPT { return m_expression.cols(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index outerStride() const EIGEN_NOEXCEPT { return m_expression.outerStride(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index innerStride() const EIGEN_NOEXCEPT { return m_expression.innerStride(); }

    EIGEN_DEVICE_FUNC
    inline ScalarWithConstIfNotLvalue* data() { return m_expression.data(); }
    EIGEN_DEVICE_FUNC
    inline const Scalar* data() const { return m_expression.data(); }

    EIGEN_DEVICE_FUNC
    inline const Scalar& coeffRef(Index rowId, Index colId) const
    {
      return m_expression.coeffRef(rowId, colId);
    }

    EIGEN_DEVICE_FUNC
    inline const Scalar& coeffRef(Index index) const
    {
      return m_expression.coeffRef(index);
    }

    template<typename Dest>
    EIGEN_DEVICE_FUNC
    inline void evalTo(Dest& dst) const { dst = m_expression; }

    EIGEN_DEVICE_FUNC
    const typename internal::remove_all<NestedExpressionType>::type&
    nestedExpression() const
    {
      return m_expression;
    }

    EIGEN_DEVICE_FUNC
    void resize(Index newSize) { m_expression.resize(newSize); }

    EIGEN_DEVICE_FUNC
    void resize(Index rows, Index cols) { m_expression.resize(rows,cols); }

  protected:
    NestedExpressionType m_expression;
};

namespace internal {
template<typename ExpressionType>
struct traits<MatrixWrapper<ExpressionType> >
 : public traits<typename remove_all<typename ExpressionType::Nested>::type >
{
  typedef MatrixXpr XprKind;
  
  enum {
    Flags0 = traits<typename remove_all<typename ExpressionType::Nested>::type >::Flags,
    LvalueBitFlag = is_lvalue<ExpressionType>::value ? LvalueBit : 0,
    Flags = (Flags0 & ~(NestByRefBit | LvalueBit)) | LvalueBitFlag
  };
};
}

template<typename ExpressionType>
class MatrixWrapper : public MatrixBase<MatrixWrapper<ExpressionType> >
{
  public:
    typedef MatrixBase<MatrixWrapper<ExpressionType> > Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(MatrixWrapper)
    EIGEN_INHERIT_ASSIGNMENT_OPERATORS(MatrixWrapper)
    typedef typename internal::remove_all<ExpressionType>::type NestedExpression;

    typedef typename internal::conditional<
                       internal::is_lvalue<ExpressionType>::value,
                       Scalar,
                       const Scalar
                     >::type ScalarWithConstIfNotLvalue;

    typedef typename internal::ref_selector<ExpressionType>::non_const_type NestedExpressionType;

    using Base::coeffRef;

    EIGEN_DEVICE_FUNC
    explicit inline MatrixWrapper(ExpressionType& matrix) : m_expression(matrix) {}

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index rows() const EIGEN_NOEXCEPT { return m_expression.rows(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index cols() const EIGEN_NOEXCEPT { return m_expression.cols(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index outerStride() const EIGEN_NOEXCEPT { return m_expression.outerStride(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index innerStride() const EIGEN_NOEXCEPT { return m_expression.innerStride(); }

    EIGEN_DEVICE_FUNC
    inline ScalarWithConstIfNotLvalue* data() { return m_expression.data(); }
    EIGEN_DEVICE_FUNC
    inline const Scalar* data() const { return m_expression.data(); }

    EIGEN_DEVICE_FUNC
    inline const Scalar& coeffRef(Index rowId, Index colId) const
    {
      return m_expression.derived().coeffRef(rowId, colId);
    }

    EIGEN_DEVICE_FUNC
    inline const Scalar& coeffRef(Index index) const
    {
      return m_expression.coeffRef(index);
    }

    EIGEN_DEVICE_FUNC
    const typename internal::remove_all<NestedExpressionType>::type&
    nestedExpression() const
    {
      return m_expression;
    }

    EIGEN_DEVICE_FUNC
    void resize(Index newSize) { m_expression.resize(newSize); }

    EIGEN_DEVICE_FUNC
    void resize(Index rows, Index cols) { m_expression.resize(rows,cols); }

  protected:
    NestedExpressionType m_expression;
};

} 

#endif 
