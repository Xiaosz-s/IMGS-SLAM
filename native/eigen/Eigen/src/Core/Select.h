#ifndef EIGEN_SELECT_H
#define EIGEN_SELECT_H

namespace Eigen {

namespace internal {
template<typename ConditionMatrixType, typename ThenMatrixType, typename ElseMatrixType>
struct traits<Select<ConditionMatrixType, ThenMatrixType, ElseMatrixType> >
 : traits<ThenMatrixType>
{
  typedef typename traits<ThenMatrixType>::Scalar Scalar;
  typedef Dense StorageKind;
  typedef typename traits<ThenMatrixType>::XprKind XprKind;
  typedef typename ConditionMatrixType::Nested ConditionMatrixNested;
  typedef typename ThenMatrixType::Nested ThenMatrixNested;
  typedef typename ElseMatrixType::Nested ElseMatrixNested;
  enum {
    RowsAtCompileTime = ConditionMatrixType::RowsAtCompileTime,
    ColsAtCompileTime = ConditionMatrixType::ColsAtCompileTime,
    MaxRowsAtCompileTime = ConditionMatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = ConditionMatrixType::MaxColsAtCompileTime,
    Flags = (unsigned int)ThenMatrixType::Flags & ElseMatrixType::Flags & RowMajorBit
  };
};
}

template<typename ConditionMatrixType, typename ThenMatrixType, typename ElseMatrixType>
class Select : public internal::dense_xpr_base< Select<ConditionMatrixType, ThenMatrixType, ElseMatrixType> >::type,
               internal::no_assignment_operator
{
  public:

    typedef typename internal::dense_xpr_base<Select>::type Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(Select)

    inline EIGEN_DEVICE_FUNC
    Select(const ConditionMatrixType& a_conditionMatrix,
           const ThenMatrixType& a_thenMatrix,
           const ElseMatrixType& a_elseMatrix)
      : m_condition(a_conditionMatrix), m_then(a_thenMatrix), m_else(a_elseMatrix)
    {
      eigen_assert(m_condition.rows() == m_then.rows() && m_condition.rows() == m_else.rows());
      eigen_assert(m_condition.cols() == m_then.cols() && m_condition.cols() == m_else.cols());
    }

    inline EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    Index rows() const EIGEN_NOEXCEPT { return m_condition.rows(); }
    inline EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    Index cols() const EIGEN_NOEXCEPT { return m_condition.cols(); }

    inline EIGEN_DEVICE_FUNC
    const Scalar coeff(Index i, Index j) const
    {
      if (m_condition.coeff(i,j))
        return m_then.coeff(i,j);
      else
        return m_else.coeff(i,j);
    }

    inline EIGEN_DEVICE_FUNC
    const Scalar coeff(Index i) const
    {
      if (m_condition.coeff(i))
        return m_then.coeff(i);
      else
        return m_else.coeff(i);
    }

    inline EIGEN_DEVICE_FUNC const ConditionMatrixType& conditionMatrix() const
    {
      return m_condition;
    }

    inline EIGEN_DEVICE_FUNC const ThenMatrixType& thenMatrix() const
    {
      return m_then;
    }

    inline EIGEN_DEVICE_FUNC const ElseMatrixType& elseMatrix() const
    {
      return m_else;
    }

  protected:
    typename ConditionMatrixType::Nested m_condition;
    typename ThenMatrixType::Nested m_then;
    typename ElseMatrixType::Nested m_else;
};

template<typename Derived>
template<typename ThenDerived,typename ElseDerived>
inline EIGEN_DEVICE_FUNC const Select<Derived,ThenDerived,ElseDerived>
DenseBase<Derived>::select(const DenseBase<ThenDerived>& thenMatrix,
                            const DenseBase<ElseDerived>& elseMatrix) const
{
  return Select<Derived,ThenDerived,ElseDerived>(derived(), thenMatrix.derived(), elseMatrix.derived());
}

template<typename Derived>
template<typename ThenDerived>
inline EIGEN_DEVICE_FUNC const Select<Derived,ThenDerived, typename ThenDerived::ConstantReturnType>
DenseBase<Derived>::select(const DenseBase<ThenDerived>& thenMatrix,
                           const typename ThenDerived::Scalar& elseScalar) const
{
  return Select<Derived,ThenDerived,typename ThenDerived::ConstantReturnType>(
    derived(), thenMatrix.derived(), ThenDerived::Constant(rows(),cols(),elseScalar));
}

template<typename Derived>
template<typename ElseDerived>
inline EIGEN_DEVICE_FUNC const Select<Derived, typename ElseDerived::ConstantReturnType, ElseDerived >
DenseBase<Derived>::select(const typename ElseDerived::Scalar& thenScalar,
                           const DenseBase<ElseDerived>& elseMatrix) const
{
  return Select<Derived,typename ElseDerived::ConstantReturnType,ElseDerived>(
    derived(), ElseDerived::Constant(rows(),cols(),thenScalar), elseMatrix.derived());
}

} 

#endif 
