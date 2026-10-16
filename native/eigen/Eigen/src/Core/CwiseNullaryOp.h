#ifndef EIGEN_CWISE_NULLARY_OP_H
#define EIGEN_CWISE_NULLARY_OP_H

namespace Eigen {

namespace internal {
template<typename NullaryOp, typename PlainObjectType>
struct traits<CwiseNullaryOp<NullaryOp, PlainObjectType> > : traits<PlainObjectType>
{
  enum {
    Flags = traits<PlainObjectType>::Flags & RowMajorBit
  };
};

} 

template<typename NullaryOp, typename PlainObjectType>
class CwiseNullaryOp : public internal::dense_xpr_base< CwiseNullaryOp<NullaryOp, PlainObjectType> >::type, internal::no_assignment_operator
{
  public:

    typedef typename internal::dense_xpr_base<CwiseNullaryOp>::type Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(CwiseNullaryOp)

    EIGEN_DEVICE_FUNC
    CwiseNullaryOp(Index rows, Index cols, const NullaryOp& func = NullaryOp())
      : m_rows(rows), m_cols(cols), m_functor(func)
    {
      eigen_assert(rows >= 0
            && (RowsAtCompileTime == Dynamic || RowsAtCompileTime == rows)
            &&  cols >= 0
            && (ColsAtCompileTime == Dynamic || ColsAtCompileTime == cols));
    }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index rows() const { return m_rows.value(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index cols() const { return m_cols.value(); }

    EIGEN_DEVICE_FUNC
    const NullaryOp& functor() const { return m_functor; }

  protected:
    const internal::variable_if_dynamic<Index, RowsAtCompileTime> m_rows;
    const internal::variable_if_dynamic<Index, ColsAtCompileTime> m_cols;
    const NullaryOp m_functor;
};

template<typename Derived>
template<typename CustomNullaryOp>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const CwiseNullaryOp<CustomNullaryOp,typename DenseBase<Derived>::PlainObject>
#else
const CwiseNullaryOp<CustomNullaryOp,PlainObject>
#endif
DenseBase<Derived>::NullaryExpr(Index rows, Index cols, const CustomNullaryOp& func)
{
  return CwiseNullaryOp<CustomNullaryOp, PlainObject>(rows, cols, func);
}

template<typename Derived>
template<typename CustomNullaryOp>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const CwiseNullaryOp<CustomNullaryOp, typename DenseBase<Derived>::PlainObject>
#else
const CwiseNullaryOp<CustomNullaryOp, PlainObject>
#endif
DenseBase<Derived>::NullaryExpr(Index size, const CustomNullaryOp& func)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  if(RowsAtCompileTime == 1) return CwiseNullaryOp<CustomNullaryOp, PlainObject>(1, size, func);
  else return CwiseNullaryOp<CustomNullaryOp, PlainObject>(size, 1, func);
}

template<typename Derived>
template<typename CustomNullaryOp>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const CwiseNullaryOp<CustomNullaryOp, typename DenseBase<Derived>::PlainObject>
#else
const CwiseNullaryOp<CustomNullaryOp, PlainObject>
#endif
DenseBase<Derived>::NullaryExpr(const CustomNullaryOp& func)
{
  return CwiseNullaryOp<CustomNullaryOp, PlainObject>(RowsAtCompileTime, ColsAtCompileTime, func);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::ConstantReturnType
DenseBase<Derived>::Constant(Index rows, Index cols, const Scalar& value)
{
  return DenseBase<Derived>::NullaryExpr(rows, cols, internal::scalar_constant_op<Scalar>(value));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::ConstantReturnType
DenseBase<Derived>::Constant(Index size, const Scalar& value)
{
  return DenseBase<Derived>::NullaryExpr(size, internal::scalar_constant_op<Scalar>(value));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::ConstantReturnType
DenseBase<Derived>::Constant(const Scalar& value)
{
  EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived)
  return DenseBase<Derived>::NullaryExpr(RowsAtCompileTime, ColsAtCompileTime, internal::scalar_constant_op<Scalar>(value));
}

template<typename Derived>
EIGEN_DEPRECATED EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::RandomAccessLinSpacedReturnType
DenseBase<Derived>::LinSpaced(Sequential_t, Index size, const Scalar& low, const Scalar& high)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return DenseBase<Derived>::NullaryExpr(size, internal::linspaced_op<Scalar>(low,high,size));
}

template<typename Derived>
EIGEN_DEPRECATED EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::RandomAccessLinSpacedReturnType
DenseBase<Derived>::LinSpaced(Sequential_t, const Scalar& low, const Scalar& high)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived)
  return DenseBase<Derived>::NullaryExpr(Derived::SizeAtCompileTime, internal::linspaced_op<Scalar>(low,high,Derived::SizeAtCompileTime));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::RandomAccessLinSpacedReturnType
DenseBase<Derived>::LinSpaced(Index size, const Scalar& low, const Scalar& high)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return DenseBase<Derived>::NullaryExpr(size, internal::linspaced_op<Scalar>(low,high,size));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::RandomAccessLinSpacedReturnType
DenseBase<Derived>::LinSpaced(const Scalar& low, const Scalar& high)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived)
  return DenseBase<Derived>::NullaryExpr(Derived::SizeAtCompileTime, internal::linspaced_op<Scalar>(low,high,Derived::SizeAtCompileTime));
}

template<typename Derived>
EIGEN_DEVICE_FUNC bool DenseBase<Derived>::isApproxToConstant
(const Scalar& val, const RealScalar& prec) const
{
  typename internal::nested_eval<Derived,1>::type self(derived());
  for(Index j = 0; j < cols(); ++j)
    for(Index i = 0; i < rows(); ++i)
      if(!internal::isApprox(self.coeff(i, j), val, prec))
        return false;
  return true;
}

template<typename Derived>
EIGEN_DEVICE_FUNC bool DenseBase<Derived>::isConstant
(const Scalar& val, const RealScalar& prec) const
{
  return isApproxToConstant(val, prec);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void DenseBase<Derived>::fill(const Scalar& val)
{
  setConstant(val);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& DenseBase<Derived>::setConstant(const Scalar& val)
{
  return derived() = Constant(rows(), cols(), val);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setConstant(Index size, const Scalar& val)
{
  resize(size);
  return setConstant(val);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setConstant(Index rows, Index cols, const Scalar& val)
{
  resize(rows, cols);
  return setConstant(val);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setConstant(NoChange_t, Index cols, const Scalar& val)
{
  return setConstant(rows(), cols, val);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setConstant(Index rows, NoChange_t, const Scalar& val)
{
  return setConstant(rows, cols(), val);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& DenseBase<Derived>::setLinSpaced(Index newSize, const Scalar& low, const Scalar& high)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return derived() = Derived::NullaryExpr(newSize, internal::linspaced_op<Scalar>(low,high,newSize));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& DenseBase<Derived>::setLinSpaced(const Scalar& low, const Scalar& high)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return setLinSpaced(size(), low, high);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::ConstantReturnType
DenseBase<Derived>::Zero(Index rows, Index cols)
{
  return Constant(rows, cols, Scalar(0));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::ConstantReturnType
DenseBase<Derived>::Zero(Index size)
{
  return Constant(size, Scalar(0));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::ConstantReturnType
DenseBase<Derived>::Zero()
{
  return Constant(Scalar(0));
}

template<typename Derived>
EIGEN_DEVICE_FUNC bool DenseBase<Derived>::isZero(const RealScalar& prec) const
{
  typename internal::nested_eval<Derived,1>::type self(derived());
  for(Index j = 0; j < cols(); ++j)
    for(Index i = 0; i < rows(); ++i)
      if(!internal::isMuchSmallerThan(self.coeff(i, j), static_cast<Scalar>(1), prec))
        return false;
  return true;
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& DenseBase<Derived>::setZero()
{
  return setConstant(Scalar(0));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setZero(Index newSize)
{
  resize(newSize);
  return setConstant(Scalar(0));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setZero(Index rows, Index cols)
{
  resize(rows, cols);
  return setConstant(Scalar(0));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setZero(NoChange_t, Index cols)
{
  return setZero(rows(), cols);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setZero(Index rows, NoChange_t)
{
  return setZero(rows, cols());
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::ConstantReturnType
DenseBase<Derived>::Ones(Index rows, Index cols)
{
  return Constant(rows, cols, Scalar(1));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::ConstantReturnType
DenseBase<Derived>::Ones(Index newSize)
{
  return Constant(newSize, Scalar(1));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename DenseBase<Derived>::ConstantReturnType
DenseBase<Derived>::Ones()
{
  return Constant(Scalar(1));
}

template<typename Derived>
EIGEN_DEVICE_FUNC bool DenseBase<Derived>::isOnes
(const RealScalar& prec) const
{
  return isApproxToConstant(Scalar(1), prec);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& DenseBase<Derived>::setOnes()
{
  return setConstant(Scalar(1));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setOnes(Index newSize)
{
  resize(newSize);
  return setConstant(Scalar(1));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setOnes(Index rows, Index cols)
{
  resize(rows, cols);
  return setConstant(Scalar(1));
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setOnes(Index rows, NoChange_t)
{
  return setOnes(rows, cols());
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setOnes(NoChange_t, Index cols)
{
  return setOnes(rows(), cols);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename MatrixBase<Derived>::IdentityReturnType
MatrixBase<Derived>::Identity(Index rows, Index cols)
{
  return DenseBase<Derived>::NullaryExpr(rows, cols, internal::scalar_identity_op<Scalar>());
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename MatrixBase<Derived>::IdentityReturnType
MatrixBase<Derived>::Identity()
{
  EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived)
  return MatrixBase<Derived>::NullaryExpr(RowsAtCompileTime, ColsAtCompileTime, internal::scalar_identity_op<Scalar>());
}

template<typename Derived>
bool MatrixBase<Derived>::isIdentity
(const RealScalar& prec) const
{
  typename internal::nested_eval<Derived,1>::type self(derived());
  for(Index j = 0; j < cols(); ++j)
  {
    for(Index i = 0; i < rows(); ++i)
    {
      if(i == j)
      {
        if(!internal::isApprox(self.coeff(i, j), static_cast<Scalar>(1), prec))
          return false;
      }
      else
      {
        if(!internal::isMuchSmallerThan(self.coeff(i, j), static_cast<RealScalar>(1), prec))
          return false;
      }
    }
  }
  return true;
}

namespace internal {

template<typename Derived, bool Big = (Derived::SizeAtCompileTime>=16)>
struct setIdentity_impl
{
  EIGEN_DEVICE_FUNC
  static EIGEN_STRONG_INLINE Derived& run(Derived& m)
  {
    return m = Derived::Identity(m.rows(), m.cols());
  }
};

template<typename Derived>
struct setIdentity_impl<Derived, true>
{
  EIGEN_DEVICE_FUNC
  static EIGEN_STRONG_INLINE Derived& run(Derived& m)
  {
    m.setZero();
    const Index size = numext::mini(m.rows(), m.cols());
    for(Index i = 0; i < size; ++i) m.coeffRef(i,i) = typename Derived::Scalar(1);
    return m;
  }
};

} 

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& MatrixBase<Derived>::setIdentity()
{
  return internal::setIdentity_impl<Derived>::run(derived());
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& MatrixBase<Derived>::setIdentity(Index rows, Index cols)
{
  derived().resize(rows, cols);
  return setIdentity();
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename MatrixBase<Derived>::BasisReturnType MatrixBase<Derived>::Unit(Index newSize, Index i)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return BasisReturnType(SquareMatrixType::Identity(newSize,newSize), i);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename MatrixBase<Derived>::BasisReturnType MatrixBase<Derived>::Unit(Index i)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return BasisReturnType(SquareMatrixType::Identity(),i);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename MatrixBase<Derived>::BasisReturnType MatrixBase<Derived>::UnitX()
{ return Derived::Unit(0); }

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename MatrixBase<Derived>::BasisReturnType MatrixBase<Derived>::UnitY()
{ return Derived::Unit(1); }

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename MatrixBase<Derived>::BasisReturnType MatrixBase<Derived>::UnitZ()
{ return Derived::Unit(2); }

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename MatrixBase<Derived>::BasisReturnType MatrixBase<Derived>::UnitW()
{ return Derived::Unit(3); }

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& MatrixBase<Derived>::setUnit(Index i)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
  eigen_assert(i<size());
  derived().setZero();
  derived().coeffRef(i) = Scalar(1);
  return derived();
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& MatrixBase<Derived>::setUnit(Index newSize, Index i)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
  eigen_assert(i<newSize);
  derived().resize(newSize);
  return setUnit(i);
}

} 

#endif 
