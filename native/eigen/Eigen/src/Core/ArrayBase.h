#ifndef EIGEN_ARRAYBASE_H
#define EIGEN_ARRAYBASE_H

namespace Eigen { 

template<typename ExpressionType> class MatrixWrapper;

template<typename Derived> class ArrayBase
  : public DenseBase<Derived>
{
  public:
#ifndef EIGEN_PARSED_BY_DOXYGEN
    
    typedef ArrayBase StorageBaseType;

    typedef ArrayBase Eigen_BaseClassForSpecializationOfGlobalMathFuncImpl;

    typedef typename internal::traits<Derived>::StorageKind StorageKind;
    typedef typename internal::traits<Derived>::Scalar Scalar;
    typedef typename internal::packet_traits<Scalar>::type PacketScalar;
    typedef typename NumTraits<Scalar>::Real RealScalar;

    typedef DenseBase<Derived> Base;
    using Base::RowsAtCompileTime;
    using Base::ColsAtCompileTime;
    using Base::SizeAtCompileTime;
    using Base::MaxRowsAtCompileTime;
    using Base::MaxColsAtCompileTime;
    using Base::MaxSizeAtCompileTime;
    using Base::IsVectorAtCompileTime;
    using Base::Flags;
    
    using Base::derived;
    using Base::const_cast_derived;
    using Base::rows;
    using Base::cols;
    using Base::size;
    using Base::coeff;
    using Base::coeffRef;
    using Base::lazyAssign;
    using Base::operator-;
    using Base::operator=;
    using Base::operator+=;
    using Base::operator-=;
    using Base::operator*=;
    using Base::operator/=;

    typedef typename Base::CoeffReturnType CoeffReturnType;

#endif 

#ifndef EIGEN_PARSED_BY_DOXYGEN
    typedef typename Base::PlainObject PlainObject;

    typedef CwiseNullaryOp<internal::scalar_constant_op<Scalar>,PlainObject> ConstantReturnType;
#endif 

#define EIGEN_CURRENT_STORAGE_BASE_CLASS Eigen::ArrayBase
#define EIGEN_DOC_UNARY_ADDONS(X,Y)
#   include "../plugins/MatrixCwiseUnaryOps.h"
#   include "../plugins/ArrayCwiseUnaryOps.h"
#   include "../plugins/CommonCwiseBinaryOps.h"
#   include "../plugins/MatrixCwiseBinaryOps.h"
#   include "../plugins/ArrayCwiseBinaryOps.h"
#   ifdef EIGEN_ARRAYBASE_PLUGIN
#     include EIGEN_ARRAYBASE_PLUGIN
#   endif
#undef EIGEN_CURRENT_STORAGE_BASE_CLASS
#undef EIGEN_DOC_UNARY_ADDONS

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator=(const ArrayBase& other)
    {
      internal::call_assignment(derived(), other.derived());
      return derived();
    }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator=(const Scalar &value)
    { Base::setConstant(value); return derived(); }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator+=(const Scalar& scalar);
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator-=(const Scalar& scalar);

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator+=(const ArrayBase<OtherDerived>& other);
    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator-=(const ArrayBase<OtherDerived>& other);

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator*=(const ArrayBase<OtherDerived>& other);

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator/=(const ArrayBase<OtherDerived>& other);

  public:
    EIGEN_DEVICE_FUNC
    ArrayBase<Derived>& array() { return *this; }
    EIGEN_DEVICE_FUNC
    const ArrayBase<Derived>& array() const { return *this; }

    EIGEN_DEVICE_FUNC
    MatrixWrapper<Derived> matrix() { return MatrixWrapper<Derived>(derived()); }
    EIGEN_DEVICE_FUNC
    const MatrixWrapper<const Derived> matrix() const { return MatrixWrapper<const Derived>(derived()); }

  protected:
    EIGEN_DEFAULT_COPY_CONSTRUCTOR(ArrayBase)
    EIGEN_DEFAULT_EMPTY_CONSTRUCTOR_AND_DESTRUCTOR(ArrayBase)

  private:
    explicit ArrayBase(Index);
    ArrayBase(Index,Index);
    template<typename OtherDerived> explicit ArrayBase(const ArrayBase<OtherDerived>&);
  protected:
    
    template<typename OtherDerived> Derived& operator+=(const MatrixBase<OtherDerived>& )
    {EIGEN_STATIC_ASSERT(std::ptrdiff_t(sizeof(typename OtherDerived::Scalar))==-1,YOU_CANNOT_MIX_ARRAYS_AND_MATRICES); return *this;}
    
    template<typename OtherDerived> Derived& operator-=(const MatrixBase<OtherDerived>& )
    {EIGEN_STATIC_ASSERT(std::ptrdiff_t(sizeof(typename OtherDerived::Scalar))==-1,YOU_CANNOT_MIX_ARRAYS_AND_MATRICES); return *this;}
};

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived &
ArrayBase<Derived>::operator-=(const ArrayBase<OtherDerived> &other)
{
  call_assignment(derived(), other.derived(), internal::sub_assign_op<Scalar,typename OtherDerived::Scalar>());
  return derived();
}

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived &
ArrayBase<Derived>::operator+=(const ArrayBase<OtherDerived>& other)
{
  call_assignment(derived(), other.derived(), internal::add_assign_op<Scalar,typename OtherDerived::Scalar>());
  return derived();
}

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived &
ArrayBase<Derived>::operator*=(const ArrayBase<OtherDerived>& other)
{
  call_assignment(derived(), other.derived(), internal::mul_assign_op<Scalar,typename OtherDerived::Scalar>());
  return derived();
}

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived &
ArrayBase<Derived>::operator/=(const ArrayBase<OtherDerived>& other)
{
  call_assignment(derived(), other.derived(), internal::div_assign_op<Scalar,typename OtherDerived::Scalar>());
  return derived();
}

} 

#endif 
