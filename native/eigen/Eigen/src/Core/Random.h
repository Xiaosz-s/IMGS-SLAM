#ifndef EIGEN_RANDOM_H
#define EIGEN_RANDOM_H

namespace Eigen { 

namespace internal {

template<typename Scalar> struct scalar_random_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_random_op)
  inline const Scalar operator() () const { return random<Scalar>(); }
};

template<typename Scalar>
struct functor_traits<scalar_random_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };

} 

template<typename Derived>
inline const typename DenseBase<Derived>::RandomReturnType
DenseBase<Derived>::Random(Index rows, Index cols)
{
  return NullaryExpr(rows, cols, internal::scalar_random_op<Scalar>());
}

template<typename Derived>
inline const typename DenseBase<Derived>::RandomReturnType
DenseBase<Derived>::Random(Index size)
{
  return NullaryExpr(size, internal::scalar_random_op<Scalar>());
}

template<typename Derived>
inline const typename DenseBase<Derived>::RandomReturnType
DenseBase<Derived>::Random()
{
  return NullaryExpr(RowsAtCompileTime, ColsAtCompileTime, internal::scalar_random_op<Scalar>());
}

template<typename Derived>
EIGEN_DEVICE_FUNC inline Derived& DenseBase<Derived>::setRandom()
{
  return *this = Random(rows(), cols());
}

template<typename Derived>
EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setRandom(Index newSize)
{
  resize(newSize);
  return setRandom();
}

template<typename Derived>
EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setRandom(Index rows, Index cols)
{
  resize(rows, cols);
  return setRandom();
}

template<typename Derived>
EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setRandom(NoChange_t, Index cols)
{
  return setRandom(rows(), cols);
}

template<typename Derived>
EIGEN_STRONG_INLINE Derived&
PlainObjectBase<Derived>::setRandom(Index rows, NoChange_t)
{
  return setRandom(rows, cols());
}

} 

#endif 
