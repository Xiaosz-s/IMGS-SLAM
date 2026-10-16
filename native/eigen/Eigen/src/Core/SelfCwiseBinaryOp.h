#ifndef EIGEN_SELFCWISEBINARYOP_H
#define EIGEN_SELFCWISEBINARYOP_H

namespace Eigen { 

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& DenseBase<Derived>::operator*=(const Scalar& other)
{
  internal::call_assignment(this->derived(), PlainObject::Constant(rows(),cols(),other), internal::mul_assign_op<Scalar,Scalar>());
  return derived();
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& ArrayBase<Derived>::operator+=(const Scalar& other)
{
  internal::call_assignment(this->derived(), PlainObject::Constant(rows(),cols(),other), internal::add_assign_op<Scalar,Scalar>());
  return derived();
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& ArrayBase<Derived>::operator-=(const Scalar& other)
{
  internal::call_assignment(this->derived(), PlainObject::Constant(rows(),cols(),other), internal::sub_assign_op<Scalar,Scalar>());
  return derived();
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived& DenseBase<Derived>::operator/=(const Scalar& other)
{
  internal::call_assignment(this->derived(), PlainObject::Constant(rows(),cols(),other), internal::div_assign_op<Scalar,Scalar>());
  return derived();
}

} 

#endif 
