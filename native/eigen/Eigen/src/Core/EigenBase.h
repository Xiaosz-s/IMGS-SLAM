#ifndef EIGEN_EIGENBASE_H
#define EIGEN_EIGENBASE_H

namespace Eigen {

template<typename Derived> struct EigenBase
{

  typedef Eigen::Index Index;

  typedef typename internal::traits<Derived>::StorageKind StorageKind;

  EIGEN_DEVICE_FUNC
  Derived& derived() { return *static_cast<Derived*>(this); }
  
  EIGEN_DEVICE_FUNC
  const Derived& derived() const { return *static_cast<const Derived*>(this); }

  EIGEN_DEVICE_FUNC
  inline Derived& const_cast_derived() const
  { return *static_cast<Derived*>(const_cast<EigenBase*>(this)); }
  EIGEN_DEVICE_FUNC
  inline const Derived& const_derived() const
  { return *static_cast<const Derived*>(this); }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  inline Index rows() const EIGEN_NOEXCEPT { return derived().rows(); }
  
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  inline Index cols() const EIGEN_NOEXCEPT { return derived().cols(); }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  inline Index size() const EIGEN_NOEXCEPT { return rows() * cols(); }

  template<typename Dest>
  EIGEN_DEVICE_FUNC
  inline void evalTo(Dest& dst) const
  { derived().evalTo(dst); }

  template<typename Dest>
  EIGEN_DEVICE_FUNC
  inline void addTo(Dest& dst) const
  {

    typename Dest::PlainObject res(rows(),cols());
    evalTo(res);
    dst += res;
  }

  template<typename Dest>
  EIGEN_DEVICE_FUNC
  inline void subTo(Dest& dst) const
  {

    typename Dest::PlainObject res(rows(),cols());
    evalTo(res);
    dst -= res;
  }

  template<typename Dest>
  EIGEN_DEVICE_FUNC inline void applyThisOnTheRight(Dest& dst) const
  {

    dst = dst * this->derived();
  }

  template<typename Dest>
  EIGEN_DEVICE_FUNC inline void applyThisOnTheLeft(Dest& dst) const
  {

    dst = this->derived() * dst;
  }

};

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC
Derived& DenseBase<Derived>::operator=(const EigenBase<OtherDerived> &other)
{
  call_assignment(derived(), other.derived());
  return derived();
}

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC
Derived& DenseBase<Derived>::operator+=(const EigenBase<OtherDerived> &other)
{
  call_assignment(derived(), other.derived(), internal::add_assign_op<Scalar,typename OtherDerived::Scalar>());
  return derived();
}

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC
Derived& DenseBase<Derived>::operator-=(const EigenBase<OtherDerived> &other)
{
  call_assignment(derived(), other.derived(), internal::sub_assign_op<Scalar,typename OtherDerived::Scalar>());
  return derived();
}

} 

#endif 
