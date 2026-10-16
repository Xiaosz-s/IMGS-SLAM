#ifndef EIGEN_DOT_H
#define EIGEN_DOT_H

namespace Eigen { 

namespace internal {

template<typename T, typename U,

         bool NeedToTranspose = T::IsVectorAtCompileTime
                && U::IsVectorAtCompileTime
                && ((int(T::RowsAtCompileTime) == 1 && int(U::ColsAtCompileTime) == 1)
                      |  
                         
                    (int(T::ColsAtCompileTime) == 1 && int(U::RowsAtCompileTime) == 1))
>
struct dot_nocheck
{
  typedef scalar_conj_product_op<typename traits<T>::Scalar,typename traits<U>::Scalar> conj_prod;
  typedef typename conj_prod::result_type ResScalar;
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE
  static ResScalar run(const MatrixBase<T>& a, const MatrixBase<U>& b)
  {
    return a.template binaryExpr<conj_prod>(b).sum();
  }
};

template<typename T, typename U>
struct dot_nocheck<T, U, true>
{
  typedef scalar_conj_product_op<typename traits<T>::Scalar,typename traits<U>::Scalar> conj_prod;
  typedef typename conj_prod::result_type ResScalar;
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE
  static ResScalar run(const MatrixBase<T>& a, const MatrixBase<U>& b)
  {
    return a.transpose().template binaryExpr<conj_prod>(b).sum();
  }
};

} 

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE
typename ScalarBinaryOpTraits<typename internal::traits<Derived>::Scalar,typename internal::traits<OtherDerived>::Scalar>::ReturnType
MatrixBase<Derived>::dot(const MatrixBase<OtherDerived>& other) const
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(OtherDerived)
  EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Derived,OtherDerived)
#if !(defined(EIGEN_NO_STATIC_ASSERT) && defined(EIGEN_NO_DEBUG))
  typedef internal::scalar_conj_product_op<Scalar,typename OtherDerived::Scalar> func;
  EIGEN_CHECK_BINARY_COMPATIBILIY(func,Scalar,typename OtherDerived::Scalar);
#endif
  
  eigen_assert(size() == other.size());

  return internal::dot_nocheck<Derived,OtherDerived>::run(*this, other);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename NumTraits<typename internal::traits<Derived>::Scalar>::Real MatrixBase<Derived>::squaredNorm() const
{
  return numext::real((*this).cwiseAbs2().sum());
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename NumTraits<typename internal::traits<Derived>::Scalar>::Real MatrixBase<Derived>::norm() const
{
  return numext::sqrt(squaredNorm());
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename MatrixBase<Derived>::PlainObject
MatrixBase<Derived>::normalized() const
{
  typedef typename internal::nested_eval<Derived,2>::type _Nested;
  _Nested n(derived());
  RealScalar z = n.squaredNorm();
  
  if(z>RealScalar(0))
    return n / numext::sqrt(z);
  else
    return n;
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void MatrixBase<Derived>::normalize()
{
  RealScalar z = squaredNorm();
  
  if(z>RealScalar(0))
    derived() /= numext::sqrt(z);
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename MatrixBase<Derived>::PlainObject
MatrixBase<Derived>::stableNormalized() const
{
  typedef typename internal::nested_eval<Derived,3>::type _Nested;
  _Nested n(derived());
  RealScalar w = n.cwiseAbs().maxCoeff();
  RealScalar z = (n/w).squaredNorm();
  if(z>RealScalar(0))
    return n / (numext::sqrt(z)*w);
  else
    return n;
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void MatrixBase<Derived>::stableNormalize()
{
  RealScalar w = cwiseAbs().maxCoeff();
  RealScalar z = (derived()/w).squaredNorm();
  if(z>RealScalar(0))
    derived() /= numext::sqrt(z)*w;
}

namespace internal {

template<typename Derived, int p>
struct lpNorm_selector
{
  typedef typename NumTraits<typename traits<Derived>::Scalar>::Real RealScalar;
  EIGEN_DEVICE_FUNC
  static inline RealScalar run(const MatrixBase<Derived>& m)
  {
    EIGEN_USING_STD(pow)
    return pow(m.cwiseAbs().array().pow(p).sum(), RealScalar(1)/p);
  }
};

template<typename Derived>
struct lpNorm_selector<Derived, 1>
{
  EIGEN_DEVICE_FUNC
  static inline typename NumTraits<typename traits<Derived>::Scalar>::Real run(const MatrixBase<Derived>& m)
  {
    return m.cwiseAbs().sum();
  }
};

template<typename Derived>
struct lpNorm_selector<Derived, 2>
{
  EIGEN_DEVICE_FUNC
  static inline typename NumTraits<typename traits<Derived>::Scalar>::Real run(const MatrixBase<Derived>& m)
  {
    return m.norm();
  }
};

template<typename Derived>
struct lpNorm_selector<Derived, Infinity>
{
  typedef typename NumTraits<typename traits<Derived>::Scalar>::Real RealScalar;
  EIGEN_DEVICE_FUNC
  static inline RealScalar run(const MatrixBase<Derived>& m)
  {
    if(Derived::SizeAtCompileTime==0 || (Derived::SizeAtCompileTime==Dynamic && m.size()==0))
      return RealScalar(0);
    return m.cwiseAbs().maxCoeff();
  }
};

} 

template<typename Derived>
template<int p>
#ifndef EIGEN_PARSED_BY_DOXYGEN
EIGEN_DEVICE_FUNC inline typename NumTraits<typename internal::traits<Derived>::Scalar>::Real
#else
EIGEN_DEVICE_FUNC MatrixBase<Derived>::RealScalar
#endif
MatrixBase<Derived>::lpNorm() const
{
  return internal::lpNorm_selector<Derived, p>::run(*this);
}

template<typename Derived>
template<typename OtherDerived>
bool MatrixBase<Derived>::isOrthogonal
(const MatrixBase<OtherDerived>& other, const RealScalar& prec) const
{
  typename internal::nested_eval<Derived,2>::type nested(derived());
  typename internal::nested_eval<OtherDerived,2>::type otherNested(other.derived());
  return numext::abs2(nested.dot(otherNested)) <= prec * prec * nested.squaredNorm() * otherNested.squaredNorm();
}

template<typename Derived>
bool MatrixBase<Derived>::isUnitary(const RealScalar& prec) const
{
  typename internal::nested_eval<Derived,1>::type self(derived());
  for(Index i = 0; i < cols(); ++i)
  {
    if(!internal::isApprox(self.col(i).squaredNorm(), static_cast<RealScalar>(1), prec))
      return false;
    for(Index j = 0; j < i; ++j)
      if(!internal::isMuchSmallerThan(self.col(i).dot(self.col(j)), static_cast<Scalar>(1), prec))
        return false;
  }
  return true;
}

} 

#endif 
