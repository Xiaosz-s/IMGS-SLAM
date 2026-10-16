#ifndef EIGEN_SPARSE_FUZZY_H
#define EIGEN_SPARSE_FUZZY_H

namespace Eigen {
  
template<typename Derived>
template<typename OtherDerived>
bool SparseMatrixBase<Derived>::isApprox(const SparseMatrixBase<OtherDerived>& other, const RealScalar &prec) const
{
  const typename internal::nested_eval<Derived,2,PlainObject>::type actualA(derived());
  typename internal::conditional<bool(IsRowMajor)==bool(OtherDerived::IsRowMajor),
    const typename internal::nested_eval<OtherDerived,2,PlainObject>::type,
    const PlainObject>::type actualB(other.derived());

  return (actualA - actualB).squaredNorm() <= prec * prec * numext::mini(actualA.squaredNorm(), actualB.squaredNorm());
}

} 

#endif 
