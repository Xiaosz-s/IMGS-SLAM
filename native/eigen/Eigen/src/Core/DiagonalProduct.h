#ifndef EIGEN_DIAGONALPRODUCT_H
#define EIGEN_DIAGONALPRODUCT_H

namespace Eigen { 

template<typename Derived>
template<typename DiagonalDerived>
EIGEN_DEVICE_FUNC inline const Product<Derived, DiagonalDerived, LazyProduct>
MatrixBase<Derived>::operator*(const DiagonalBase<DiagonalDerived> &a_diagonal) const
{
  return Product<Derived, DiagonalDerived, LazyProduct>(derived(),a_diagonal.derived());
}

} 

#endif 
