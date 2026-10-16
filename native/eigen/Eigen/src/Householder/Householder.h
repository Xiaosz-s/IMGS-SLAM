#ifndef EIGEN_HOUSEHOLDER_H
#define EIGEN_HOUSEHOLDER_H

namespace Eigen { 

namespace internal {
template<int n> struct decrement_size
{
  enum {
    ret = n==Dynamic ? n : n-1
  };
};
}

template<typename Derived>
EIGEN_DEVICE_FUNC
void MatrixBase<Derived>::makeHouseholderInPlace(Scalar& tau, RealScalar& beta)
{
  VectorBlock<Derived, internal::decrement_size<Base::SizeAtCompileTime>::ret> essentialPart(derived(), 1, size()-1);
  makeHouseholder(essentialPart, tau, beta);
}

template<typename Derived>
template<typename EssentialPart>
EIGEN_DEVICE_FUNC
void MatrixBase<Derived>::makeHouseholder(
  EssentialPart& essential,
  Scalar& tau,
  RealScalar& beta) const
{
  using std::sqrt;
  using numext::conj;
  
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(EssentialPart)
  VectorBlock<const Derived, EssentialPart::SizeAtCompileTime> tail(derived(), 1, size()-1);
  
  RealScalar tailSqNorm = size()==1 ? RealScalar(0) : tail.squaredNorm();
  Scalar c0 = coeff(0);
  const RealScalar tol = (std::numeric_limits<RealScalar>::min)();

  if(tailSqNorm <= tol && numext::abs2(numext::imag(c0))<=tol)
  {
    tau = RealScalar(0);
    beta = numext::real(c0);
    essential.setZero();
  }
  else
  {
    beta = sqrt(numext::abs2(c0) + tailSqNorm);
    if (numext::real(c0)>=RealScalar(0))
      beta = -beta;
    essential = tail / (c0 - beta);
    tau = conj((beta - c0) / beta);
  }
}

template<typename Derived>
template<typename EssentialPart>
EIGEN_DEVICE_FUNC
void MatrixBase<Derived>::applyHouseholderOnTheLeft(
  const EssentialPart& essential,
  const Scalar& tau,
  Scalar* workspace)
{
  if(rows() == 1)
  {
    *this *= Scalar(1)-tau;
  }
  else if(tau!=Scalar(0))
  {
    Map<typename internal::plain_row_type<PlainObject>::type> tmp(workspace,cols());
    Block<Derived, EssentialPart::SizeAtCompileTime, Derived::ColsAtCompileTime> bottom(derived(), 1, 0, rows()-1, cols());
    tmp.noalias() = essential.adjoint() * bottom;
    tmp += this->row(0);
    this->row(0) -= tau * tmp;
    bottom.noalias() -= tau * essential * tmp;
  }
}

template<typename Derived>
template<typename EssentialPart>
EIGEN_DEVICE_FUNC
void MatrixBase<Derived>::applyHouseholderOnTheRight(
  const EssentialPart& essential,
  const Scalar& tau,
  Scalar* workspace)
{
  if(cols() == 1)
  {
    *this *= Scalar(1)-tau;
  }
  else if(tau!=Scalar(0))
  {
    Map<typename internal::plain_col_type<PlainObject>::type> tmp(workspace,rows());
    Block<Derived, Derived::RowsAtCompileTime, EssentialPart::SizeAtCompileTime> right(derived(), 0, 1, rows(), cols()-1);
    tmp.noalias() = right * essential;
    tmp += this->col(0);
    this->col(0) -= tau * tmp;
    right.noalias() -= tau * tmp * essential.adjoint();
  }
}

} 

#endif 
