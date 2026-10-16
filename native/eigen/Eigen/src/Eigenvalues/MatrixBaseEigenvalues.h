#ifndef EIGEN_MATRIXBASEEIGENVALUES_H
#define EIGEN_MATRIXBASEEIGENVALUES_H

namespace Eigen { 

namespace internal {

template<typename Derived, bool IsComplex>
struct eigenvalues_selector
{
  
  static inline typename MatrixBase<Derived>::EigenvaluesReturnType const
  run(const MatrixBase<Derived>& m)
  {
    typedef typename Derived::PlainObject PlainObject;
    PlainObject m_eval(m);
    return ComplexEigenSolver<PlainObject>(m_eval, false).eigenvalues();
  }
};

template<typename Derived>
struct eigenvalues_selector<Derived, false>
{
  static inline typename MatrixBase<Derived>::EigenvaluesReturnType const
  run(const MatrixBase<Derived>& m)
  {
    typedef typename Derived::PlainObject PlainObject;
    PlainObject m_eval(m);
    return EigenSolver<PlainObject>(m_eval, false).eigenvalues();
  }
};

} 

template<typename Derived>
inline typename MatrixBase<Derived>::EigenvaluesReturnType
MatrixBase<Derived>::eigenvalues() const
{
  return internal::eigenvalues_selector<Derived, NumTraits<Scalar>::IsComplex>::run(derived());
}

template<typename MatrixType, unsigned int UpLo> 
EIGEN_DEVICE_FUNC inline typename SelfAdjointView<MatrixType, UpLo>::EigenvaluesReturnType
SelfAdjointView<MatrixType, UpLo>::eigenvalues() const
{
  PlainObject thisAsMatrix(*this);
  return SelfAdjointEigenSolver<PlainObject>(thisAsMatrix, false).eigenvalues();
}

template<typename Derived>
inline typename MatrixBase<Derived>::RealScalar
MatrixBase<Derived>::operatorNorm() const
{
  using std::sqrt;
  typename Derived::PlainObject m_eval(derived());

  return sqrt((m_eval*m_eval.adjoint())
                 .eval()
		 .template selfadjointView<Lower>()
		 .eigenvalues()
		 .maxCoeff()
		 );
}

template<typename MatrixType, unsigned int UpLo>
EIGEN_DEVICE_FUNC inline typename SelfAdjointView<MatrixType, UpLo>::RealScalar
SelfAdjointView<MatrixType, UpLo>::operatorNorm() const
{
  return eigenvalues().cwiseAbs().maxCoeff();
}

} 

#endif
