#ifndef EIGEN_LEAST_SQUARE_CONJUGATE_GRADIENT_H
#define EIGEN_LEAST_SQUARE_CONJUGATE_GRADIENT_H

namespace Eigen { 

namespace internal {

template<typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
EIGEN_DONT_INLINE
void least_square_conjugate_gradient(const MatrixType& mat, const Rhs& rhs, Dest& x,
                                     const Preconditioner& precond, Index& iters,
                                     typename Dest::RealScalar& tol_error)
{
  using std::sqrt;
  using std::abs;
  typedef typename Dest::RealScalar RealScalar;
  typedef typename Dest::Scalar Scalar;
  typedef Matrix<Scalar,Dynamic,1> VectorType;
  
  RealScalar tol = tol_error;
  Index maxIters = iters;
  
  Index m = mat.rows(), n = mat.cols();

  VectorType residual        = rhs - mat * x;
  VectorType normal_residual = mat.adjoint() * residual;

  RealScalar rhsNorm2 = (mat.adjoint()*rhs).squaredNorm();
  if(rhsNorm2 == 0) 
  {
    x.setZero();
    iters = 0;
    tol_error = 0;
    return;
  }
  RealScalar threshold = tol*tol*rhsNorm2;
  RealScalar residualNorm2 = normal_residual.squaredNorm();
  if (residualNorm2 < threshold)
  {
    iters = 0;
    tol_error = sqrt(residualNorm2 / rhsNorm2);
    return;
  }
  
  VectorType p(n);
  p = precond.solve(normal_residual);                         

  VectorType z(n), tmp(m);
  RealScalar absNew = numext::real(normal_residual.dot(p));  
  Index i = 0;
  while(i < maxIters)
  {
    tmp.noalias() = mat * p;

    Scalar alpha = absNew / tmp.squaredNorm();      
    x += alpha * p;                                 
    residual -= alpha * tmp;                        
    normal_residual = mat.adjoint() * residual;     
    
    residualNorm2 = normal_residual.squaredNorm();
    if(residualNorm2 < threshold)
      break;
    
    z = precond.solve(normal_residual);             

    RealScalar absOld = absNew;
    absNew = numext::real(normal_residual.dot(z));  
    RealScalar beta = absNew / absOld;              
    p = z + beta * p;                               
    i++;
  }
  tol_error = sqrt(residualNorm2 / rhsNorm2);
  iters = i;
}

}

template< typename MatrixType_,
          typename Preconditioner_ = LeastSquareDiagonalPreconditioner<typename MatrixType_::Scalar> >
class LeastSquaresConjugateGradient;

namespace internal {

template< typename MatrixType_, typename Preconditioner_>
struct traits<LeastSquaresConjugateGradient<MatrixType_,Preconditioner_> >
{
  typedef MatrixType_ MatrixType;
  typedef Preconditioner_ Preconditioner;
};

}

template< typename MatrixType_, typename Preconditioner_>
class LeastSquaresConjugateGradient : public IterativeSolverBase<LeastSquaresConjugateGradient<MatrixType_,Preconditioner_> >
{
  typedef IterativeSolverBase<LeastSquaresConjugateGradient> Base;
  using Base::matrix;
  using Base::m_error;
  using Base::m_iterations;
  using Base::m_info;
  using Base::m_isInitialized;
public:
  typedef MatrixType_ MatrixType;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;
  typedef Preconditioner_ Preconditioner;

public:

  LeastSquaresConjugateGradient() : Base() {}

  template<typename MatrixDerived>
  explicit LeastSquaresConjugateGradient(const EigenBase<MatrixDerived>& A) : Base(A.derived()) {}

  ~LeastSquaresConjugateGradient() {}

  template<typename Rhs,typename Dest>
  void _solve_vector_with_guess_impl(const Rhs& b, Dest& x) const
  {
    m_iterations = Base::maxIterations();
    m_error = Base::m_tolerance;

    internal::least_square_conjugate_gradient(matrix(), b, x, Base::m_preconditioner, m_iterations, m_error);
    m_info = m_error <= Base::m_tolerance ? Success : NoConvergence;
  }

};

} 

#endif 
