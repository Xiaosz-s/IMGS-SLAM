#ifndef EIGEN_BICGSTAB_H
#define EIGEN_BICGSTAB_H

namespace Eigen { 

namespace internal {

template<typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
bool bicgstab(const MatrixType& mat, const Rhs& rhs, Dest& x,
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

  Index n = mat.cols();
  VectorType r  = rhs - mat * x;
  VectorType r0 = r;
  
  RealScalar r0_sqnorm = r0.squaredNorm();
  RealScalar rhs_sqnorm = rhs.squaredNorm();
  if(rhs_sqnorm == 0)
  {
    x.setZero();
    return true;
  }
  Scalar rho    = 1;
  Scalar alpha  = 1;
  Scalar w      = 1;
  
  VectorType v = VectorType::Zero(n), p = VectorType::Zero(n);
  VectorType y(n),  z(n);
  VectorType kt(n), ks(n);

  VectorType s(n), t(n);

  RealScalar tol2 = tol*tol*rhs_sqnorm;
  RealScalar eps2 = NumTraits<Scalar>::epsilon()*NumTraits<Scalar>::epsilon();
  Index i = 0;
  Index restarts = 0;

  while ( r.squaredNorm() > tol2 && i<maxIters )
  {
    Scalar rho_old = rho;

    rho = r0.dot(r);
    if (abs(rho) < eps2*r0_sqnorm)
    {

      r  = rhs - mat * x;
      r0 = r;
      rho = r0_sqnorm = r.squaredNorm();
      if(restarts++ == 0)
        i = 0;
    }
    Scalar beta = (rho/rho_old) * (alpha / w);
    p = r + beta * (p - w * v);
    
    y = precond.solve(p);
    
    v.noalias() = mat * y;

    alpha = rho / r0.dot(v);
    s = r - alpha * v;

    z = precond.solve(s);
    t.noalias() = mat * z;

    RealScalar tmp = t.squaredNorm();
    if(tmp>RealScalar(0))
      w = t.dot(s) / tmp;
    else
      w = Scalar(0);
    x += alpha * y + w * z;
    r = s - w * t;
    ++i;
  }
  tol_error = sqrt(r.squaredNorm()/rhs_sqnorm);
  iters = i;
  return true; 
}

}

template< typename MatrixType_,
          typename Preconditioner_ = DiagonalPreconditioner<typename MatrixType_::Scalar> >
class BiCGSTAB;

namespace internal {

template< typename MatrixType_, typename Preconditioner_>
struct traits<BiCGSTAB<MatrixType_,Preconditioner_> >
{
  typedef MatrixType_ MatrixType;
  typedef Preconditioner_ Preconditioner;
};

}

template< typename MatrixType_, typename Preconditioner_>
class BiCGSTAB : public IterativeSolverBase<BiCGSTAB<MatrixType_,Preconditioner_> >
{
  typedef IterativeSolverBase<BiCGSTAB> Base;
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

  BiCGSTAB() : Base() {}

  template<typename MatrixDerived>
  explicit BiCGSTAB(const EigenBase<MatrixDerived>& A) : Base(A.derived()) {}

  ~BiCGSTAB() {}

  template<typename Rhs,typename Dest>
  void _solve_vector_with_guess_impl(const Rhs& b, Dest& x) const
  {    
    m_iterations = Base::maxIterations();
    m_error = Base::m_tolerance;
    
    bool ret = internal::bicgstab(matrix(), b, x, Base::m_preconditioner, m_iterations, m_error);

    m_info = (!ret) ? NumericalIssue
           : m_error <= Base::m_tolerance ? Success
           : NoConvergence;
  }

protected:

};

} 

#endif 
