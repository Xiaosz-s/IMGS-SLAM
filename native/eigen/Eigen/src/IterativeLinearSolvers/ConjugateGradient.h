#ifndef EIGEN_CONJUGATE_GRADIENT_H
#define EIGEN_CONJUGATE_GRADIENT_H

namespace Eigen { 

namespace internal {

template<typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
EIGEN_DONT_INLINE
void conjugate_gradient(const MatrixType& mat, const Rhs& rhs, Dest& x,
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

  VectorType residual = rhs - mat * x; 

  RealScalar rhsNorm2 = rhs.squaredNorm();
  if(rhsNorm2 == 0) 
  {
    x.setZero();
    iters = 0;
    tol_error = 0;
    return;
  }
  const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();
  RealScalar threshold = numext::maxi(RealScalar(tol*tol*rhsNorm2),considerAsZero);
  RealScalar residualNorm2 = residual.squaredNorm();
  if (residualNorm2 < threshold)
  {
    iters = 0;
    tol_error = sqrt(residualNorm2 / rhsNorm2);
    return;
  }

  VectorType p(n);
  p = precond.solve(residual);      

  VectorType z(n), tmp(n);
  RealScalar absNew = numext::real(residual.dot(p));  
  Index i = 0;
  while(i < maxIters)
  {
    tmp.noalias() = mat * p;                    

    Scalar alpha = absNew / p.dot(tmp);         
    x += alpha * p;                             
    residual -= alpha * tmp;                    
    
    residualNorm2 = residual.squaredNorm();
    if(residualNorm2 < threshold)
      break;
    
    z = precond.solve(residual);                

    RealScalar absOld = absNew;
    absNew = numext::real(residual.dot(z));     
    RealScalar beta = absNew / absOld;          
    p = z + beta * p;                           
    i++;
  }
  tol_error = sqrt(residualNorm2 / rhsNorm2);
  iters = i;
}

}

template< typename MatrixType_, int UpLo_=Lower,
          typename Preconditioner_ = DiagonalPreconditioner<typename MatrixType_::Scalar> >
class ConjugateGradient;

namespace internal {

template< typename MatrixType_, int UpLo_, typename Preconditioner_>
struct traits<ConjugateGradient<MatrixType_,UpLo_,Preconditioner_> >
{
  typedef MatrixType_ MatrixType;
  typedef Preconditioner_ Preconditioner;
};

}

template< typename MatrixType_, int UpLo_, typename Preconditioner_>
class ConjugateGradient : public IterativeSolverBase<ConjugateGradient<MatrixType_,UpLo_,Preconditioner_> >
{
  typedef IterativeSolverBase<ConjugateGradient> Base;
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

  enum {
    UpLo = UpLo_
  };

public:

  ConjugateGradient() : Base() {}

  template<typename MatrixDerived>
  explicit ConjugateGradient(const EigenBase<MatrixDerived>& A) : Base(A.derived()) {}

  ~ConjugateGradient() {}

  template<typename Rhs,typename Dest>
  void _solve_vector_with_guess_impl(const Rhs& b, Dest& x) const
  {
    typedef typename Base::MatrixWrapper MatrixWrapper;
    typedef typename Base::ActualMatrixType ActualMatrixType;
    enum {
      TransposeInput  =   (!MatrixWrapper::MatrixFree)
                      &&  (UpLo==(Lower|Upper))
                      &&  (!MatrixType::IsRowMajor)
                      &&  (!NumTraits<Scalar>::IsComplex)
    };
    typedef typename internal::conditional<TransposeInput,Transpose<const ActualMatrixType>, ActualMatrixType const&>::type RowMajorWrapper;
    EIGEN_STATIC_ASSERT(EIGEN_IMPLIES(MatrixWrapper::MatrixFree,UpLo==(Lower|Upper)),MATRIX_FREE_CONJUGATE_GRADIENT_IS_COMPATIBLE_WITH_UPPER_UNION_LOWER_MODE_ONLY);
    typedef typename internal::conditional<UpLo==(Lower|Upper),
                                           RowMajorWrapper,
                                           typename MatrixWrapper::template ConstSelfAdjointViewReturnType<UpLo>::Type
                                          >::type SelfAdjointWrapper;

    m_iterations = Base::maxIterations();
    m_error = Base::m_tolerance;

    RowMajorWrapper row_mat(matrix());
    internal::conjugate_gradient(SelfAdjointWrapper(row_mat), b, x, Base::m_preconditioner, m_iterations, m_error);
    m_info = m_error <= Base::m_tolerance ? Success : NoConvergence;
  }

protected:

};

} 

#endif 
