#ifndef EIGEN_ITERATIVE_SOLVER_BASE_H
#define EIGEN_ITERATIVE_SOLVER_BASE_H

namespace Eigen {

namespace internal {

template<typename MatrixType>
struct is_ref_compatible_impl
{
private:
  template <typename T0>
  struct any_conversion
  {
    template <typename T> any_conversion(const volatile T&);
    template <typename T> any_conversion(T&);
  };
  struct yes {int a[1];};
  struct no  {int a[2];};

  template<typename T>
  static yes test(const Ref<const T>&, int);
  template<typename T>
  static no  test(any_conversion<T>, ...);

public:
  static MatrixType ms_from;
  enum { value = sizeof(test<MatrixType>(ms_from, 0))==sizeof(yes) };
};

template<typename MatrixType>
struct is_ref_compatible
{
  enum { value = is_ref_compatible_impl<typename remove_all<MatrixType>::type>::value };
};

template<typename MatrixType, bool MatrixFree = !internal::is_ref_compatible<MatrixType>::value>
class generic_matrix_wrapper;

template<typename MatrixType>
class generic_matrix_wrapper<MatrixType,false>
{
public:
  typedef Ref<const MatrixType> ActualMatrixType;
  template<int UpLo> struct ConstSelfAdjointViewReturnType {
    typedef typename ActualMatrixType::template ConstSelfAdjointViewReturnType<UpLo>::Type Type;
  };

  enum {
    MatrixFree = false
  };

  generic_matrix_wrapper()
    : m_dummy(0,0), m_matrix(m_dummy)
  {}

  template<typename InputType>
  generic_matrix_wrapper(const InputType &mat)
    : m_matrix(mat)
  {}

  const ActualMatrixType& matrix() const
  {
    return m_matrix;
  }

  template<typename MatrixDerived>
  void grab(const EigenBase<MatrixDerived> &mat)
  {
    m_matrix.~Ref<const MatrixType>();
    ::new (&m_matrix) Ref<const MatrixType>(mat.derived());
  }

  void grab(const Ref<const MatrixType> &mat)
  {
    if(&(mat.derived()) != &m_matrix)
    {
      m_matrix.~Ref<const MatrixType>();
      ::new (&m_matrix) Ref<const MatrixType>(mat);
    }
  }

protected:
  MatrixType m_dummy; 
  ActualMatrixType m_matrix;
};

template<typename MatrixType>
class generic_matrix_wrapper<MatrixType,true>
{
public:
  typedef MatrixType ActualMatrixType;
  template<int UpLo> struct ConstSelfAdjointViewReturnType
  {
    typedef ActualMatrixType Type;
  };

  enum {
    MatrixFree = true
  };

  generic_matrix_wrapper()
    : mp_matrix(0)
  {}

  generic_matrix_wrapper(const MatrixType &mat)
    : mp_matrix(&mat)
  {}

  const ActualMatrixType& matrix() const
  {
    return *mp_matrix;
  }

  void grab(const MatrixType &mat)
  {
    mp_matrix = &mat;
  }

protected:
  const ActualMatrixType *mp_matrix;
};

}

template< typename Derived>
class IterativeSolverBase : public SparseSolverBase<Derived>
{
protected:
  typedef SparseSolverBase<Derived> Base;
  using Base::m_isInitialized;

public:
  typedef typename internal::traits<Derived>::MatrixType MatrixType;
  typedef typename internal::traits<Derived>::Preconditioner Preconditioner;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::StorageIndex StorageIndex;
  typedef typename MatrixType::RealScalar RealScalar;

  enum {
    ColsAtCompileTime = MatrixType::ColsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
  };

public:

  using Base::derived;

  IterativeSolverBase()
  {
    init();
  }

  template<typename MatrixDerived>
  explicit IterativeSolverBase(const EigenBase<MatrixDerived>& A)
    : m_matrixWrapper(A.derived())
  {
    init();
    compute(matrix());
  }

  ~IterativeSolverBase() {}

  template<typename MatrixDerived>
  Derived& analyzePattern(const EigenBase<MatrixDerived>& A)
  {
    grab(A.derived());
    m_preconditioner.analyzePattern(matrix());
    m_isInitialized = true;
    m_analysisIsOk = true;
    m_info = m_preconditioner.info();
    return derived();
  }

  template<typename MatrixDerived>
  Derived& factorize(const EigenBase<MatrixDerived>& A)
  {
    eigen_assert(m_analysisIsOk && "You must first call analyzePattern()");
    grab(A.derived());
    m_preconditioner.factorize(matrix());
    m_factorizationIsOk = true;
    m_info = m_preconditioner.info();
    return derived();
  }

  template<typename MatrixDerived>
  Derived& compute(const EigenBase<MatrixDerived>& A)
  {
    grab(A.derived());
    m_preconditioner.compute(matrix());
    m_isInitialized = true;
    m_analysisIsOk = true;
    m_factorizationIsOk = true;
    m_info = m_preconditioner.info();
    return derived();
  }

  EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return matrix().rows(); }

  EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return matrix().cols(); }

  RealScalar tolerance() const { return m_tolerance; }

  Derived& setTolerance(const RealScalar& tolerance)
  {
    m_tolerance = tolerance;
    return derived();
  }

  Preconditioner& preconditioner() { return m_preconditioner; }

  const Preconditioner& preconditioner() const { return m_preconditioner; }

  Index maxIterations() const
  {
    return (m_maxIterations<0) ? 2*matrix().cols() : m_maxIterations;
  }

  Derived& setMaxIterations(Index maxIters)
  {
    m_maxIterations = maxIters;
    return derived();
  }

  Index iterations() const
  {
    eigen_assert(m_isInitialized && "ConjugateGradient is not initialized.");
    return m_iterations;
  }

  RealScalar error() const
  {
    eigen_assert(m_isInitialized && "ConjugateGradient is not initialized.");
    return m_error;
  }

  template<typename Rhs,typename Guess>
  inline const SolveWithGuess<Derived, Rhs, Guess>
  solveWithGuess(const MatrixBase<Rhs>& b, const Guess& x0) const
  {
    eigen_assert(m_isInitialized && "Solver is not initialized.");
    eigen_assert(derived().rows()==b.rows() && "solve(): invalid number of rows of the right hand side matrix b");
    return SolveWithGuess<Derived, Rhs, Guess>(derived(), b.derived(), x0);
  }

  ComputationInfo info() const
  {
    eigen_assert(m_isInitialized && "IterativeSolverBase is not initialized.");
    return m_info;
  }

  template<typename Rhs, typename DestDerived>
  void _solve_with_guess_impl(const Rhs& b, SparseMatrixBase<DestDerived> &aDest) const
  {
    eigen_assert(rows()==b.rows());

    Index rhsCols = b.cols();
    Index size = b.rows();
    DestDerived& dest(aDest.derived());
    typedef typename DestDerived::Scalar DestScalar;
    Eigen::Matrix<DestScalar,Dynamic,1> tb(size);
    Eigen::Matrix<DestScalar,Dynamic,1> tx(cols());

    typename DestDerived::PlainObject tmp(cols(),rhsCols);
    ComputationInfo global_info = Success;
    for(Index k=0; k<rhsCols; ++k)
    {
      tb = b.col(k);
      tx = dest.col(k);
      derived()._solve_vector_with_guess_impl(tb,tx);
      tmp.col(k) = tx.sparseView(0);

      if(m_info==NumericalIssue)
        global_info = NumericalIssue;
      else if(m_info==NoConvergence)
        global_info = NoConvergence;
    }
    m_info = global_info;
    dest.swap(tmp);
  }

  template<typename Rhs, typename DestDerived>
  typename internal::enable_if<Rhs::ColsAtCompileTime!=1 && DestDerived::ColsAtCompileTime!=1>::type
  _solve_with_guess_impl(const Rhs& b, MatrixBase<DestDerived> &aDest) const
  {
    eigen_assert(rows()==b.rows());

    Index rhsCols = b.cols();
    DestDerived& dest(aDest.derived());
    ComputationInfo global_info = Success;
    for(Index k=0; k<rhsCols; ++k)
    {
      typename DestDerived::ColXpr xk(dest,k);
      typename Rhs::ConstColXpr bk(b,k);
      derived()._solve_vector_with_guess_impl(bk,xk);

      if(m_info==NumericalIssue)
        global_info = NumericalIssue;
      else if(m_info==NoConvergence)
        global_info = NoConvergence;
    }
    m_info = global_info;
  }

  template<typename Rhs, typename DestDerived>
  typename internal::enable_if<Rhs::ColsAtCompileTime==1 || DestDerived::ColsAtCompileTime==1>::type
  _solve_with_guess_impl(const Rhs& b, MatrixBase<DestDerived> &dest) const
  {
    derived()._solve_vector_with_guess_impl(b,dest.derived());
  }

  template<typename Rhs,typename Dest>
  void _solve_impl(const Rhs& b, Dest& x) const
  {
    x.setZero();
    derived()._solve_with_guess_impl(b,x);
  }

protected:
  void init()
  {
    m_isInitialized = false;
    m_analysisIsOk = false;
    m_factorizationIsOk = false;
    m_maxIterations = -1;
    m_tolerance = NumTraits<Scalar>::epsilon();
  }

  typedef internal::generic_matrix_wrapper<MatrixType> MatrixWrapper;
  typedef typename MatrixWrapper::ActualMatrixType ActualMatrixType;

  const ActualMatrixType& matrix() const
  {
    return m_matrixWrapper.matrix();
  }

  template<typename InputType>
  void grab(const InputType &A)
  {
    m_matrixWrapper.grab(A);
  }

  MatrixWrapper m_matrixWrapper;
  Preconditioner m_preconditioner;

  Index m_maxIterations;
  RealScalar m_tolerance;

  mutable RealScalar m_error;
  mutable Index m_iterations;
  mutable ComputationInfo m_info;
  mutable bool m_analysisIsOk, m_factorizationIsOk;
};

} 

#endif 
