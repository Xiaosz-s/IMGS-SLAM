#ifndef EIGEN_COMPLEX_EIGEN_SOLVER_H
#define EIGEN_COMPLEX_EIGEN_SOLVER_H

#include "./ComplexSchur.h"

namespace Eigen { 

template<typename MatrixType_> class ComplexEigenSolver
{
  public:

    typedef MatrixType_ MatrixType;

    enum {
      RowsAtCompileTime = MatrixType::RowsAtCompileTime,
      ColsAtCompileTime = MatrixType::ColsAtCompileTime,
      Options = MatrixType::Options,
      MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };

    typedef typename MatrixType::Scalar Scalar;
    typedef typename NumTraits<Scalar>::Real RealScalar;
    typedef Eigen::Index Index; 

    typedef std::complex<RealScalar> ComplexScalar;

    typedef Matrix<ComplexScalar, ColsAtCompileTime, 1, Options&(~RowMajor), MaxColsAtCompileTime, 1> EigenvalueType;

    typedef Matrix<ComplexScalar, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> EigenvectorType;

    ComplexEigenSolver()
            : m_eivec(),
              m_eivalues(),
              m_schur(),
              m_isInitialized(false),
              m_eigenvectorsOk(false),
              m_matX()
    {}

    explicit ComplexEigenSolver(Index size)
            : m_eivec(size, size),
              m_eivalues(size),
              m_schur(size),
              m_isInitialized(false),
              m_eigenvectorsOk(false),
              m_matX(size, size)
    {}

    template<typename InputType>
    explicit ComplexEigenSolver(const EigenBase<InputType>& matrix, bool computeEigenvectors = true)
            : m_eivec(matrix.rows(),matrix.cols()),
              m_eivalues(matrix.cols()),
              m_schur(matrix.rows()),
              m_isInitialized(false),
              m_eigenvectorsOk(false),
              m_matX(matrix.rows(),matrix.cols())
    {
      compute(matrix.derived(), computeEigenvectors);
    }

    const EigenvectorType& eigenvectors() const
    {
      eigen_assert(m_isInitialized && "ComplexEigenSolver is not initialized.");
      eigen_assert(m_eigenvectorsOk && "The eigenvectors have not been computed together with the eigenvalues.");
      return m_eivec;
    }

    const EigenvalueType& eigenvalues() const
    {
      eigen_assert(m_isInitialized && "ComplexEigenSolver is not initialized.");
      return m_eivalues;
    }

    template<typename InputType>
    ComplexEigenSolver& compute(const EigenBase<InputType>& matrix, bool computeEigenvectors = true);

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "ComplexEigenSolver is not initialized.");
      return m_schur.info();
    }

    ComplexEigenSolver& setMaxIterations(Index maxIters)
    {
      m_schur.setMaxIterations(maxIters);
      return *this;
    }

    Index getMaxIterations()
    {
      return m_schur.getMaxIterations();
    }

  protected:
    
    static void check_template_parameters()
    {
      EIGEN_STATIC_ASSERT_NON_INTEGER(Scalar);
    }
    
    EigenvectorType m_eivec;
    EigenvalueType m_eivalues;
    ComplexSchur<MatrixType> m_schur;
    bool m_isInitialized;
    bool m_eigenvectorsOk;
    EigenvectorType m_matX;

  private:
    void doComputeEigenvectors(RealScalar matrixnorm);
    void sortEigenvalues(bool computeEigenvectors);
};

template<typename MatrixType>
template<typename InputType>
ComplexEigenSolver<MatrixType>& 
ComplexEigenSolver<MatrixType>::compute(const EigenBase<InputType>& matrix, bool computeEigenvectors)
{
  check_template_parameters();

  eigen_assert(matrix.cols() == matrix.rows());

  m_schur.compute(matrix.derived(), computeEigenvectors);

  if(m_schur.info() == Success)
  {
    m_eivalues = m_schur.matrixT().diagonal();
    if(computeEigenvectors)
      doComputeEigenvectors(m_schur.matrixT().norm());
    sortEigenvalues(computeEigenvectors);
  }

  m_isInitialized = true;
  m_eigenvectorsOk = computeEigenvectors;
  return *this;
}

template<typename MatrixType>
void ComplexEigenSolver<MatrixType>::doComputeEigenvectors(RealScalar matrixnorm)
{
  const Index n = m_eivalues.size();

  matrixnorm = numext::maxi(matrixnorm,(std::numeric_limits<RealScalar>::min)());

  m_matX = EigenvectorType::Zero(n, n);
  for(Index k=n-1 ; k>=0 ; k--)
  {
    m_matX.coeffRef(k,k) = ComplexScalar(1.0,0.0);
    
    for(Index i=k-1 ; i>=0 ; i--)
    {
      m_matX.coeffRef(i,k) = -m_schur.matrixT().coeff(i,k);
      if(k-i-1>0)
        m_matX.coeffRef(i,k) -= (m_schur.matrixT().row(i).segment(i+1,k-i-1) * m_matX.col(k).segment(i+1,k-i-1)).value();
      ComplexScalar z = m_schur.matrixT().coeff(i,i) - m_schur.matrixT().coeff(k,k);
      if(z==ComplexScalar(0))
      {

        numext::real_ref(z) = NumTraits<RealScalar>::epsilon() * matrixnorm;
      }
      m_matX.coeffRef(i,k) = m_matX.coeff(i,k) / z;
    }
  }

  m_eivec.noalias() = m_schur.matrixU() * m_matX;
  
  for(Index k=0 ; k<n ; k++)
  {
    m_eivec.col(k).normalize();
  }
}

template<typename MatrixType>
void ComplexEigenSolver<MatrixType>::sortEigenvalues(bool computeEigenvectors)
{
  const Index n =  m_eivalues.size();
  for (Index i=0; i<n; i++)
  {
    Index k;
    m_eivalues.cwiseAbs().tail(n-i).minCoeff(&k);
    if (k != 0)
    {
      k += i;
      std::swap(m_eivalues[k],m_eivalues[i]);
      if(computeEigenvectors)
	m_eivec.col(i).swap(m_eivec.col(k));
    }
  }
}

} 

#endif 
