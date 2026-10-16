#ifndef EIGEN_GENERALIZEDEIGENSOLVER_H
#define EIGEN_GENERALIZEDEIGENSOLVER_H

#include "./RealQZ.h"

namespace Eigen { 

template<typename MatrixType_> class GeneralizedEigenSolver
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

    typedef Matrix<Scalar, ColsAtCompileTime, 1, Options & ~RowMajor, MaxColsAtCompileTime, 1> VectorType;

    typedef Matrix<ComplexScalar, ColsAtCompileTime, 1, Options & ~RowMajor, MaxColsAtCompileTime, 1> ComplexVectorType;

    typedef CwiseBinaryOp<internal::scalar_quotient_op<ComplexScalar,Scalar>,ComplexVectorType,VectorType> EigenvalueType;

    typedef Matrix<ComplexScalar, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> EigenvectorsType;

    GeneralizedEigenSolver()
      : m_eivec(),
        m_alphas(),
        m_betas(),
        m_valuesOkay(false),
        m_vectorsOkay(false),
        m_realQZ()
    {}

    explicit GeneralizedEigenSolver(Index size)
      : m_eivec(size, size),
        m_alphas(size),
        m_betas(size),
        m_valuesOkay(false),
        m_vectorsOkay(false),
        m_realQZ(size),
        m_tmp(size)
    {}

    GeneralizedEigenSolver(const MatrixType& A, const MatrixType& B, bool computeEigenvectors = true)
      : m_eivec(A.rows(), A.cols()),
        m_alphas(A.cols()),
        m_betas(A.cols()),
        m_valuesOkay(false),
        m_vectorsOkay(false),
        m_realQZ(A.cols()),
        m_tmp(A.cols())
    {
      compute(A, B, computeEigenvectors);
    }

    EigenvectorsType eigenvectors() const {
      eigen_assert(m_vectorsOkay && "Eigenvectors for GeneralizedEigenSolver were not calculated.");
      return m_eivec;
    }

    EigenvalueType eigenvalues() const
    {
      eigen_assert(m_valuesOkay && "GeneralizedEigenSolver is not initialized.");
      return EigenvalueType(m_alphas,m_betas);
    }

    ComplexVectorType alphas() const
    {
      eigen_assert(m_valuesOkay && "GeneralizedEigenSolver is not initialized.");
      return m_alphas;
    }

    VectorType betas() const
    {
      eigen_assert(m_valuesOkay && "GeneralizedEigenSolver is not initialized.");
      return m_betas;
    }

    GeneralizedEigenSolver& compute(const MatrixType& A, const MatrixType& B, bool computeEigenvectors = true);

    ComputationInfo info() const
    {
      eigen_assert(m_valuesOkay && "EigenSolver is not initialized.");
      return m_realQZ.info();
    }

    GeneralizedEigenSolver& setMaxIterations(Index maxIters)
    {
      m_realQZ.setMaxIterations(maxIters);
      return *this;
    }

  protected:
    
    static void check_template_parameters()
    {
      EIGEN_STATIC_ASSERT_NON_INTEGER(Scalar);
      EIGEN_STATIC_ASSERT(!NumTraits<Scalar>::IsComplex, NUMERIC_TYPE_MUST_BE_REAL);
    }
    
    EigenvectorsType m_eivec;
    ComplexVectorType m_alphas;
    VectorType m_betas;
    bool m_valuesOkay, m_vectorsOkay;
    RealQZ<MatrixType> m_realQZ;
    ComplexVectorType m_tmp;
};

template<typename MatrixType>
GeneralizedEigenSolver<MatrixType>&
GeneralizedEigenSolver<MatrixType>::compute(const MatrixType& A, const MatrixType& B, bool computeEigenvectors)
{
  check_template_parameters();
  
  using std::sqrt;
  using std::abs;
  eigen_assert(A.cols() == A.rows() && B.cols() == A.rows() && B.cols() == B.rows());
  Index size = A.cols();
  m_valuesOkay = false;
  m_vectorsOkay = false;

  m_realQZ.compute(A, B, computeEigenvectors);
  if (m_realQZ.info() == Success)
  {
    
    m_alphas.resize(size);
    m_betas.resize(size);
    if (computeEigenvectors)
    {
      m_eivec.resize(size,size);
      m_tmp.resize(size);
    }

    Map<VectorType> v(reinterpret_cast<Scalar*>(m_tmp.data()), size);
    ComplexVectorType &cv = m_tmp;
    const MatrixType &mS = m_realQZ.matrixS();
    const MatrixType &mT = m_realQZ.matrixT();

    Index i = 0;
    while (i < size)
    {
      if (i == size - 1 || mS.coeff(i+1, i) == Scalar(0))
      {
        
        m_alphas.coeffRef(i) = mS.diagonal().coeff(i);
        m_betas.coeffRef(i)  = mT.diagonal().coeff(i);
        if (computeEigenvectors)
        {
          v.setConstant(Scalar(0.0));
          v.coeffRef(i) = Scalar(1.0);
          
          if(abs(m_betas.coeffRef(i)) >= (std::numeric_limits<RealScalar>::min)())
          {
            
            const Scalar alpha = real(m_alphas.coeffRef(i));
            const Scalar beta = m_betas.coeffRef(i);
            for (Index j = i-1; j >= 0; j--)
            {
              const Index st = j+1;
              const Index sz = i-j;
              if (j > 0 && mS.coeff(j, j-1) != Scalar(0))
              {
                
                Matrix<Scalar, 2, 1> rhs = (alpha*mT.template block<2,Dynamic>(j-1,st,2,sz) - beta*mS.template block<2,Dynamic>(j-1,st,2,sz)) .lazyProduct( v.segment(st,sz) );
                Matrix<Scalar, 2, 2> lhs = beta * mS.template block<2,2>(j-1,j-1) - alpha * mT.template block<2,2>(j-1,j-1);
                v.template segment<2>(j-1) = lhs.partialPivLu().solve(rhs);
                j--;
              }
              else
              {
                v.coeffRef(j) = -v.segment(st,sz).transpose().cwiseProduct(beta*mS.block(j,st,1,sz) - alpha*mT.block(j,st,1,sz)).sum() / (beta*mS.coeffRef(j,j) - alpha*mT.coeffRef(j,j));
              }
            }
          }
          m_eivec.col(i).real().noalias() = m_realQZ.matrixZ().transpose() * v;
          m_eivec.col(i).real().normalize();
          m_eivec.col(i).imag().setConstant(0);
        }
        ++i;
      }
      else
      {

        RealScalar a = mT.diagonal().coeff(i),
                   b = mT.diagonal().coeff(i+1);
        const RealScalar beta = m_betas.coeffRef(i) = m_betas.coeffRef(i+1) = a*b;

        Matrix<RealScalar,2,2> S2 = mS.template block<2,2>(i,i) * Matrix<Scalar,2,1>(b,a).asDiagonal();

        Scalar p = Scalar(0.5) * (S2.coeff(0,0) - S2.coeff(1,1));
        Scalar z = sqrt(abs(p * p + S2.coeff(1,0) * S2.coeff(0,1)));
        const ComplexScalar alpha = ComplexScalar(S2.coeff(1,1) + p, (beta > 0) ? z : -z);
        m_alphas.coeffRef(i)   = conj(alpha);
        m_alphas.coeffRef(i+1) = alpha;

        if (computeEigenvectors) {
          
          cv.setZero();
          cv.coeffRef(i+1) = Scalar(1.0);
          
          cv.coeffRef(i) = -(static_cast<Scalar>(beta*mS.coeffRef(i,i+1)) - alpha*mT.coeffRef(i,i+1))
                          / (static_cast<Scalar>(beta*mS.coeffRef(i,i))   - alpha*mT.coeffRef(i,i));
          for (Index j = i-1; j >= 0; j--)
          {
            const Index st = j+1;
            const Index sz = i+1-j;
            if (j > 0 && mS.coeff(j, j-1) != Scalar(0))
            {
              
              Matrix<ComplexScalar, 2, 1> rhs = (alpha*mT.template block<2,Dynamic>(j-1,st,2,sz) - beta*mS.template block<2,Dynamic>(j-1,st,2,sz)) .lazyProduct( cv.segment(st,sz) );
              Matrix<ComplexScalar, 2, 2> lhs = beta * mS.template block<2,2>(j-1,j-1) - alpha * mT.template block<2,2>(j-1,j-1);
              cv.template segment<2>(j-1) = lhs.partialPivLu().solve(rhs);
              j--;
            } else {
              cv.coeffRef(j) =  cv.segment(st,sz).transpose().cwiseProduct(beta*mS.block(j,st,1,sz) - alpha*mT.block(j,st,1,sz)).sum()
                              / (alpha*mT.coeffRef(j,j) - static_cast<Scalar>(beta*mS.coeffRef(j,j)));
            }
          }
          m_eivec.col(i+1).noalias() = (m_realQZ.matrixZ().transpose() * cv);
          m_eivec.col(i+1).normalize();
          m_eivec.col(i) = m_eivec.col(i+1).conjugate();
        }
        i += 2;
      }
    }

    m_valuesOkay = true;
    m_vectorsOkay = computeEigenvectors;
  }
  return *this;
}

} 

#endif 
