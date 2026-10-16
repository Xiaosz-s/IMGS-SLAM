#ifndef EIGEN_COMPLEX_SCHUR_H
#define EIGEN_COMPLEX_SCHUR_H

#include "./HessenbergDecomposition.h"

namespace Eigen { 

namespace internal {
template<typename MatrixType, bool IsComplex> struct complex_schur_reduce_to_hessenberg;
}

template<typename MatrixType_> class ComplexSchur
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

    typedef Matrix<ComplexScalar, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> ComplexMatrixType;

    explicit ComplexSchur(Index size = RowsAtCompileTime==Dynamic ? 1 : RowsAtCompileTime)
      : m_matT(size,size),
        m_matU(size,size),
        m_hess(size),
        m_isInitialized(false),
        m_matUisUptodate(false),
        m_maxIters(-1)
    {}

    template<typename InputType>
    explicit ComplexSchur(const EigenBase<InputType>& matrix, bool computeU = true)
      : m_matT(matrix.rows(),matrix.cols()),
        m_matU(matrix.rows(),matrix.cols()),
        m_hess(matrix.rows()),
        m_isInitialized(false),
        m_matUisUptodate(false),
        m_maxIters(-1)
    {
      compute(matrix.derived(), computeU);
    }

    const ComplexMatrixType& matrixU() const
    {
      eigen_assert(m_isInitialized && "ComplexSchur is not initialized.");
      eigen_assert(m_matUisUptodate && "The matrix U has not been computed during the ComplexSchur decomposition.");
      return m_matU;
    }

    const ComplexMatrixType& matrixT() const
    {
      eigen_assert(m_isInitialized && "ComplexSchur is not initialized.");
      return m_matT;
    }

    template<typename InputType>
    ComplexSchur& compute(const EigenBase<InputType>& matrix, bool computeU = true);

    template<typename HessMatrixType, typename OrthMatrixType>
    ComplexSchur& computeFromHessenberg(const HessMatrixType& matrixH, const OrthMatrixType& matrixQ,  bool computeU=true);

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "ComplexSchur is not initialized.");
      return m_info;
    }

    ComplexSchur& setMaxIterations(Index maxIters)
    {
      m_maxIters = maxIters;
      return *this;
    }

    Index getMaxIterations()
    {
      return m_maxIters;
    }

    static const int m_maxIterationsPerRow = 30;

  protected:
    ComplexMatrixType m_matT, m_matU;
    HessenbergDecomposition<MatrixType> m_hess;
    ComputationInfo m_info;
    bool m_isInitialized;
    bool m_matUisUptodate;
    Index m_maxIters;

  private:  
    bool subdiagonalEntryIsNeglegible(Index i);
    ComplexScalar computeShift(Index iu, Index iter);
    void reduceToTriangularForm(bool computeU);
    friend struct internal::complex_schur_reduce_to_hessenberg<MatrixType, NumTraits<Scalar>::IsComplex>;
};

template<typename MatrixType>
inline bool ComplexSchur<MatrixType>::subdiagonalEntryIsNeglegible(Index i)
{
  RealScalar d = numext::norm1(m_matT.coeff(i,i)) + numext::norm1(m_matT.coeff(i+1,i+1));
  RealScalar sd = numext::norm1(m_matT.coeff(i+1,i));
  if (internal::isMuchSmallerThan(sd, d, NumTraits<RealScalar>::epsilon()))
  {
    m_matT.coeffRef(i+1,i) = ComplexScalar(0);
    return true;
  }
  return false;
}

template<typename MatrixType>
typename ComplexSchur<MatrixType>::ComplexScalar ComplexSchur<MatrixType>::computeShift(Index iu, Index iter)
{
  using std::abs;
  if (iter == 10 || iter == 20) 
  {
    
    return abs(numext::real(m_matT.coeff(iu,iu-1))) + abs(numext::real(m_matT.coeff(iu-1,iu-2)));
  }

  Matrix<ComplexScalar,2,2> t = m_matT.template block<2,2>(iu-1,iu-1);
  RealScalar normt = t.cwiseAbs().sum();
  t /= normt;     

  ComplexScalar b = t.coeff(0,1) * t.coeff(1,0);
  ComplexScalar c = t.coeff(0,0) - t.coeff(1,1);
  ComplexScalar disc = sqrt(c*c + RealScalar(4)*b);
  ComplexScalar det = t.coeff(0,0) * t.coeff(1,1) - b;
  ComplexScalar trace = t.coeff(0,0) + t.coeff(1,1);
  ComplexScalar eival1 = (trace + disc) / RealScalar(2);
  ComplexScalar eival2 = (trace - disc) / RealScalar(2);
  RealScalar eival1_norm = numext::norm1(eival1);
  RealScalar eival2_norm = numext::norm1(eival2);

  if(eival1_norm > eival2_norm)
    eival2 = det / eival1;
  else if(eival2_norm!=RealScalar(0))
    eival1 = det / eival2;

  if(numext::norm1(eival1-t.coeff(1,1)) < numext::norm1(eival2-t.coeff(1,1)))
    return normt * eival1;
  else
    return normt * eival2;
}

template<typename MatrixType>
template<typename InputType>
ComplexSchur<MatrixType>& ComplexSchur<MatrixType>::compute(const EigenBase<InputType>& matrix, bool computeU)
{
  m_matUisUptodate = false;
  eigen_assert(matrix.cols() == matrix.rows());

  if(matrix.cols() == 1)
  {
    m_matT = matrix.derived().template cast<ComplexScalar>();
    if(computeU)  m_matU = ComplexMatrixType::Identity(1,1);
    m_info = Success;
    m_isInitialized = true;
    m_matUisUptodate = computeU;
    return *this;
  }

  internal::complex_schur_reduce_to_hessenberg<MatrixType, NumTraits<Scalar>::IsComplex>::run(*this, matrix.derived(), computeU);
  computeFromHessenberg(m_matT, m_matU, computeU);
  return *this;
}

template<typename MatrixType>
template<typename HessMatrixType, typename OrthMatrixType>
ComplexSchur<MatrixType>& ComplexSchur<MatrixType>::computeFromHessenberg(const HessMatrixType& matrixH, const OrthMatrixType& matrixQ, bool computeU)
{
  m_matT = matrixH;
  if(computeU)
    m_matU = matrixQ;
  reduceToTriangularForm(computeU);
  return *this;
}
namespace internal {

template<typename MatrixType, bool IsComplex>
struct complex_schur_reduce_to_hessenberg
{
  
  static void run(ComplexSchur<MatrixType>& _this, const MatrixType& matrix, bool computeU)
  {
    _this.m_hess.compute(matrix);
    _this.m_matT = _this.m_hess.matrixH();
    if(computeU)  _this.m_matU = _this.m_hess.matrixQ();
  }
};

template<typename MatrixType>
struct complex_schur_reduce_to_hessenberg<MatrixType, false>
{
  static void run(ComplexSchur<MatrixType>& _this, const MatrixType& matrix, bool computeU)
  {
    typedef typename ComplexSchur<MatrixType>::ComplexScalar ComplexScalar;

    _this.m_hess.compute(matrix);
    _this.m_matT = _this.m_hess.matrixH().template cast<ComplexScalar>();
    if(computeU)  
    {
      
      MatrixType Q = _this.m_hess.matrixQ(); 
      _this.m_matU = Q.template cast<ComplexScalar>();
    }
  }
};

} 

template<typename MatrixType>
void ComplexSchur<MatrixType>::reduceToTriangularForm(bool computeU)
{  
  Index maxIters = m_maxIters;
  if (maxIters == -1)
    maxIters = m_maxIterationsPerRow * m_matT.rows();

  Index iu = m_matT.cols() - 1;
  Index il;
  Index iter = 0; 
  Index totalIter = 0; 

  while(true)
  {
    
    while(iu > 0)
    {
      if(!subdiagonalEntryIsNeglegible(iu-1)) break;
      iter = 0;
      --iu;
    }

    if(iu==0) break;

    iter++;
    totalIter++;
    if(totalIter > maxIters) break;

    il = iu-1;
    while(il > 0 && !subdiagonalEntryIsNeglegible(il-1))
    {
      --il;
    }

    ComplexScalar shift = computeShift(iu, iter);
    JacobiRotation<ComplexScalar> rot;
    rot.makeGivens(m_matT.coeff(il,il) - shift, m_matT.coeff(il+1,il));
    m_matT.rightCols(m_matT.cols()-il).applyOnTheLeft(il, il+1, rot.adjoint());
    m_matT.topRows((std::min)(il+2,iu)+1).applyOnTheRight(il, il+1, rot);
    if(computeU) m_matU.applyOnTheRight(il, il+1, rot);

    for(Index i=il+1 ; i<iu ; i++)
    {
      rot.makeGivens(m_matT.coeffRef(i,i-1), m_matT.coeffRef(i+1,i-1), &m_matT.coeffRef(i,i-1));
      m_matT.coeffRef(i+1,i-1) = ComplexScalar(0);
      m_matT.rightCols(m_matT.cols()-i).applyOnTheLeft(i, i+1, rot.adjoint());
      m_matT.topRows((std::min)(i+2,iu)+1).applyOnTheRight(i, i+1, rot);
      if(computeU) m_matU.applyOnTheRight(i, i+1, rot);
    }
  }

  if(totalIter <= maxIters)
    m_info = Success;
  else
    m_info = NoConvergence;

  m_isInitialized = true;
  m_matUisUptodate = computeU;
}

} 

#endif 
