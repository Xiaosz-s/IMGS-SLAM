#ifndef EIGEN_REAL_SCHUR_H
#define EIGEN_REAL_SCHUR_H

#include "./HessenbergDecomposition.h"

namespace Eigen { 

template<typename MatrixType_> class RealSchur
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
    typedef std::complex<typename NumTraits<Scalar>::Real> ComplexScalar;
    typedef Eigen::Index Index; 

    typedef Matrix<ComplexScalar, ColsAtCompileTime, 1, Options & ~RowMajor, MaxColsAtCompileTime, 1> EigenvalueType;
    typedef Matrix<Scalar, ColsAtCompileTime, 1, Options & ~RowMajor, MaxColsAtCompileTime, 1> ColumnVectorType;

    explicit RealSchur(Index size = RowsAtCompileTime==Dynamic ? 1 : RowsAtCompileTime)
            : m_matT(size, size),
              m_matU(size, size),
              m_workspaceVector(size),
              m_hess(size),
              m_isInitialized(false),
              m_matUisUptodate(false),
              m_maxIters(-1)
    { }

    template<typename InputType>
    explicit RealSchur(const EigenBase<InputType>& matrix, bool computeU = true)
            : m_matT(matrix.rows(),matrix.cols()),
              m_matU(matrix.rows(),matrix.cols()),
              m_workspaceVector(matrix.rows()),
              m_hess(matrix.rows()),
              m_isInitialized(false),
              m_matUisUptodate(false),
              m_maxIters(-1)
    {
      compute(matrix.derived(), computeU);
    }

    const MatrixType& matrixU() const
    {
      eigen_assert(m_isInitialized && "RealSchur is not initialized.");
      eigen_assert(m_matUisUptodate && "The matrix U has not been computed during the RealSchur decomposition.");
      return m_matU;
    }

    const MatrixType& matrixT() const
    {
      eigen_assert(m_isInitialized && "RealSchur is not initialized.");
      return m_matT;
    }

    template<typename InputType>
    RealSchur& compute(const EigenBase<InputType>& matrix, bool computeU = true);

    template<typename HessMatrixType, typename OrthMatrixType>
    RealSchur& computeFromHessenberg(const HessMatrixType& matrixH, const OrthMatrixType& matrixQ,  bool computeU);

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "RealSchur is not initialized.");
      return m_info;
    }

    RealSchur& setMaxIterations(Index maxIters)
    {
      m_maxIters = maxIters;
      return *this;
    }

    Index getMaxIterations()
    {
      return m_maxIters;
    }

    static const int m_maxIterationsPerRow = 40;

  private:
    
    MatrixType m_matT;
    MatrixType m_matU;
    ColumnVectorType m_workspaceVector;
    HessenbergDecomposition<MatrixType> m_hess;
    ComputationInfo m_info;
    bool m_isInitialized;
    bool m_matUisUptodate;
    Index m_maxIters;

    typedef Matrix<Scalar,3,1> Vector3s;

    Scalar computeNormOfT();
    Index findSmallSubdiagEntry(Index iu, const Scalar& considerAsZero);
    void splitOffTwoRows(Index iu, bool computeU, const Scalar& exshift);
    void computeShift(Index iu, Index iter, Scalar& exshift, Vector3s& shiftInfo);
    void initFrancisQRStep(Index il, Index iu, const Vector3s& shiftInfo, Index& im, Vector3s& firstHouseholderVector);
    void performFrancisQRStep(Index il, Index im, Index iu, bool computeU, const Vector3s& firstHouseholderVector, Scalar* workspace);
};

template<typename MatrixType>
template<typename InputType>
RealSchur<MatrixType>& RealSchur<MatrixType>::compute(const EigenBase<InputType>& matrix, bool computeU)
{
  const Scalar considerAsZero = (std::numeric_limits<Scalar>::min)();

  eigen_assert(matrix.cols() == matrix.rows());
  Index maxIters = m_maxIters;
  if (maxIters == -1)
    maxIters = m_maxIterationsPerRow * matrix.rows();

  Scalar scale = matrix.derived().cwiseAbs().maxCoeff();
  if(scale<considerAsZero)
  {
    m_matT.setZero(matrix.rows(),matrix.cols());
    if(computeU)
      m_matU.setIdentity(matrix.rows(),matrix.cols());
    m_info = Success;
    m_isInitialized = true;
    m_matUisUptodate = computeU;
    return *this;
  }

  m_hess.compute(matrix.derived()/scale);

  m_workspaceVector.resize(matrix.cols());
  if(computeU)
    m_hess.matrixQ().evalTo(m_matU, m_workspaceVector);
  computeFromHessenberg(m_hess.matrixH(), m_matU, computeU);

  m_matT *= scale;
  
  return *this;
}
template<typename MatrixType>
template<typename HessMatrixType, typename OrthMatrixType>
RealSchur<MatrixType>& RealSchur<MatrixType>::computeFromHessenberg(const HessMatrixType& matrixH, const OrthMatrixType& matrixQ,  bool computeU)
{
  using std::abs;

  m_matT = matrixH;
  m_workspaceVector.resize(m_matT.cols());
  if(computeU && !internal::is_same_dense(m_matU,matrixQ))
    m_matU = matrixQ;
  
  Index maxIters = m_maxIters;
  if (maxIters == -1)
    maxIters = m_maxIterationsPerRow * matrixH.rows();
  Scalar* workspace = &m_workspaceVector.coeffRef(0);

  Index iu = m_matT.cols() - 1;
  Index iter = 0;      
  Index totalIter = 0; 
  Scalar exshift(0);   
  Scalar norm = computeNormOfT();

  Scalar considerAsZero = numext::maxi<Scalar>( norm * numext::abs2(NumTraits<Scalar>::epsilon()),
                                                (std::numeric_limits<Scalar>::min)() );

  if(norm!=Scalar(0))
  {
    while (iu >= 0)
    {
      Index il = findSmallSubdiagEntry(iu,considerAsZero);

      if (il == iu) 
      {
        m_matT.coeffRef(iu,iu) = m_matT.coeff(iu,iu) + exshift;
        if (iu > 0)
          m_matT.coeffRef(iu, iu-1) = Scalar(0);
        iu--;
        iter = 0;
      }
      else if (il == iu-1) 
      {
        splitOffTwoRows(iu, computeU, exshift);
        iu -= 2;
        iter = 0;
      }
      else 
      {
        
        Vector3s firstHouseholderVector = Vector3s::Zero(), shiftInfo;
        computeShift(iu, iter, exshift, shiftInfo);
        iter = iter + 1;
        totalIter = totalIter + 1;
        if (totalIter > maxIters) break;
        Index im;
        initFrancisQRStep(il, iu, shiftInfo, im, firstHouseholderVector);
        performFrancisQRStep(il, im, iu, computeU, firstHouseholderVector, workspace);
      }
    }
  }
  if(totalIter <= maxIters)
    m_info = Success;
  else
    m_info = NoConvergence;

  m_isInitialized = true;
  m_matUisUptodate = computeU;
  return *this;
}

template<typename MatrixType>
inline typename MatrixType::Scalar RealSchur<MatrixType>::computeNormOfT()
{
  const Index size = m_matT.cols();

  Scalar norm(0);
  for (Index j = 0; j < size; ++j)
    norm += m_matT.col(j).segment(0, (std::min)(size,j+2)).cwiseAbs().sum();
  return norm;
}

template<typename MatrixType>
inline Index RealSchur<MatrixType>::findSmallSubdiagEntry(Index iu, const Scalar& considerAsZero)
{
  using std::abs;
  Index res = iu;
  while (res > 0)
  {
    Scalar s = abs(m_matT.coeff(res-1,res-1)) + abs(m_matT.coeff(res,res));

    s = numext::maxi<Scalar>(s * NumTraits<Scalar>::epsilon(), considerAsZero);
    
    if (abs(m_matT.coeff(res,res-1)) <= s)
      break;
    res--;
  }
  return res;
}

template<typename MatrixType>
inline void RealSchur<MatrixType>::splitOffTwoRows(Index iu, bool computeU, const Scalar& exshift)
{
  using std::sqrt;
  using std::abs;
  const Index size = m_matT.cols();

  Scalar p = Scalar(0.5) * (m_matT.coeff(iu-1,iu-1) - m_matT.coeff(iu,iu));
  Scalar q = p * p + m_matT.coeff(iu,iu-1) * m_matT.coeff(iu-1,iu);   
  m_matT.coeffRef(iu,iu) += exshift;
  m_matT.coeffRef(iu-1,iu-1) += exshift;

  if (q >= Scalar(0)) 
  {
    Scalar z = sqrt(abs(q));
    JacobiRotation<Scalar> rot;
    if (p >= Scalar(0))
      rot.makeGivens(p + z, m_matT.coeff(iu, iu-1));
    else
      rot.makeGivens(p - z, m_matT.coeff(iu, iu-1));

    m_matT.rightCols(size-iu+1).applyOnTheLeft(iu-1, iu, rot.adjoint());
    m_matT.topRows(iu+1).applyOnTheRight(iu-1, iu, rot);
    m_matT.coeffRef(iu, iu-1) = Scalar(0); 
    if (computeU)
      m_matU.applyOnTheRight(iu-1, iu, rot);
  }

  if (iu > 1) 
    m_matT.coeffRef(iu-1, iu-2) = Scalar(0);
}

template<typename MatrixType>
inline void RealSchur<MatrixType>::computeShift(Index iu, Index iter, Scalar& exshift, Vector3s& shiftInfo)
{
  using std::sqrt;
  using std::abs;
  shiftInfo.coeffRef(0) = m_matT.coeff(iu,iu);
  shiftInfo.coeffRef(1) = m_matT.coeff(iu-1,iu-1);
  shiftInfo.coeffRef(2) = m_matT.coeff(iu,iu-1) * m_matT.coeff(iu-1,iu);

  if (iter == 10)
  {
    exshift += shiftInfo.coeff(0);
    for (Index i = 0; i <= iu; ++i)
      m_matT.coeffRef(i,i) -= shiftInfo.coeff(0);
    Scalar s = abs(m_matT.coeff(iu,iu-1)) + abs(m_matT.coeff(iu-1,iu-2));
    shiftInfo.coeffRef(0) = Scalar(0.75) * s;
    shiftInfo.coeffRef(1) = Scalar(0.75) * s;
    shiftInfo.coeffRef(2) = Scalar(-0.4375) * s * s;
  }

  if (iter == 30)
  {
    Scalar s = (shiftInfo.coeff(1) - shiftInfo.coeff(0)) / Scalar(2.0);
    s = s * s + shiftInfo.coeff(2);
    if (s > Scalar(0))
    {
      s = sqrt(s);
      if (shiftInfo.coeff(1) < shiftInfo.coeff(0))
        s = -s;
      s = s + (shiftInfo.coeff(1) - shiftInfo.coeff(0)) / Scalar(2.0);
      s = shiftInfo.coeff(0) - shiftInfo.coeff(2) / s;
      exshift += s;
      for (Index i = 0; i <= iu; ++i)
        m_matT.coeffRef(i,i) -= s;
      shiftInfo.setConstant(Scalar(0.964));
    }
  }
}

template<typename MatrixType>
inline void RealSchur<MatrixType>::initFrancisQRStep(Index il, Index iu, const Vector3s& shiftInfo, Index& im, Vector3s& firstHouseholderVector)
{
  using std::abs;
  Vector3s& v = firstHouseholderVector; 

  for (im = iu-2; im >= il; --im)
  {
    const Scalar Tmm = m_matT.coeff(im,im);
    const Scalar r = shiftInfo.coeff(0) - Tmm;
    const Scalar s = shiftInfo.coeff(1) - Tmm;
    v.coeffRef(0) = (r * s - shiftInfo.coeff(2)) / m_matT.coeff(im+1,im) + m_matT.coeff(im,im+1);
    v.coeffRef(1) = m_matT.coeff(im+1,im+1) - Tmm - r - s;
    v.coeffRef(2) = m_matT.coeff(im+2,im+1);
    if (im == il) {
      break;
    }
    const Scalar lhs = m_matT.coeff(im,im-1) * (abs(v.coeff(1)) + abs(v.coeff(2)));
    const Scalar rhs = v.coeff(0) * (abs(m_matT.coeff(im-1,im-1)) + abs(Tmm) + abs(m_matT.coeff(im+1,im+1)));
    if (abs(lhs) < NumTraits<Scalar>::epsilon() * rhs)
      break;
  }
}

template<typename MatrixType>
inline void RealSchur<MatrixType>::performFrancisQRStep(Index il, Index im, Index iu, bool computeU, const Vector3s& firstHouseholderVector, Scalar* workspace)
{
  eigen_assert(im >= il);
  eigen_assert(im <= iu-2);

  const Index size = m_matT.cols();

  for (Index k = im; k <= iu-2; ++k)
  {
    bool firstIteration = (k == im);

    Vector3s v;
    if (firstIteration)
      v = firstHouseholderVector;
    else
      v = m_matT.template block<3,1>(k,k-1);

    Scalar tau, beta;
    Matrix<Scalar, 2, 1> ess;
    v.makeHouseholder(ess, tau, beta);
    
    if (beta != Scalar(0)) 
    {
      if (firstIteration && k > il)
        m_matT.coeffRef(k,k-1) = -m_matT.coeff(k,k-1);
      else if (!firstIteration)
        m_matT.coeffRef(k,k-1) = beta;

      m_matT.block(k, k, 3, size-k).applyHouseholderOnTheLeft(ess, tau, workspace);
      m_matT.block(0, k, (std::min)(iu,k+3) + 1, 3).applyHouseholderOnTheRight(ess, tau, workspace);
      if (computeU)
        m_matU.block(0, k, size, 3).applyHouseholderOnTheRight(ess, tau, workspace);
    }
  }

  Matrix<Scalar, 2, 1> v = m_matT.template block<2,1>(iu-1, iu-2);
  Scalar tau, beta;
  Matrix<Scalar, 1, 1> ess;
  v.makeHouseholder(ess, tau, beta);

  if (beta != Scalar(0)) 
  {
    m_matT.coeffRef(iu-1, iu-2) = beta;
    m_matT.block(iu-1, iu-1, 2, size-iu+1).applyHouseholderOnTheLeft(ess, tau, workspace);
    m_matT.block(0, iu-1, iu+1, 2).applyHouseholderOnTheRight(ess, tau, workspace);
    if (computeU)
      m_matU.block(0, iu-1, size, 2).applyHouseholderOnTheRight(ess, tau, workspace);
  }

  for (Index i = im+2; i <= iu; ++i)
  {
    m_matT.coeffRef(i,i-2) = Scalar(0);
    if (i > im+2)
      m_matT.coeffRef(i,i-3) = Scalar(0);
  }
}

} 

#endif 
