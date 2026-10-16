#ifndef EIGEN_SELFADJOINTEIGENSOLVER_H
#define EIGEN_SELFADJOINTEIGENSOLVER_H

#include "./Tridiagonalization.h"

namespace Eigen { 

template<typename MatrixType_>
class GeneralizedSelfAdjointEigenSolver;

namespace internal {
template<typename SolverType,int Size,bool IsComplex> struct direct_selfadjoint_eigenvalues;

template<typename MatrixType, typename DiagType, typename SubDiagType>
EIGEN_DEVICE_FUNC
ComputationInfo computeFromTridiagonal_impl(DiagType& diag, SubDiagType& subdiag, const Index maxIterations, bool computeEigenvectors, MatrixType& eivec);
}

template<typename MatrixType_> class SelfAdjointEigenSolver
{
  public:

    typedef MatrixType_ MatrixType;
    enum {
      Size = MatrixType::RowsAtCompileTime,
      ColsAtCompileTime = MatrixType::ColsAtCompileTime,
      Options = MatrixType::Options,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };

    typedef typename MatrixType::Scalar Scalar;
    typedef Eigen::Index Index; 
    
    typedef Matrix<Scalar,Size,Size,ColMajor,MaxColsAtCompileTime,MaxColsAtCompileTime> EigenvectorsType;

    typedef typename NumTraits<Scalar>::Real RealScalar;
    
    friend struct internal::direct_selfadjoint_eigenvalues<SelfAdjointEigenSolver,Size,NumTraits<Scalar>::IsComplex>;

    typedef typename internal::plain_col_type<MatrixType, RealScalar>::type RealVectorType;
    typedef Tridiagonalization<MatrixType> TridiagonalizationType;
    typedef typename TridiagonalizationType::SubDiagonalType SubDiagonalType;

    EIGEN_DEVICE_FUNC
    SelfAdjointEigenSolver()
        : m_eivec(),
          m_eivalues(),
          m_subdiag(),
          m_hcoeffs(),
          m_info(InvalidInput),
          m_isInitialized(false),
          m_eigenvectorsOk(false)
    { }

    EIGEN_DEVICE_FUNC
    explicit SelfAdjointEigenSolver(Index size)
        : m_eivec(size, size),
          m_eivalues(size),
          m_subdiag(size > 1 ? size - 1 : 1),
          m_hcoeffs(size > 1 ? size - 1 : 1),
          m_isInitialized(false),
          m_eigenvectorsOk(false)
    {}

    template<typename InputType>
    EIGEN_DEVICE_FUNC
    explicit SelfAdjointEigenSolver(const EigenBase<InputType>& matrix, int options = ComputeEigenvectors)
      : m_eivec(matrix.rows(), matrix.cols()),
        m_eivalues(matrix.cols()),
        m_subdiag(matrix.rows() > 1 ? matrix.rows() - 1 : 1),
        m_hcoeffs(matrix.cols() > 1 ? matrix.cols() - 1 : 1),
        m_isInitialized(false),
        m_eigenvectorsOk(false)
    {
      compute(matrix.derived(), options);
    }

    template<typename InputType>
    EIGEN_DEVICE_FUNC
    SelfAdjointEigenSolver& compute(const EigenBase<InputType>& matrix, int options = ComputeEigenvectors);

    EIGEN_DEVICE_FUNC
    SelfAdjointEigenSolver& computeDirect(const MatrixType& matrix, int options = ComputeEigenvectors);

    SelfAdjointEigenSolver& computeFromTridiagonal(const RealVectorType& diag, const SubDiagonalType& subdiag , int options=ComputeEigenvectors);

    EIGEN_DEVICE_FUNC
    const EigenvectorsType& eigenvectors() const
    {
      eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
      eigen_assert(m_eigenvectorsOk && "The eigenvectors have not been computed together with the eigenvalues.");
      return m_eivec;
    }

    EIGEN_DEVICE_FUNC
    const RealVectorType& eigenvalues() const
    {
      eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
      return m_eivalues;
    }

    EIGEN_DEVICE_FUNC
    MatrixType operatorSqrt() const
    {
      eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
      eigen_assert(m_eigenvectorsOk && "The eigenvectors have not been computed together with the eigenvalues.");
      return m_eivec * m_eivalues.cwiseSqrt().asDiagonal() * m_eivec.adjoint();
    }

    EIGEN_DEVICE_FUNC
    MatrixType operatorInverseSqrt() const
    {
      eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
      eigen_assert(m_eigenvectorsOk && "The eigenvectors have not been computed together with the eigenvalues.");
      return m_eivec * m_eivalues.cwiseInverse().cwiseSqrt().asDiagonal() * m_eivec.adjoint();
    }

    EIGEN_DEVICE_FUNC
    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
      return m_info;
    }

    static const int m_maxIterations = 30;

  protected:
    static EIGEN_DEVICE_FUNC
    void check_template_parameters()
    {
      EIGEN_STATIC_ASSERT_NON_INTEGER(Scalar);
    }
    
    EigenvectorsType m_eivec;
    RealVectorType m_eivalues;
    typename TridiagonalizationType::SubDiagonalType m_subdiag;
    typename TridiagonalizationType::CoeffVectorType m_hcoeffs;
    ComputationInfo m_info;
    bool m_isInitialized;
    bool m_eigenvectorsOk;
};

namespace internal {

template<int StorageOrder,typename RealScalar, typename Scalar, typename Index>
EIGEN_DEVICE_FUNC
static void tridiagonal_qr_step(RealScalar* diag, RealScalar* subdiag, Index start, Index end, Scalar* matrixQ, Index n);
}

template<typename MatrixType>
template<typename InputType>
EIGEN_DEVICE_FUNC
SelfAdjointEigenSolver<MatrixType>& SelfAdjointEigenSolver<MatrixType>
::compute(const EigenBase<InputType>& a_matrix, int options)
{
  check_template_parameters();
  
  const InputType &matrix(a_matrix.derived());
  
  EIGEN_USING_STD(abs);
  eigen_assert(matrix.cols() == matrix.rows());
  eigen_assert((options&~(EigVecMask|GenEigMask))==0
          && (options&EigVecMask)!=EigVecMask
          && "invalid option parameter");
  bool computeEigenvectors = (options&ComputeEigenvectors)==ComputeEigenvectors;
  Index n = matrix.cols();
  m_eivalues.resize(n,1);

  if(n==1)
  {
    m_eivec = matrix;
    m_eivalues.coeffRef(0,0) = numext::real(m_eivec.coeff(0,0));
    if(computeEigenvectors)
      m_eivec.setOnes(n,n);
    m_info = Success;
    m_isInitialized = true;
    m_eigenvectorsOk = computeEigenvectors;
    return *this;
  }

  RealVectorType& diag = m_eivalues;
  EigenvectorsType& mat = m_eivec;

  mat = matrix.template triangularView<Lower>();
  RealScalar scale = mat.cwiseAbs().maxCoeff();
  if(scale==RealScalar(0)) scale = RealScalar(1);
  mat.template triangularView<Lower>() /= scale;
  m_subdiag.resize(n-1);
  m_hcoeffs.resize(n-1);
  internal::tridiagonalization_inplace(mat, diag, m_subdiag, m_hcoeffs, computeEigenvectors);

  m_info = internal::computeFromTridiagonal_impl(diag, m_subdiag, m_maxIterations, computeEigenvectors, m_eivec);

  m_eivalues *= scale;

  m_isInitialized = true;
  m_eigenvectorsOk = computeEigenvectors;
  return *this;
}

template<typename MatrixType>
SelfAdjointEigenSolver<MatrixType>& SelfAdjointEigenSolver<MatrixType>
::computeFromTridiagonal(const RealVectorType& diag, const SubDiagonalType& subdiag , int options)
{
  
  bool computeEigenvectors = (options&ComputeEigenvectors)==ComputeEigenvectors;

  m_eivalues = diag;
  m_subdiag = subdiag;
  if (computeEigenvectors)
  {
    m_eivec.setIdentity(diag.size(), diag.size());
  }
  m_info = internal::computeFromTridiagonal_impl(m_eivalues, m_subdiag, m_maxIterations, computeEigenvectors, m_eivec);

  m_isInitialized = true;
  m_eigenvectorsOk = computeEigenvectors;
  return *this;
}

namespace internal {

template<typename MatrixType, typename DiagType, typename SubDiagType>
EIGEN_DEVICE_FUNC
ComputationInfo computeFromTridiagonal_impl(DiagType& diag, SubDiagType& subdiag, const Index maxIterations, bool computeEigenvectors, MatrixType& eivec)
{
  ComputationInfo info;
  typedef typename MatrixType::Scalar Scalar;

  Index n = diag.size();
  Index end = n-1;
  Index start = 0;
  Index iter = 0; 
  
  typedef typename DiagType::RealScalar RealScalar;
  const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();
  const RealScalar precision_inv = RealScalar(1)/NumTraits<RealScalar>::epsilon();
  while (end>0)
  {
    for (Index i = start; i<end; ++i) {
      if (numext::abs(subdiag[i]) < considerAsZero) {
        subdiag[i] = RealScalar(0);
      } else {

        const RealScalar scaled_subdiag = precision_inv * subdiag[i];
        if (scaled_subdiag * scaled_subdiag <= (numext::abs(diag[i])+numext::abs(diag[i+1]))) {
          subdiag[i] = RealScalar(0);
        }
      }
    }

    while (end>0 && subdiag[end-1]==RealScalar(0))
    {
      end--;
    }
    if (end<=0)
      break;

    iter++;
    if(iter > maxIterations * n) break;

    start = end - 1;
    while (start>0 && subdiag[start-1]!=0)
      start--;

    internal::tridiagonal_qr_step<MatrixType::Flags&RowMajorBit ? RowMajor : ColMajor>(diag.data(), subdiag.data(), start, end, computeEigenvectors ? eivec.data() : (Scalar*)0, n);
  }
  if (iter <= maxIterations * n)
    info = Success;
  else
    info = NoConvergence;

  if (info == Success)
  {
    for (Index i = 0; i < n-1; ++i)
    {
      Index k;
      diag.segment(i,n-i).minCoeff(&k);
      if (k > 0)
      {
        numext::swap(diag[i], diag[k+i]);
        if(computeEigenvectors)
          eivec.col(i).swap(eivec.col(k+i));
      }
    }
  }
  return info;
}
  
template<typename SolverType,int Size,bool IsComplex> struct direct_selfadjoint_eigenvalues
{
  EIGEN_DEVICE_FUNC
  static inline void run(SolverType& eig, const typename SolverType::MatrixType& A, int options)
  { eig.compute(A,options); }
};

template<typename SolverType> struct direct_selfadjoint_eigenvalues<SolverType,3,false>
{
  typedef typename SolverType::MatrixType MatrixType;
  typedef typename SolverType::RealVectorType VectorType;
  typedef typename SolverType::Scalar Scalar;
  typedef typename SolverType::EigenvectorsType EigenvectorsType;

  EIGEN_DEVICE_FUNC
  static inline void computeRoots(const MatrixType& m, VectorType& roots)
  {
    EIGEN_USING_STD(sqrt)
    EIGEN_USING_STD(atan2)
    EIGEN_USING_STD(cos)
    EIGEN_USING_STD(sin)
    const Scalar s_inv3 = Scalar(1)/Scalar(3);
    const Scalar s_sqrt3 = sqrt(Scalar(3));

    Scalar c0 = m(0,0)*m(1,1)*m(2,2) + Scalar(2)*m(1,0)*m(2,0)*m(2,1) - m(0,0)*m(2,1)*m(2,1) - m(1,1)*m(2,0)*m(2,0) - m(2,2)*m(1,0)*m(1,0);
    Scalar c1 = m(0,0)*m(1,1) - m(1,0)*m(1,0) + m(0,0)*m(2,2) - m(2,0)*m(2,0) + m(1,1)*m(2,2) - m(2,1)*m(2,1);
    Scalar c2 = m(0,0) + m(1,1) + m(2,2);

    Scalar c2_over_3 = c2*s_inv3;
    Scalar a_over_3 = (c2*c2_over_3 - c1)*s_inv3;
    a_over_3 = numext::maxi(a_over_3, Scalar(0));

    Scalar half_b = Scalar(0.5)*(c0 + c2_over_3*(Scalar(2)*c2_over_3*c2_over_3 - c1));

    Scalar q = a_over_3*a_over_3*a_over_3 - half_b*half_b;
    q = numext::maxi(q, Scalar(0));

    Scalar rho = sqrt(a_over_3);
    Scalar theta = atan2(sqrt(q),half_b)*s_inv3;  
    Scalar cos_theta = cos(theta);
    Scalar sin_theta = sin(theta);
    
    roots(0) = c2_over_3 - rho*(cos_theta + s_sqrt3*sin_theta); 
    roots(1) = c2_over_3 - rho*(cos_theta - s_sqrt3*sin_theta); 
    roots(2) = c2_over_3 + Scalar(2)*rho*cos_theta;
  }

  EIGEN_DEVICE_FUNC
  static inline bool extract_kernel(MatrixType& mat, Ref<VectorType> res, Ref<VectorType> representative)
  {
    EIGEN_USING_STD(abs);
    EIGEN_USING_STD(sqrt);
    Index i0;
    
    mat.diagonal().cwiseAbs().maxCoeff(&i0);

    representative = mat.col(i0);
    Scalar n0, n1;
    VectorType c0, c1;
    n0 = (c0 = representative.cross(mat.col((i0+1)%3))).squaredNorm();
    n1 = (c1 = representative.cross(mat.col((i0+2)%3))).squaredNorm();
    if(n0>n1) res = c0/sqrt(n0);
    else      res = c1/sqrt(n1);

    return true;
  }

  EIGEN_DEVICE_FUNC
  static inline void run(SolverType& solver, const MatrixType& mat, int options)
  {
    eigen_assert(mat.cols() == 3 && mat.cols() == mat.rows());
    eigen_assert((options&~(EigVecMask|GenEigMask))==0
            && (options&EigVecMask)!=EigVecMask
            && "invalid option parameter");
    bool computeEigenvectors = (options&ComputeEigenvectors)==ComputeEigenvectors;
    
    EigenvectorsType& eivecs = solver.m_eivec;
    VectorType& eivals = solver.m_eivalues;

    Scalar shift = mat.trace() / Scalar(3);
    
    MatrixType scaledMat = mat.template selfadjointView<Lower>();
    scaledMat.diagonal().array() -= shift;
    Scalar scale = scaledMat.cwiseAbs().maxCoeff();
    if(scale > 0) scaledMat /= scale;   

    computeRoots(scaledMat,eivals);

    if(computeEigenvectors)
    {
      if((eivals(2)-eivals(0))<=Eigen::NumTraits<Scalar>::epsilon())
      {
        
        eivecs.setIdentity();
      }
      else
      {
        MatrixType tmp;
        tmp = scaledMat;

        Scalar d0 = eivals(2) - eivals(1);
        Scalar d1 = eivals(1) - eivals(0);
        Index k(0), l(2);
        if(d0 > d1)
        {
          numext::swap(k,l);
          d0 = d1;
        }

        {
          tmp.diagonal().array () -= eivals(k);
          
          extract_kernel(tmp, eivecs.col(k), eivecs.col(l));
        }

        if(d0<=2*Eigen::NumTraits<Scalar>::epsilon()*d1)
        {

          eivecs.col(l) -= eivecs.col(k).dot(eivecs.col(l))*eivecs.col(l);
          eivecs.col(l).normalize();
        }
        else
        {
          tmp = scaledMat;
          tmp.diagonal().array () -= eivals(l);

          VectorType dummy;
          extract_kernel(tmp, eivecs.col(l), dummy);
        }

        eivecs.col(1) = eivecs.col(2).cross(eivecs.col(0)).normalized();
      }
    }

    eivals *= scale;
    eivals.array() += shift;
    
    solver.m_info = Success;
    solver.m_isInitialized = true;
    solver.m_eigenvectorsOk = computeEigenvectors;
  }
};

template<typename SolverType> 
struct direct_selfadjoint_eigenvalues<SolverType,2,false>
{
  typedef typename SolverType::MatrixType MatrixType;
  typedef typename SolverType::RealVectorType VectorType;
  typedef typename SolverType::Scalar Scalar;
  typedef typename SolverType::EigenvectorsType EigenvectorsType;
  
  EIGEN_DEVICE_FUNC
  static inline void computeRoots(const MatrixType& m, VectorType& roots)
  {
    EIGEN_USING_STD(sqrt);
    const Scalar t0 = Scalar(0.5) * sqrt( numext::abs2(m(0,0)-m(1,1)) + Scalar(4)*numext::abs2(m(1,0)));
    const Scalar t1 = Scalar(0.5) * (m(0,0) + m(1,1));
    roots(0) = t1 - t0;
    roots(1) = t1 + t0;
  }
  
  EIGEN_DEVICE_FUNC
  static inline void run(SolverType& solver, const MatrixType& mat, int options)
  {
    EIGEN_USING_STD(sqrt);
    EIGEN_USING_STD(abs);
    
    eigen_assert(mat.cols() == 2 && mat.cols() == mat.rows());
    eigen_assert((options&~(EigVecMask|GenEigMask))==0
            && (options&EigVecMask)!=EigVecMask
            && "invalid option parameter");
    bool computeEigenvectors = (options&ComputeEigenvectors)==ComputeEigenvectors;
    
    EigenvectorsType& eivecs = solver.m_eivec;
    VectorType& eivals = solver.m_eivalues;

    Scalar shift = mat.trace() / Scalar(2);
    MatrixType scaledMat = mat;
    scaledMat.coeffRef(0,1) = mat.coeff(1,0);
    scaledMat.diagonal().array() -= shift;
    Scalar scale = scaledMat.cwiseAbs().maxCoeff();
    if(scale > Scalar(0))
      scaledMat /= scale;

    computeRoots(scaledMat,eivals);

    if(computeEigenvectors)
    {
      if((eivals(1)-eivals(0))<=abs(eivals(1))*Eigen::NumTraits<Scalar>::epsilon())
      {
        eivecs.setIdentity();
      }
      else
      {
        scaledMat.diagonal().array () -= eivals(1);
        Scalar a2 = numext::abs2(scaledMat(0,0));
        Scalar c2 = numext::abs2(scaledMat(1,1));
        Scalar b2 = numext::abs2(scaledMat(1,0));
        if(a2>c2)
        {
          eivecs.col(1) << -scaledMat(1,0), scaledMat(0,0);
          eivecs.col(1) /= sqrt(a2+b2);
        }
        else
        {
          eivecs.col(1) << -scaledMat(1,1), scaledMat(1,0);
          eivecs.col(1) /= sqrt(c2+b2);
        }

        eivecs.col(0) << eivecs.col(1).unitOrthogonal();
      }
    }

    eivals *= scale;
    eivals.array() += shift;

    solver.m_info = Success;
    solver.m_isInitialized = true;
    solver.m_eigenvectorsOk = computeEigenvectors;
  }
};

}

template<typename MatrixType>
EIGEN_DEVICE_FUNC
SelfAdjointEigenSolver<MatrixType>& SelfAdjointEigenSolver<MatrixType>
::computeDirect(const MatrixType& matrix, int options)
{
  internal::direct_selfadjoint_eigenvalues<SelfAdjointEigenSolver,Size,NumTraits<Scalar>::IsComplex>::run(*this,matrix,options);
  return *this;
}

namespace internal {

template<int StorageOrder,typename RealScalar, typename Scalar, typename Index>
EIGEN_DEVICE_FUNC
static void tridiagonal_qr_step(RealScalar* diag, RealScalar* subdiag, Index start, Index end, Scalar* matrixQ, Index n)
{
  
  RealScalar td = (diag[end-1] - diag[end])*RealScalar(0.5);
  RealScalar e = subdiag[end-1];

  RealScalar mu = diag[end];
  if(td==RealScalar(0)) {
    mu -= numext::abs(e);
  } else if (e != RealScalar(0)) {
    const RealScalar e2 = numext::abs2(e);
    const RealScalar h = numext::hypot(td,e);
    if(e2 == RealScalar(0)) {
      mu -= e / ((td + (td>RealScalar(0) ? h : -h)) / e);
    } else {
      mu -= e2 / (td + (td>RealScalar(0) ? h : -h)); 
    }
  }

  RealScalar x = diag[start] - mu;
  RealScalar z = subdiag[start];

  for (Index k = start; k < end && z != RealScalar(0); ++k)
  {
    JacobiRotation<RealScalar> rot;
    rot.makeGivens(x, z);

    RealScalar sdk = rot.s() * diag[k] + rot.c() * subdiag[k];
    RealScalar dkp1 = rot.s() * subdiag[k] + rot.c() * diag[k+1];

    diag[k] = rot.c() * (rot.c() * diag[k] - rot.s() * subdiag[k]) - rot.s() * (rot.c() * subdiag[k] - rot.s() * diag[k+1]);
    diag[k+1] = rot.s() * sdk + rot.c() * dkp1;
    subdiag[k] = rot.c() * sdk - rot.s() * dkp1;
    
    if (k > start)
      subdiag[k - 1] = rot.c() * subdiag[k-1] - rot.s() * z;

    x = subdiag[k];
    if (k < end - 1)
    {
      z = -rot.s() * subdiag[k+1];
      subdiag[k + 1] = rot.c() * subdiag[k+1];
    }

    if (matrixQ)
    {
      
      Map<Matrix<Scalar,Dynamic,Dynamic,StorageOrder> > q(matrixQ,n,n);
      q.applyOnTheRight(k,k+1,rot);
    }
  }
}

} 

} 

#endif 
