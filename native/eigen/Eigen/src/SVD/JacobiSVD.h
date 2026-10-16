#ifndef EIGEN_JACOBISVD_H
#define EIGEN_JACOBISVD_H

namespace Eigen { 

namespace internal {

template<typename MatrixType, int QRPreconditioner,
         bool IsComplex = NumTraits<typename MatrixType::Scalar>::IsComplex>
struct svd_precondition_2x2_block_to_be_real {};

enum { PreconditionIfMoreColsThanRows, PreconditionIfMoreRowsThanCols };

template<typename MatrixType, int QRPreconditioner, int Case>
struct qr_preconditioner_should_do_anything
{
  enum { a = MatrixType::RowsAtCompileTime != Dynamic &&
             MatrixType::ColsAtCompileTime != Dynamic &&
             MatrixType::ColsAtCompileTime <= MatrixType::RowsAtCompileTime,
         b = MatrixType::RowsAtCompileTime != Dynamic &&
             MatrixType::ColsAtCompileTime != Dynamic &&
             MatrixType::RowsAtCompileTime <= MatrixType::ColsAtCompileTime,
         ret = !( (QRPreconditioner == NoQRPreconditioner) ||
                  (Case == PreconditionIfMoreColsThanRows && bool(a)) ||
                  (Case == PreconditionIfMoreRowsThanCols && bool(b)) )
  };
};

template<typename MatrixType, int QRPreconditioner, int Case,
         bool DoAnything = qr_preconditioner_should_do_anything<MatrixType, QRPreconditioner, Case>::ret
> struct qr_preconditioner_impl {};

template<typename MatrixType, int QRPreconditioner, int Case>
class qr_preconditioner_impl<MatrixType, QRPreconditioner, Case, false>
{
public:
  void allocate(const JacobiSVD<MatrixType, QRPreconditioner>&) {}
  bool run(JacobiSVD<MatrixType, QRPreconditioner>&, const MatrixType&)
  {
    return false;
  }
};

template<typename MatrixType>
class qr_preconditioner_impl<MatrixType, FullPivHouseholderQRPreconditioner, PreconditionIfMoreRowsThanCols, true>
{
public:
  typedef typename MatrixType::Scalar Scalar;
  enum
  {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime
  };
  typedef Matrix<Scalar, 1, RowsAtCompileTime, RowMajor, 1, MaxRowsAtCompileTime> WorkspaceType;

  void allocate(const JacobiSVD<MatrixType, FullPivHouseholderQRPreconditioner>& svd)
  {
    if (svd.rows() != m_qr.rows() || svd.cols() != m_qr.cols())
    {
      m_qr.~QRType();
      ::new (&m_qr) QRType(svd.rows(), svd.cols());
    }
    if (svd.m_computeFullU) m_workspace.resize(svd.rows());
  }

  bool run(JacobiSVD<MatrixType, FullPivHouseholderQRPreconditioner>& svd, const MatrixType& matrix)
  {
    if(matrix.rows() > matrix.cols())
    {
      m_qr.compute(matrix);
      svd.m_workMatrix = m_qr.matrixQR().block(0,0,matrix.cols(),matrix.cols()).template triangularView<Upper>();
      if(svd.m_computeFullU) m_qr.matrixQ().evalTo(svd.m_matrixU, m_workspace);
      if(svd.computeV()) svd.m_matrixV = m_qr.colsPermutation();
      return true;
    }
    return false;
  }
private:
  typedef FullPivHouseholderQR<MatrixType> QRType;
  QRType m_qr;
  WorkspaceType m_workspace;
};

template<typename MatrixType>
class qr_preconditioner_impl<MatrixType, FullPivHouseholderQRPreconditioner, PreconditionIfMoreColsThanRows, true>
{
public:
  typedef typename MatrixType::Scalar Scalar;
  enum
  {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime,
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime,
    Options = MatrixType::Options
  };

  typedef typename internal::make_proper_matrix_type<
    Scalar, ColsAtCompileTime, RowsAtCompileTime, Options, MaxColsAtCompileTime, MaxRowsAtCompileTime
  >::type TransposeTypeWithSameStorageOrder;

  void allocate(const JacobiSVD<MatrixType, FullPivHouseholderQRPreconditioner>& svd)
  {
    if (svd.cols() != m_qr.rows() || svd.rows() != m_qr.cols())
    {
      m_qr.~QRType();
      ::new (&m_qr) QRType(svd.cols(), svd.rows());
    }
    m_adjoint.resize(svd.cols(), svd.rows());
    if (svd.m_computeFullV) m_workspace.resize(svd.cols());
  }

  bool run(JacobiSVD<MatrixType, FullPivHouseholderQRPreconditioner>& svd, const MatrixType& matrix)
  {
    if(matrix.cols() > matrix.rows())
    {
      m_adjoint = matrix.adjoint();
      m_qr.compute(m_adjoint);
      svd.m_workMatrix = m_qr.matrixQR().block(0,0,matrix.rows(),matrix.rows()).template triangularView<Upper>().adjoint();
      if(svd.m_computeFullV) m_qr.matrixQ().evalTo(svd.m_matrixV, m_workspace);
      if(svd.computeU()) svd.m_matrixU = m_qr.colsPermutation();
      return true;
    }
    else return false;
  }
private:
  typedef FullPivHouseholderQR<TransposeTypeWithSameStorageOrder> QRType;
  QRType m_qr;
  TransposeTypeWithSameStorageOrder m_adjoint;
  typename internal::plain_row_type<MatrixType>::type m_workspace;
};

template<typename MatrixType>
class qr_preconditioner_impl<MatrixType, ColPivHouseholderQRPreconditioner, PreconditionIfMoreRowsThanCols, true>
{
public:
  void allocate(const JacobiSVD<MatrixType, ColPivHouseholderQRPreconditioner>& svd)
  {
    if (svd.rows() != m_qr.rows() || svd.cols() != m_qr.cols())
    {
      m_qr.~QRType();
      ::new (&m_qr) QRType(svd.rows(), svd.cols());
    }
    if (svd.m_computeFullU) m_workspace.resize(svd.rows());
    else if (svd.m_computeThinU) m_workspace.resize(svd.cols());
  }

  bool run(JacobiSVD<MatrixType, ColPivHouseholderQRPreconditioner>& svd, const MatrixType& matrix)
  {
    if(matrix.rows() > matrix.cols())
    {
      m_qr.compute(matrix);
      svd.m_workMatrix = m_qr.matrixQR().block(0,0,matrix.cols(),matrix.cols()).template triangularView<Upper>();
      if(svd.m_computeFullU) m_qr.householderQ().evalTo(svd.m_matrixU, m_workspace);
      else if(svd.m_computeThinU)
      {
        svd.m_matrixU.setIdentity(matrix.rows(), matrix.cols());
        m_qr.householderQ().applyThisOnTheLeft(svd.m_matrixU, m_workspace);
      }
      if(svd.computeV()) svd.m_matrixV = m_qr.colsPermutation();
      return true;
    }
    return false;
  }

private:
  typedef ColPivHouseholderQR<MatrixType> QRType;
  QRType m_qr;
  typename internal::plain_col_type<MatrixType>::type m_workspace;
};

template<typename MatrixType>
class qr_preconditioner_impl<MatrixType, ColPivHouseholderQRPreconditioner, PreconditionIfMoreColsThanRows, true>
{
public:
  typedef typename MatrixType::Scalar Scalar;
  enum
  {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime,
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime,
    Options = MatrixType::Options
  };

  typedef typename internal::make_proper_matrix_type<
    Scalar, ColsAtCompileTime, RowsAtCompileTime, Options, MaxColsAtCompileTime, MaxRowsAtCompileTime
  >::type TransposeTypeWithSameStorageOrder;

  void allocate(const JacobiSVD<MatrixType, ColPivHouseholderQRPreconditioner>& svd)
  {
    if (svd.cols() != m_qr.rows() || svd.rows() != m_qr.cols())
    {
      m_qr.~QRType();
      ::new (&m_qr) QRType(svd.cols(), svd.rows());
    }
    if (svd.m_computeFullV) m_workspace.resize(svd.cols());
    else if (svd.m_computeThinV) m_workspace.resize(svd.rows());
    m_adjoint.resize(svd.cols(), svd.rows());
  }

  bool run(JacobiSVD<MatrixType, ColPivHouseholderQRPreconditioner>& svd, const MatrixType& matrix)
  {
    if(matrix.cols() > matrix.rows())
    {
      m_adjoint = matrix.adjoint();
      m_qr.compute(m_adjoint);

      svd.m_workMatrix = m_qr.matrixQR().block(0,0,matrix.rows(),matrix.rows()).template triangularView<Upper>().adjoint();
      if(svd.m_computeFullV) m_qr.householderQ().evalTo(svd.m_matrixV, m_workspace);
      else if(svd.m_computeThinV)
      {
        svd.m_matrixV.setIdentity(matrix.cols(), matrix.rows());
        m_qr.householderQ().applyThisOnTheLeft(svd.m_matrixV, m_workspace);
      }
      if(svd.computeU()) svd.m_matrixU = m_qr.colsPermutation();
      return true;
    }
    else return false;
  }

private:
  typedef ColPivHouseholderQR<TransposeTypeWithSameStorageOrder> QRType;
  QRType m_qr;
  TransposeTypeWithSameStorageOrder m_adjoint;
  typename internal::plain_row_type<MatrixType>::type m_workspace;
};

template<typename MatrixType>
class qr_preconditioner_impl<MatrixType, HouseholderQRPreconditioner, PreconditionIfMoreRowsThanCols, true>
{
public:
  void allocate(const JacobiSVD<MatrixType, HouseholderQRPreconditioner>& svd)
  {
    if (svd.rows() != m_qr.rows() || svd.cols() != m_qr.cols())
    {
      m_qr.~QRType();
      ::new (&m_qr) QRType(svd.rows(), svd.cols());
    }
    if (svd.m_computeFullU) m_workspace.resize(svd.rows());
    else if (svd.m_computeThinU) m_workspace.resize(svd.cols());
  }

  bool run(JacobiSVD<MatrixType, HouseholderQRPreconditioner>& svd, const MatrixType& matrix)
  {
    if(matrix.rows() > matrix.cols())
    {
      m_qr.compute(matrix);
      svd.m_workMatrix = m_qr.matrixQR().block(0,0,matrix.cols(),matrix.cols()).template triangularView<Upper>();
      if(svd.m_computeFullU) m_qr.householderQ().evalTo(svd.m_matrixU, m_workspace);
      else if(svd.m_computeThinU)
      {
        svd.m_matrixU.setIdentity(matrix.rows(), matrix.cols());
        m_qr.householderQ().applyThisOnTheLeft(svd.m_matrixU, m_workspace);
      }
      if(svd.computeV()) svd.m_matrixV.setIdentity(matrix.cols(), matrix.cols());
      return true;
    }
    return false;
  }
private:
  typedef HouseholderQR<MatrixType> QRType;
  QRType m_qr;
  typename internal::plain_col_type<MatrixType>::type m_workspace;
};

template<typename MatrixType>
class qr_preconditioner_impl<MatrixType, HouseholderQRPreconditioner, PreconditionIfMoreColsThanRows, true>
{
public:
  typedef typename MatrixType::Scalar Scalar;
  enum
  {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime,
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime,
    Options = MatrixType::Options
  };

  typedef typename internal::make_proper_matrix_type<
    Scalar, ColsAtCompileTime, RowsAtCompileTime, Options, MaxColsAtCompileTime, MaxRowsAtCompileTime
  >::type TransposeTypeWithSameStorageOrder;

  void allocate(const JacobiSVD<MatrixType, HouseholderQRPreconditioner>& svd)
  {
    if (svd.cols() != m_qr.rows() || svd.rows() != m_qr.cols())
    {
      m_qr.~QRType();
      ::new (&m_qr) QRType(svd.cols(), svd.rows());
    }
    if (svd.m_computeFullV) m_workspace.resize(svd.cols());
    else if (svd.m_computeThinV) m_workspace.resize(svd.rows());
    m_adjoint.resize(svd.cols(), svd.rows());
  }

  bool run(JacobiSVD<MatrixType, HouseholderQRPreconditioner>& svd, const MatrixType& matrix)
  {
    if(matrix.cols() > matrix.rows())
    {
      m_adjoint = matrix.adjoint();
      m_qr.compute(m_adjoint);

      svd.m_workMatrix = m_qr.matrixQR().block(0,0,matrix.rows(),matrix.rows()).template triangularView<Upper>().adjoint();
      if(svd.m_computeFullV) m_qr.householderQ().evalTo(svd.m_matrixV, m_workspace);
      else if(svd.m_computeThinV)
      {
        svd.m_matrixV.setIdentity(matrix.cols(), matrix.rows());
        m_qr.householderQ().applyThisOnTheLeft(svd.m_matrixV, m_workspace);
      }
      if(svd.computeU()) svd.m_matrixU.setIdentity(matrix.rows(), matrix.rows());
      return true;
    }
    else return false;
  }

private:
  typedef HouseholderQR<TransposeTypeWithSameStorageOrder> QRType;
  QRType m_qr;
  TransposeTypeWithSameStorageOrder m_adjoint;
  typename internal::plain_row_type<MatrixType>::type m_workspace;
};

template<typename MatrixType, int QRPreconditioner>
struct svd_precondition_2x2_block_to_be_real<MatrixType, QRPreconditioner, false>
{
  typedef JacobiSVD<MatrixType, QRPreconditioner> SVD;
  typedef typename MatrixType::RealScalar RealScalar;
  static bool run(typename SVD::WorkMatrixType&, SVD&, Index, Index, RealScalar&) { return true; }
};

template<typename MatrixType, int QRPreconditioner>
struct svd_precondition_2x2_block_to_be_real<MatrixType, QRPreconditioner, true>
{
  typedef JacobiSVD<MatrixType, QRPreconditioner> SVD;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;
  static bool run(typename SVD::WorkMatrixType& work_matrix, SVD& svd, Index p, Index q, RealScalar& maxDiagEntry)
  {
    using std::sqrt;
    using std::abs;
    Scalar z;
    JacobiRotation<Scalar> rot;
    RealScalar n = sqrt(numext::abs2(work_matrix.coeff(p,p)) + numext::abs2(work_matrix.coeff(q,p)));

    const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();
    const RealScalar precision = NumTraits<Scalar>::epsilon();

    if(n==0)
    {
      
      work_matrix.coeffRef(p,p) = work_matrix.coeffRef(q,p) = Scalar(0);

      if(abs(numext::imag(work_matrix.coeff(p,q)))>considerAsZero)
      {
        
        z = abs(work_matrix.coeff(p,q)) / work_matrix.coeff(p,q);
        work_matrix.row(p) *= z;
        if(svd.computeU()) svd.m_matrixU.col(p) *= conj(z);
      }
      if(abs(numext::imag(work_matrix.coeff(q,q)))>considerAsZero)
      {
        z = abs(work_matrix.coeff(q,q)) / work_matrix.coeff(q,q);
        work_matrix.row(q) *= z;
        if(svd.computeU()) svd.m_matrixU.col(q) *= conj(z);
      }
      
    }
    else
    {
      rot.c() = conj(work_matrix.coeff(p,p)) / n;
      rot.s() = work_matrix.coeff(q,p) / n;
      work_matrix.applyOnTheLeft(p,q,rot);
      if(svd.computeU()) svd.m_matrixU.applyOnTheRight(p,q,rot.adjoint());
      if(abs(numext::imag(work_matrix.coeff(p,q)))>considerAsZero)
      {
        z = abs(work_matrix.coeff(p,q)) / work_matrix.coeff(p,q);
        work_matrix.col(q) *= z;
        if(svd.computeV()) svd.m_matrixV.col(q) *= z;
      }
      if(abs(numext::imag(work_matrix.coeff(q,q)))>considerAsZero)
      {
        z = abs(work_matrix.coeff(q,q)) / work_matrix.coeff(q,q);
        work_matrix.row(q) *= z;
        if(svd.computeU()) svd.m_matrixU.col(q) *= conj(z);
      }
    }

    maxDiagEntry = numext::maxi<RealScalar>(maxDiagEntry,numext::maxi<RealScalar>(abs(work_matrix.coeff(p,p)), abs(work_matrix.coeff(q,q))));
    
    RealScalar threshold = numext::maxi<RealScalar>(considerAsZero, precision * maxDiagEntry);
    return abs(work_matrix.coeff(p,q))>threshold || abs(work_matrix.coeff(q,p)) > threshold;
  }
};

template<typename MatrixType_, int QRPreconditioner>
struct traits<JacobiSVD<MatrixType_,QRPreconditioner> >
        : traits<MatrixType_>
{
  typedef MatrixType_ MatrixType;
};

} 

template<typename MatrixType_, int QRPreconditioner> class JacobiSVD
 : public SVDBase<JacobiSVD<MatrixType_,QRPreconditioner> >
{
    typedef SVDBase<JacobiSVD> Base;
  public:

    typedef MatrixType_ MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename NumTraits<typename MatrixType::Scalar>::Real RealScalar;
    enum {
      RowsAtCompileTime = MatrixType::RowsAtCompileTime,
      ColsAtCompileTime = MatrixType::ColsAtCompileTime,
      DiagSizeAtCompileTime = EIGEN_SIZE_MIN_PREFER_DYNAMIC(RowsAtCompileTime,ColsAtCompileTime),
      MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime,
      MaxDiagSizeAtCompileTime = EIGEN_SIZE_MIN_PREFER_FIXED(MaxRowsAtCompileTime,MaxColsAtCompileTime),
      MatrixOptions = MatrixType::Options
    };

    typedef typename Base::MatrixUType MatrixUType;
    typedef typename Base::MatrixVType MatrixVType;
    typedef typename Base::SingularValuesType SingularValuesType;
    
    typedef typename internal::plain_row_type<MatrixType>::type RowType;
    typedef typename internal::plain_col_type<MatrixType>::type ColType;
    typedef Matrix<Scalar, DiagSizeAtCompileTime, DiagSizeAtCompileTime,
                   MatrixOptions, MaxDiagSizeAtCompileTime, MaxDiagSizeAtCompileTime>
            WorkMatrixType;

    JacobiSVD()
    {}

    JacobiSVD(Index rows, Index cols, unsigned int computationOptions = 0)
    {
      allocate(rows, cols, computationOptions);
    }

    explicit JacobiSVD(const MatrixType& matrix, unsigned int computationOptions = 0)
    {
      compute(matrix, computationOptions);
    }

    JacobiSVD& compute(const MatrixType& matrix, unsigned int computationOptions);

    JacobiSVD& compute(const MatrixType& matrix)
    {
      return compute(matrix, m_computationOptions);
    }

    using Base::computeU;
    using Base::computeV;
    using Base::rows;
    using Base::cols;
    using Base::rank;

  private:
    void allocate(Index rows, Index cols, unsigned int computationOptions);

  protected:
    using Base::m_matrixU;
    using Base::m_matrixV;
    using Base::m_singularValues;
    using Base::m_info;
    using Base::m_isInitialized;
    using Base::m_isAllocated;
    using Base::m_usePrescribedThreshold;
    using Base::m_computeFullU;
    using Base::m_computeThinU;
    using Base::m_computeFullV;
    using Base::m_computeThinV;
    using Base::m_computationOptions;
    using Base::m_nonzeroSingularValues;
    using Base::m_rows;
    using Base::m_cols;
    using Base::m_diagSize;
    using Base::m_prescribedThreshold;
    WorkMatrixType m_workMatrix;

    template<typename MatrixType__, int QRPreconditioner_, bool IsComplex_>
    friend struct internal::svd_precondition_2x2_block_to_be_real;
    template<typename MatrixType__, int QRPreconditioner_, int Case_, bool DoAnything_>
    friend struct internal::qr_preconditioner_impl;

    internal::qr_preconditioner_impl<MatrixType, QRPreconditioner, internal::PreconditionIfMoreColsThanRows> m_qr_precond_morecols;
    internal::qr_preconditioner_impl<MatrixType, QRPreconditioner, internal::PreconditionIfMoreRowsThanCols> m_qr_precond_morerows;
    MatrixType m_scaledMatrix;
};

template<typename MatrixType, int QRPreconditioner>
void JacobiSVD<MatrixType, QRPreconditioner>::allocate(Eigen::Index rows, Eigen::Index cols, unsigned int computationOptions)
{
  eigen_assert(rows >= 0 && cols >= 0);

  if (m_isAllocated &&
      rows == m_rows &&
      cols == m_cols &&
      computationOptions == m_computationOptions)
  {
    return;
  }

  m_rows = rows;
  m_cols = cols;
  m_info = Success;
  m_isInitialized = false;
  m_isAllocated = true;
  m_computationOptions = computationOptions;
  m_computeFullU = (computationOptions & ComputeFullU) != 0;
  m_computeThinU = (computationOptions & ComputeThinU) != 0;
  m_computeFullV = (computationOptions & ComputeFullV) != 0;
  m_computeThinV = (computationOptions & ComputeThinV) != 0;
  eigen_assert(!(m_computeFullU && m_computeThinU) && "JacobiSVD: you can't ask for both full and thin U");
  eigen_assert(!(m_computeFullV && m_computeThinV) && "JacobiSVD: you can't ask for both full and thin V");
  eigen_assert(EIGEN_IMPLIES(m_computeThinU || m_computeThinV, MatrixType::ColsAtCompileTime==Dynamic) &&
              "JacobiSVD: thin U and V are only available when your matrix has a dynamic number of columns.");
  if (QRPreconditioner == FullPivHouseholderQRPreconditioner)
  {
      eigen_assert(!(m_computeThinU || m_computeThinV) &&
              "JacobiSVD: can't compute thin U or thin V with the FullPivHouseholderQR preconditioner. "
              "Use the ColPivHouseholderQR preconditioner instead.");
  }
  m_diagSize = (std::min)(m_rows, m_cols);
  m_singularValues.resize(m_diagSize);
  if(RowsAtCompileTime==Dynamic)
    m_matrixU.resize(m_rows, m_computeFullU ? m_rows
                            : m_computeThinU ? m_diagSize
                            : 0);
  if(ColsAtCompileTime==Dynamic)
    m_matrixV.resize(m_cols, m_computeFullV ? m_cols
                            : m_computeThinV ? m_diagSize
                            : 0);
  m_workMatrix.resize(m_diagSize, m_diagSize);
  
  if(m_cols>m_rows)   m_qr_precond_morecols.allocate(*this);
  if(m_rows>m_cols)   m_qr_precond_morerows.allocate(*this);
  if(m_rows!=m_cols)  m_scaledMatrix.resize(rows,cols);
}

template<typename MatrixType, int QRPreconditioner>
JacobiSVD<MatrixType, QRPreconditioner>&
JacobiSVD<MatrixType, QRPreconditioner>::compute(const MatrixType& matrix, unsigned int computationOptions)
{
  using std::abs;
  allocate(matrix.rows(), matrix.cols(), computationOptions);

  const RealScalar precision = RealScalar(2) * NumTraits<Scalar>::epsilon();

  const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();

  RealScalar scale = matrix.cwiseAbs().template maxCoeff<PropagateNaN>();
  if (!(numext::isfinite)(scale)) {
    m_isInitialized = true;
    m_info = InvalidInput;
    return *this;
  }
  if(scale==RealScalar(0)) scale = RealScalar(1);

  if(m_rows!=m_cols)
  {
    m_scaledMatrix = matrix / scale;
    m_qr_precond_morecols.run(*this, m_scaledMatrix);
    m_qr_precond_morerows.run(*this, m_scaledMatrix);
  }
  else
  {
    m_workMatrix = matrix.block(0,0,m_diagSize,m_diagSize) / scale;
    if(m_computeFullU) m_matrixU.setIdentity(m_rows,m_rows);
    if(m_computeThinU) m_matrixU.setIdentity(m_rows,m_diagSize);
    if(m_computeFullV) m_matrixV.setIdentity(m_cols,m_cols);
    if(m_computeThinV) m_matrixV.setIdentity(m_cols, m_diagSize);
  }

  RealScalar maxDiagEntry = m_workMatrix.cwiseAbs().diagonal().maxCoeff();

  bool finished = false;
  while(!finished)
  {
    finished = true;

    for(Index p = 1; p < m_diagSize; ++p)
    {
      for(Index q = 0; q < p; ++q)
      {

        RealScalar threshold = numext::maxi<RealScalar>(considerAsZero, precision * maxDiagEntry);
        if(abs(m_workMatrix.coeff(p,q))>threshold || abs(m_workMatrix.coeff(q,p)) > threshold)
        {
          finished = false;

          if(internal::svd_precondition_2x2_block_to_be_real<MatrixType, QRPreconditioner>::run(m_workMatrix, *this, p, q, maxDiagEntry))
          {
            JacobiRotation<RealScalar> j_left, j_right;
            internal::real_2x2_jacobi_svd(m_workMatrix, p, q, &j_left, &j_right);

            m_workMatrix.applyOnTheLeft(p,q,j_left);
            if(computeU()) m_matrixU.applyOnTheRight(p,q,j_left.transpose());

            m_workMatrix.applyOnTheRight(p,q,j_right);
            if(computeV()) m_matrixV.applyOnTheRight(p,q,j_right);

            maxDiagEntry = numext::maxi<RealScalar>(maxDiagEntry,numext::maxi<RealScalar>(abs(m_workMatrix.coeff(p,p)), abs(m_workMatrix.coeff(q,q))));
          }
        }
      }
    }
  }

  for(Index i = 0; i < m_diagSize; ++i)
  {

    if(NumTraits<Scalar>::IsComplex && abs(numext::imag(m_workMatrix.coeff(i,i)))>considerAsZero)
    {
      RealScalar a = abs(m_workMatrix.coeff(i,i));
      m_singularValues.coeffRef(i) = abs(a);
      if(computeU()) m_matrixU.col(i) *= m_workMatrix.coeff(i,i)/a;
    }
    else
    {
      
      RealScalar a = numext::real(m_workMatrix.coeff(i,i));
      m_singularValues.coeffRef(i) = abs(a);
      if(computeU() && (a<RealScalar(0))) m_matrixU.col(i) = -m_matrixU.col(i);
    }
  }
  
  m_singularValues *= scale;

  m_nonzeroSingularValues = m_diagSize;
  for(Index i = 0; i < m_diagSize; i++)
  {
    Index pos;
    RealScalar maxRemainingSingularValue = m_singularValues.tail(m_diagSize-i).maxCoeff(&pos);
    if(maxRemainingSingularValue == RealScalar(0))
    {
      m_nonzeroSingularValues = i;
      break;
    }
    if(pos)
    {
      pos += i;
      std::swap(m_singularValues.coeffRef(i), m_singularValues.coeffRef(pos));
      if(computeU()) m_matrixU.col(pos).swap(m_matrixU.col(i));
      if(computeV()) m_matrixV.col(pos).swap(m_matrixV.col(i));
    }
  }

  m_isInitialized = true;
  return *this;
}

template<typename Derived>
JacobiSVD<typename MatrixBase<Derived>::PlainObject>
MatrixBase<Derived>::jacobiSvd(unsigned int computationOptions) const
{
  return JacobiSVD<PlainObject>(*this, computationOptions);
}

} 

#endif 
