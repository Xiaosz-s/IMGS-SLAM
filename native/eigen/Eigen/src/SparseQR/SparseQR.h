#ifndef EIGEN_SPARSE_QR_H
#define EIGEN_SPARSE_QR_H

namespace Eigen {

template<typename MatrixType, typename OrderingType> class SparseQR;
template<typename SparseQRType> struct SparseQRMatrixQReturnType;
template<typename SparseQRType> struct SparseQRMatrixQTransposeReturnType;
template<typename SparseQRType, typename Derived> struct SparseQR_QProduct;
namespace internal {
  template <typename SparseQRType> struct traits<SparseQRMatrixQReturnType<SparseQRType> >
  {
    typedef typename SparseQRType::MatrixType ReturnType;
    typedef typename ReturnType::StorageIndex StorageIndex;
    typedef typename ReturnType::StorageKind StorageKind;
    enum {
      RowsAtCompileTime = Dynamic,
      ColsAtCompileTime = Dynamic
    };
  };
  template <typename SparseQRType> struct traits<SparseQRMatrixQTransposeReturnType<SparseQRType> >
  {
    typedef typename SparseQRType::MatrixType ReturnType;
  };
  template <typename SparseQRType, typename Derived> struct traits<SparseQR_QProduct<SparseQRType, Derived> >
  {
    typedef typename Derived::PlainObject ReturnType;
  };
} 

template<typename MatrixType_, typename OrderingType_>
class SparseQR : public SparseSolverBase<SparseQR<MatrixType_,OrderingType_> >
{
  protected:
    typedef SparseSolverBase<SparseQR<MatrixType_,OrderingType_> > Base;
    using Base::m_isInitialized;
  public:
    using Base::_solve_impl;
    typedef MatrixType_ MatrixType;
    typedef OrderingType_ OrderingType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename MatrixType::StorageIndex StorageIndex;
    typedef SparseMatrix<Scalar,ColMajor,StorageIndex> QRMatrixType;
    typedef Matrix<StorageIndex, Dynamic, 1> IndexVector;
    typedef Matrix<Scalar, Dynamic, 1> ScalarVector;
    typedef PermutationMatrix<Dynamic, Dynamic, StorageIndex> PermutationType;

    enum {
      ColsAtCompileTime = MatrixType::ColsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };
    
  public:
    SparseQR () :  m_analysisIsok(false), m_lastError(""), m_useDefaultThreshold(true),m_isQSorted(false),m_isEtreeOk(false)
    { }

    explicit SparseQR(const MatrixType& mat) : m_analysisIsok(false), m_lastError(""), m_useDefaultThreshold(true),m_isQSorted(false),m_isEtreeOk(false)
    {
      compute(mat);
    }

    void compute(const MatrixType& mat)
    {
      analyzePattern(mat);
      factorize(mat);
    }
    void analyzePattern(const MatrixType& mat);
    void factorize(const MatrixType& mat);

    inline Index rows() const { return m_pmat.rows(); }

    inline Index cols() const { return m_pmat.cols();}

    const QRMatrixType& matrixR() const { return m_R; }

    Index rank() const
    {
      eigen_assert(m_isInitialized && "The factorization should be called first, use compute()");
      return m_nonzeropivots; 
    }

    SparseQRMatrixQReturnType<SparseQR> matrixQ() const 
    { return SparseQRMatrixQReturnType<SparseQR>(*this); }

    const PermutationType& colsPermutation() const
    { 
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_outputPerm_c;
    }

    std::string lastErrorMessage() const { return m_lastError; }

    template<typename Rhs, typename Dest>
    bool _solve_impl(const MatrixBase<Rhs> &B, MatrixBase<Dest> &dest) const
    {
      eigen_assert(m_isInitialized && "The factorization should be called first, use compute()");
      eigen_assert(this->rows() == B.rows() && "SparseQR::solve() : invalid number of rows in the right hand side matrix");

      Index rank = this->rank();

      typename Dest::PlainObject y, b;
      y = this->matrixQ().adjoint() * B;
      b = y;

      y.resize((std::max<Index>)(cols(),y.rows()),y.cols());
      y.topRows(rank) = this->matrixR().topLeftCorner(rank, rank).template triangularView<Upper>().solve(b.topRows(rank));
      y.bottomRows(y.rows()-rank).setZero();

      if (m_perm_c.size())  dest = colsPermutation() * y.topRows(cols());
      else                  dest = y.topRows(cols());
      
      m_info = Success;
      return true;
    }

    void setPivotThreshold(const RealScalar& threshold)
    {
      m_useDefaultThreshold = false;
      m_threshold = threshold;
    }

    template<typename Rhs>
    inline const Solve<SparseQR, Rhs> solve(const MatrixBase<Rhs>& B) const 
    {
      eigen_assert(m_isInitialized && "The factorization should be called first, use compute()");
      eigen_assert(this->rows() == B.rows() && "SparseQR::solve() : invalid number of rows in the right hand side matrix");
      return Solve<SparseQR, Rhs>(*this, B.derived());
    }
    template<typename Rhs>
    inline const Solve<SparseQR, Rhs> solve(const SparseMatrixBase<Rhs>& B) const
    {
          eigen_assert(m_isInitialized && "The factorization should be called first, use compute()");
          eigen_assert(this->rows() == B.rows() && "SparseQR::solve() : invalid number of rows in the right hand side matrix");
          return Solve<SparseQR, Rhs>(*this, B.derived());
    }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_info;
    }

    inline void _sort_matrix_Q()
    {
      if(this->m_isQSorted) return;
      
      SparseMatrix<Scalar, RowMajor, Index> mQrm(this->m_Q);
      this->m_Q = mQrm;
      this->m_isQSorted = true;
    }

  protected:
    bool m_analysisIsok;
    bool m_factorizationIsok;
    mutable ComputationInfo m_info;
    std::string m_lastError;
    QRMatrixType m_pmat;            
    QRMatrixType m_R;               
    QRMatrixType m_Q;               
    ScalarVector m_hcoeffs;         
    PermutationType m_perm_c;       
    PermutationType m_pivotperm;    
    PermutationType m_outputPerm_c; 
    RealScalar m_threshold;         
    bool m_useDefaultThreshold;     
    Index m_nonzeropivots;          
    IndexVector m_etree;            
    IndexVector m_firstRowElt;      
    bool m_isQSorted;               
    bool m_isEtreeOk;               
    
    template <typename, typename > friend struct SparseQR_QProduct;
    
};

template <typename MatrixType, typename OrderingType>
void SparseQR<MatrixType,OrderingType>::analyzePattern(const MatrixType& mat)
{
  eigen_assert(mat.isCompressed() && "SparseQR requires a sparse matrix in compressed mode. Call .makeCompressed() before passing it to SparseQR");
  
  typename internal::conditional<MatrixType::IsRowMajor,QRMatrixType,const MatrixType&>::type matCpy(mat);
  
  OrderingType ord; 
  ord(matCpy, m_perm_c); 
  Index n = mat.cols();
  Index m = mat.rows();
  Index diagSize = (std::min)(m,n);
  
  if (!m_perm_c.size())
  {
    m_perm_c.resize(n);
    m_perm_c.indices().setLinSpaced(n, 0,StorageIndex(n-1));
  }

  m_outputPerm_c = m_perm_c.inverse();
  internal::coletree(matCpy, m_etree, m_firstRowElt, m_outputPerm_c.indices().data());
  m_isEtreeOk = true;
  
  m_R.resize(m, n);
  m_Q.resize(m, diagSize);

  m_R.reserve(2*mat.nonZeros()); 
  m_Q.reserve(2*mat.nonZeros());
  m_hcoeffs.resize(diagSize);
  m_analysisIsok = true;
}

template <typename MatrixType, typename OrderingType>
void SparseQR<MatrixType,OrderingType>::factorize(const MatrixType& mat)
{
  using std::abs;
  
  eigen_assert(m_analysisIsok && "analyzePattern() should be called before this step");
  StorageIndex m = StorageIndex(mat.rows());
  StorageIndex n = StorageIndex(mat.cols());
  StorageIndex diagSize = (std::min)(m,n);
  IndexVector mark((std::max)(m,n)); mark.setConstant(-1);  
  IndexVector Ridx(n), Qidx(m);                             
  Index nzcolR, nzcolQ;                                     
  ScalarVector tval(m);                                     
  RealScalar pivotThreshold = m_threshold;
  
  m_R.setZero();
  m_Q.setZero();
  m_pmat = mat;
  if(!m_isEtreeOk)
  {
    m_outputPerm_c = m_perm_c.inverse();
    internal::coletree(m_pmat, m_etree, m_firstRowElt, m_outputPerm_c.indices().data());
    m_isEtreeOk = true;
  }

  m_pmat.uncompress(); 

  {

    IndexVector originalOuterIndicesCpy;
    const StorageIndex *originalOuterIndices = mat.outerIndexPtr();
    if(MatrixType::IsRowMajor)
    {
      originalOuterIndicesCpy = IndexVector::Map(m_pmat.outerIndexPtr(),n+1);
      originalOuterIndices = originalOuterIndicesCpy.data();
    }
    
    for (int i = 0; i < n; i++)
    {
      Index p = m_perm_c.size() ? m_perm_c.indices()(i) : i;
      m_pmat.outerIndexPtr()[p] = originalOuterIndices[i]; 
      m_pmat.innerNonZeroPtr()[p] = originalOuterIndices[i+1] - originalOuterIndices[i]; 
    }
  }

  if(m_useDefaultThreshold) 
  {
    RealScalar max2Norm = 0.0;
    for (int j = 0; j < n; j++) max2Norm = numext::maxi(max2Norm, m_pmat.col(j).norm());
    if(max2Norm==RealScalar(0))
      max2Norm = RealScalar(1);
    pivotThreshold = 20 * (m + n) * max2Norm * NumTraits<RealScalar>::epsilon();
  }

  m_pivotperm.setIdentity(n);
  
  StorageIndex nonzeroCol = 0; 
  m_Q.startVec(0);

  for (StorageIndex col = 0; col < n; ++col)
  {
    mark.setConstant(-1);
    m_R.startVec(col);
    mark(nonzeroCol) = col;
    Qidx(0) = nonzeroCol;
    nzcolR = 0; nzcolQ = 1;
    bool found_diag = nonzeroCol>=m;
    tval.setZero(); 

    for (typename QRMatrixType::InnerIterator itp(m_pmat, col); itp || !found_diag; ++itp)
    {
      StorageIndex curIdx = nonzeroCol;
      if(itp) curIdx = StorageIndex(itp.row());
      if(curIdx == nonzeroCol) found_diag = true;

      StorageIndex st = m_firstRowElt(curIdx); 
      if (st < 0 )
      {
        m_lastError = "Empty row found during numerical factorization";
        m_info = InvalidInput;
        return;
      }

      Index bi = nzcolR;
      for (; mark(st) != col; st = m_etree(st))
      {
        Ridx(nzcolR) = st;  
        mark(st) = col;     
        nzcolR++;
      }

      Index nt = nzcolR-bi;
      for(Index i = 0; i < nt/2; i++) std::swap(Ridx(bi+i), Ridx(nzcolR-i-1));

      if(itp) tval(curIdx) = itp.value();
      else    tval(curIdx) = Scalar(0);

      if(curIdx > nonzeroCol && mark(curIdx) != col ) 
      {
        Qidx(nzcolQ) = curIdx;  
        mark(curIdx) = col;     
        nzcolQ++;
      }
    }

    for (Index i = nzcolR-1; i >= 0; i--)
    {
      Index curIdx = Ridx(i);

      Scalar tdot(0);

      tdot = m_Q.col(curIdx).dot(tval);

      tdot *= m_hcoeffs(curIdx);

      for (typename QRMatrixType::InnerIterator itq(m_Q, curIdx); itq; ++itq)
        tval(itq.row()) -= itq.value() * tdot;

      if(m_etree(Ridx(i)) == nonzeroCol)
      {
        for (typename QRMatrixType::InnerIterator itq(m_Q, curIdx); itq; ++itq)
        {
          StorageIndex iQ = StorageIndex(itq.row());
          if (mark(iQ) != col)
          {
            Qidx(nzcolQ++) = iQ;  
            mark(iQ) = col;       
          }
        }
      }
    } 
    
    Scalar tau = RealScalar(0);
    RealScalar beta = 0;
    
    if(nonzeroCol < diagSize)
    {

      Scalar c0 = nzcolQ ? tval(Qidx(0)) : Scalar(0);

      RealScalar sqrNorm = 0.;
      for (Index itq = 1; itq < nzcolQ; ++itq) sqrNorm += numext::abs2(tval(Qidx(itq)));
      if(sqrNorm == RealScalar(0) && numext::imag(c0) == RealScalar(0))
      {
        beta = numext::real(c0);
        tval(Qidx(0)) = 1;
      }
      else
      {
        using std::sqrt;
        beta = sqrt(numext::abs2(c0) + sqrNorm);
        if(numext::real(c0) >= RealScalar(0))
          beta = -beta;
        tval(Qidx(0)) = 1;
        for (Index itq = 1; itq < nzcolQ; ++itq)
          tval(Qidx(itq)) /= (c0 - beta);
        tau = numext::conj((beta-c0) / beta);
          
      }
    }

    for (Index  i = nzcolR-1; i >= 0; i--)
    {
      Index curIdx = Ridx(i);
      if(curIdx < nonzeroCol) 
      {
        m_R.insertBackByOuterInnerUnordered(col, curIdx) = tval(curIdx);
        tval(curIdx) = Scalar(0.);
      }
    }

    if(nonzeroCol < diagSize && abs(beta) >= pivotThreshold)
    {
      m_R.insertBackByOuterInner(col, nonzeroCol) = beta;
      
      m_hcoeffs(nonzeroCol) = tau;
      
      for (Index itq = 0; itq < nzcolQ; ++itq)
      {
        Index iQ = Qidx(itq);
        m_Q.insertBackByOuterInnerUnordered(nonzeroCol,iQ) = tval(iQ);
        tval(iQ) = Scalar(0.);
      }
      nonzeroCol++;
      if(nonzeroCol<diagSize)
        m_Q.startVec(nonzeroCol);
    }
    else
    {
      
      for (Index j = nonzeroCol; j < n-1; j++) 
        std::swap(m_pivotperm.indices()(j), m_pivotperm.indices()[j+1]);

      internal::coletree(m_pmat, m_etree, m_firstRowElt, m_pivotperm.indices().data());
      m_isEtreeOk = false;
    }
  }
  
  m_hcoeffs.tail(diagSize-nonzeroCol).setZero();

  m_Q.finalize();
  m_Q.makeCompressed();
  m_R.finalize();
  m_R.makeCompressed();
  m_isQSorted = false;

  m_nonzeropivots = nonzeroCol;
  
  if(nonzeroCol<n)
  {
    
    QRMatrixType tempR(m_R);
    m_R = tempR * m_pivotperm;

    m_outputPerm_c = m_outputPerm_c * m_pivotperm;
  }
  
  m_isInitialized = true; 
  m_factorizationIsok = true;
  m_info = Success;
}

template <typename SparseQRType, typename Derived>
struct SparseQR_QProduct : ReturnByValue<SparseQR_QProduct<SparseQRType, Derived> >
{
  typedef typename SparseQRType::QRMatrixType MatrixType;
  typedef typename SparseQRType::Scalar Scalar;
  
  SparseQR_QProduct(const SparseQRType& qr, const Derived& other, bool transpose) : 
  m_qr(qr),m_other(other),m_transpose(transpose) {}
  inline Index rows() const { return m_qr.matrixQ().rows(); }
  inline Index cols() const { return m_other.cols(); }

  template<typename DesType>
  void evalTo(DesType& res) const
  {
    Index m = m_qr.rows();
    Index n = m_qr.cols();
    Index diagSize = (std::min)(m,n);
    res = m_other;
    if (m_transpose)
    {
      eigen_assert(m_qr.m_Q.rows() == m_other.rows() && "Non conforming object sizes");
      
      for(Index j = 0; j < res.cols(); j++){
        for (Index k = 0; k < diagSize; k++)
        {
          Scalar tau = Scalar(0);
          tau = m_qr.m_Q.col(k).dot(res.col(j));
          if(tau==Scalar(0)) continue;
          tau = tau * m_qr.m_hcoeffs(k);
          res.col(j) -= tau * m_qr.m_Q.col(k);
        }
      }
    }
    else
    {
      eigen_assert(m_qr.matrixQ().cols() == m_other.rows() && "Non conforming object sizes");

      res.conservativeResize(rows(), cols());

      for(Index j = 0; j < res.cols(); j++)
      {
        Index start_k = internal::is_identity<Derived>::value ? numext::mini(j,diagSize-1) : diagSize-1;
        for (Index k = start_k; k >=0; k--)
        {
          Scalar tau = Scalar(0);
          tau = m_qr.m_Q.col(k).dot(res.col(j));
          if(tau==Scalar(0)) continue;
          tau = tau * numext::conj(m_qr.m_hcoeffs(k));
          res.col(j) -= tau * m_qr.m_Q.col(k);
        }
      }
    }
  }
  
  const SparseQRType& m_qr;
  const Derived& m_other;
  bool m_transpose; 
};

template<typename SparseQRType>
struct SparseQRMatrixQReturnType : public EigenBase<SparseQRMatrixQReturnType<SparseQRType> >
{  
  typedef typename SparseQRType::Scalar Scalar;
  typedef Matrix<Scalar,Dynamic,Dynamic> DenseMatrix;
  enum {
    RowsAtCompileTime = Dynamic,
    ColsAtCompileTime = Dynamic
  };
  explicit SparseQRMatrixQReturnType(const SparseQRType& qr) : m_qr(qr) {}
  template<typename Derived>
  SparseQR_QProduct<SparseQRType, Derived> operator*(const MatrixBase<Derived>& other)
  {
    return SparseQR_QProduct<SparseQRType,Derived>(m_qr,other.derived(),false);
  }
  
  SparseQRMatrixQTransposeReturnType<SparseQRType> adjoint() const
  {
    return SparseQRMatrixQTransposeReturnType<SparseQRType>(m_qr);
  }
  inline Index rows() const { return m_qr.rows(); }
  inline Index cols() const { return m_qr.rows(); }
  
  SparseQRMatrixQTransposeReturnType<SparseQRType> transpose() const
  {
    return SparseQRMatrixQTransposeReturnType<SparseQRType>(m_qr);
  }
  const SparseQRType& m_qr;
};

template<typename SparseQRType>
struct SparseQRMatrixQTransposeReturnType
{
  explicit SparseQRMatrixQTransposeReturnType(const SparseQRType& qr) : m_qr(qr) {}
  template<typename Derived>
  SparseQR_QProduct<SparseQRType,Derived> operator*(const MatrixBase<Derived>& other)
  {
    return SparseQR_QProduct<SparseQRType,Derived>(m_qr,other.derived(), true);
  }
  const SparseQRType& m_qr;
};

namespace internal {
  
template<typename SparseQRType>
struct evaluator_traits<SparseQRMatrixQReturnType<SparseQRType> >
{
  typedef typename SparseQRType::MatrixType MatrixType;
  typedef typename storage_kind_to_evaluator_kind<typename MatrixType::StorageKind>::Kind Kind;
  typedef SparseShape Shape;
};

template< typename DstXprType, typename SparseQRType>
struct Assignment<DstXprType, SparseQRMatrixQReturnType<SparseQRType>, internal::assign_op<typename DstXprType::Scalar,typename DstXprType::Scalar>, Sparse2Sparse>
{
  typedef SparseQRMatrixQReturnType<SparseQRType> SrcXprType;
  typedef typename DstXprType::Scalar Scalar;
  typedef typename DstXprType::StorageIndex StorageIndex;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<Scalar,Scalar> &)
  {
    typename DstXprType::PlainObject idMat(src.rows(), src.cols());
    idMat.setIdentity();
    
    const_cast<SparseQRType *>(&src.m_qr)->_sort_matrix_Q();
    dst = SparseQR_QProduct<SparseQRType, DstXprType>(src.m_qr, idMat, false);
  }
};

template< typename DstXprType, typename SparseQRType>
struct Assignment<DstXprType, SparseQRMatrixQReturnType<SparseQRType>, internal::assign_op<typename DstXprType::Scalar,typename DstXprType::Scalar>, Sparse2Dense>
{
  typedef SparseQRMatrixQReturnType<SparseQRType> SrcXprType;
  typedef typename DstXprType::Scalar Scalar;
  typedef typename DstXprType::StorageIndex StorageIndex;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<Scalar,Scalar> &)
  {
    dst = src.m_qr.matrixQ() * DstXprType::Identity(src.m_qr.rows(), src.m_qr.rows());
  }
};

} 

} 

#endif
