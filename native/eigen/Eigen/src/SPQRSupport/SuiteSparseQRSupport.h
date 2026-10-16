#ifndef EIGEN_SUITESPARSEQRSUPPORT_H
#define EIGEN_SUITESPARSEQRSUPPORT_H

namespace Eigen {
  
  template<typename MatrixType> class SPQR; 
  template<typename SPQRType> struct SPQRMatrixQReturnType; 
  template<typename SPQRType> struct SPQRMatrixQTransposeReturnType; 
  template <typename SPQRType, typename Derived> struct SPQR_QProduct;
  namespace internal {
    template <typename SPQRType> struct traits<SPQRMatrixQReturnType<SPQRType> >
    {
      typedef typename SPQRType::MatrixType ReturnType;
    };
    template <typename SPQRType> struct traits<SPQRMatrixQTransposeReturnType<SPQRType> >
    {
      typedef typename SPQRType::MatrixType ReturnType;
    };
    template <typename SPQRType, typename Derived> struct traits<SPQR_QProduct<SPQRType, Derived> >
    {
      typedef typename Derived::PlainObject ReturnType;
    };
  } 

template<typename MatrixType_>
class SPQR : public SparseSolverBase<SPQR<MatrixType_> >
{
  protected:
    typedef SparseSolverBase<SPQR<MatrixType_> > Base;
    using Base::m_isInitialized;
  public:
    typedef typename MatrixType_::Scalar Scalar;
    typedef typename MatrixType_::RealScalar RealScalar;
    typedef SuiteSparse_long StorageIndex ;
    typedef SparseMatrix<Scalar, ColMajor, StorageIndex> MatrixType;
    typedef Map<PermutationMatrix<Dynamic, Dynamic, StorageIndex> > PermutationType;
    enum {
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };
  public:
    SPQR() 
      : m_analysisIsOk(false),
        m_factorizationIsOk(false),
        m_isRUpToDate(false),
        m_ordering(SPQR_ORDERING_DEFAULT),
        m_allow_tol(SPQR_DEFAULT_TOL),
        m_tolerance (NumTraits<Scalar>::epsilon()),
        m_cR(0),
        m_E(0),
        m_H(0),
        m_HPinv(0),
        m_HTau(0),
        m_useDefaultThreshold(true)
    { 
      cholmod_l_start(&m_cc);
    }
    
    explicit SPQR(const MatrixType_& matrix)
      : m_analysisIsOk(false),
        m_factorizationIsOk(false),
        m_isRUpToDate(false),
        m_ordering(SPQR_ORDERING_DEFAULT),
        m_allow_tol(SPQR_DEFAULT_TOL),
        m_tolerance (NumTraits<Scalar>::epsilon()),
        m_cR(0),
        m_E(0),
        m_H(0),
        m_HPinv(0),
        m_HTau(0),
        m_useDefaultThreshold(true)
    {
      cholmod_l_start(&m_cc);
      compute(matrix);
    }
    
    ~SPQR()
    {
      SPQR_free();
      cholmod_l_finish(&m_cc);
    }
    void SPQR_free()
    {
      cholmod_l_free_sparse(&m_H, &m_cc);
      cholmod_l_free_sparse(&m_cR, &m_cc);
      cholmod_l_free_dense(&m_HTau, &m_cc);
      std::free(m_E);
      std::free(m_HPinv);
    }

    void compute(const MatrixType_& matrix)
    {
      if(m_isInitialized) SPQR_free();

      MatrixType mat(matrix);

      RealScalar pivotThreshold = m_tolerance;
      if(m_useDefaultThreshold) 
      {
        RealScalar max2Norm = 0.0;
        for (int j = 0; j < mat.cols(); j++) max2Norm = numext::maxi(max2Norm, mat.col(j).norm());
        if(max2Norm==RealScalar(0))
          max2Norm = RealScalar(1);
        pivotThreshold = 20 * (mat.rows() + mat.cols()) * max2Norm * NumTraits<RealScalar>::epsilon();
      }
      cholmod_sparse A; 
      A = viewAsCholmod(mat);
      m_rows = matrix.rows();
      Index col = matrix.cols();
      m_rank = SuiteSparseQR<Scalar>(m_ordering, pivotThreshold, col, &A, 
                             &m_cR, &m_E, &m_H, &m_HPinv, &m_HTau, &m_cc);

      if (!m_cR)
      {
        m_info = NumericalIssue;
        m_isInitialized = false;
        return;
      }
      m_info = Success;
      m_isInitialized = true;
      m_isRUpToDate = false;
    }

    inline Index rows() const {return m_rows; }

    inline Index cols() const { return m_cR->ncol; }
    
    template<typename Rhs, typename Dest>
    void _solve_impl(const MatrixBase<Rhs> &b, MatrixBase<Dest> &dest) const
    {
      eigen_assert(m_isInitialized && " The QR factorization should be computed first, call compute()");
      eigen_assert(b.cols()==1 && "This method is for vectors only");

      typename Dest::PlainObject y, y2;
      y = matrixQ().transpose() * b;

      Index rk = this->rank();
      y2 = y;
      y.resize((std::max)(cols(),Index(y.rows())),y.cols());
      y.topRows(rk) = this->matrixR().topLeftCorner(rk, rk).template triangularView<Upper>().solve(y2.topRows(rk));

      for(Index i = 0; i < rk; ++i) dest.row(m_E[i]) = y.row(i);
      for(Index i = rk; i < cols(); ++i) dest.row(m_E[i]).setZero();

      m_info = Success;
    }

    const MatrixType matrixR() const
    {
      eigen_assert(m_isInitialized && " The QR factorization should be computed first, call compute()");
      if(!m_isRUpToDate) {
        m_R = viewAsEigen<Scalar,ColMajor, typename MatrixType::StorageIndex>(*m_cR);
        m_isRUpToDate = true;
      }
      return m_R;
    }
    
    SPQRMatrixQReturnType<SPQR> matrixQ() const
    {
      return SPQRMatrixQReturnType<SPQR>(*this);
    }
    
    PermutationType colsPermutation() const
    { 
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return PermutationType(m_E, m_cR->ncol);
    }

    Index rank() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_cc.SPQR_istat[4];
    }
    
    void setSPQROrdering(int ord) { m_ordering = ord;}
    
    void setPivotThreshold(const RealScalar& tol)
    {
      m_useDefaultThreshold = false;
      m_tolerance = tol;
    }

    cholmod_common *cholmodCommon() const { return &m_cc; }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_info;
    }
  protected:
    bool m_analysisIsOk;
    bool m_factorizationIsOk;
    mutable bool m_isRUpToDate;
    mutable ComputationInfo m_info;
    int m_ordering; 
    int m_allow_tol; 
    RealScalar m_tolerance; 
    mutable cholmod_sparse *m_cR; 
    mutable MatrixType m_R; 
    mutable StorageIndex *m_E; 
    mutable cholmod_sparse *m_H;  
    mutable StorageIndex *m_HPinv; 
    mutable cholmod_dense *m_HTau; 
    mutable Index m_rank; 
    mutable cholmod_common m_cc; 
    bool m_useDefaultThreshold;     
    Index m_rows;
    template<typename ,typename > friend struct SPQR_QProduct;
};

template <typename SPQRType, typename Derived>
struct SPQR_QProduct : ReturnByValue<SPQR_QProduct<SPQRType,Derived> >
{
  typedef typename SPQRType::Scalar Scalar;
  typedef typename SPQRType::StorageIndex StorageIndex;
  
  SPQR_QProduct(const SPQRType& spqr, const Derived& other, bool transpose) : m_spqr(spqr),m_other(other),m_transpose(transpose) {}
  
  inline Index rows() const { return m_transpose ? m_spqr.rows() : m_spqr.cols(); }
  inline Index cols() const { return m_other.cols(); }
  
  template<typename ResType>
  void evalTo(ResType& res) const
  {
    cholmod_dense y_cd;
    cholmod_dense *x_cd; 
    int method = m_transpose ? SPQR_QTX : SPQR_QX; 
    cholmod_common *cc = m_spqr.cholmodCommon();
    y_cd = viewAsCholmod(m_other.const_cast_derived());
    x_cd = SuiteSparseQR_qmult<Scalar>(method, m_spqr.m_H, m_spqr.m_HTau, m_spqr.m_HPinv, &y_cd, cc);
    res = Matrix<Scalar,ResType::RowsAtCompileTime,ResType::ColsAtCompileTime>::Map(reinterpret_cast<Scalar*>(x_cd->x), x_cd->nrow, x_cd->ncol);
    cholmod_l_free_dense(&x_cd, cc);
  }
  const SPQRType& m_spqr; 
  const Derived& m_other; 
  bool m_transpose; 
  
};
template<typename SPQRType>
struct SPQRMatrixQReturnType{
  
  SPQRMatrixQReturnType(const SPQRType& spqr) : m_spqr(spqr) {}
  template<typename Derived>
  SPQR_QProduct<SPQRType, Derived> operator*(const MatrixBase<Derived>& other)
  {
    return SPQR_QProduct<SPQRType,Derived>(m_spqr,other.derived(),false);
  }
  SPQRMatrixQTransposeReturnType<SPQRType> adjoint() const
  {
    return SPQRMatrixQTransposeReturnType<SPQRType>(m_spqr);
  }
  
  SPQRMatrixQTransposeReturnType<SPQRType> transpose() const
  {
    return SPQRMatrixQTransposeReturnType<SPQRType>(m_spqr);
  }
  const SPQRType& m_spqr;
};

template<typename SPQRType>
struct SPQRMatrixQTransposeReturnType{
  SPQRMatrixQTransposeReturnType(const SPQRType& spqr) : m_spqr(spqr) {}
  template<typename Derived>
  SPQR_QProduct<SPQRType,Derived> operator*(const MatrixBase<Derived>& other)
  {
    return SPQR_QProduct<SPQRType,Derived>(m_spqr,other.derived(), true);
  }
  const SPQRType& m_spqr;
};

}
#endif
