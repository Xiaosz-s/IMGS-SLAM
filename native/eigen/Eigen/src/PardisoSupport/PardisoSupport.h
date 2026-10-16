#ifndef EIGEN_PARDISOSUPPORT_H
#define EIGEN_PARDISOSUPPORT_H

namespace Eigen { 

template<typename MatrixType_> class PardisoLU;
template<typename MatrixType_, int Options=Upper> class PardisoLLT;
template<typename MatrixType_, int Options=Upper> class PardisoLDLT;

namespace internal
{
  template<typename IndexType>
  struct pardiso_run_selector
  {
    static IndexType run( _MKL_DSS_HANDLE_t pt, IndexType maxfct, IndexType mnum, IndexType type, IndexType phase, IndexType n, void *a,
                      IndexType *ia, IndexType *ja, IndexType *perm, IndexType nrhs, IndexType *iparm, IndexType msglvl, void *b, void *x)
    {
      IndexType error = 0;
      ::pardiso(pt, &maxfct, &mnum, &type, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x, &error);
      return error;
    }
  };
  template<>
  struct pardiso_run_selector<long long int>
  {
    typedef long long int IndexType;
    static IndexType run( _MKL_DSS_HANDLE_t pt, IndexType maxfct, IndexType mnum, IndexType type, IndexType phase, IndexType n, void *a,
                      IndexType *ia, IndexType *ja, IndexType *perm, IndexType nrhs, IndexType *iparm, IndexType msglvl, void *b, void *x)
    {
      IndexType error = 0;
      ::pardiso_64(pt, &maxfct, &mnum, &type, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x, &error);
      return error;
    }
  };

  template<class Pardiso> struct pardiso_traits;

  template<typename MatrixType_>
  struct pardiso_traits< PardisoLU<MatrixType_> >
  {
    typedef MatrixType_ MatrixType;
    typedef typename MatrixType_::Scalar Scalar;
    typedef typename MatrixType_::RealScalar RealScalar;
    typedef typename MatrixType_::StorageIndex StorageIndex;
  };

  template<typename MatrixType_, int Options>
  struct pardiso_traits< PardisoLLT<MatrixType_, Options> >
  {
    typedef MatrixType_ MatrixType;
    typedef typename MatrixType_::Scalar Scalar;
    typedef typename MatrixType_::RealScalar RealScalar;
    typedef typename MatrixType_::StorageIndex StorageIndex;
  };

  template<typename MatrixType_, int Options>
  struct pardiso_traits< PardisoLDLT<MatrixType_, Options> >
  {
    typedef MatrixType_ MatrixType;
    typedef typename MatrixType_::Scalar Scalar;
    typedef typename MatrixType_::RealScalar RealScalar;
    typedef typename MatrixType_::StorageIndex StorageIndex;
  };

} 

template<class Derived>
class PardisoImpl : public SparseSolverBase<Derived>
{
  protected:
    typedef SparseSolverBase<Derived> Base;
    using Base::derived;
    using Base::m_isInitialized;
    
    typedef internal::pardiso_traits<Derived> Traits;
  public:
    using Base::_solve_impl;
    
    typedef typename Traits::MatrixType MatrixType;
    typedef typename Traits::Scalar Scalar;
    typedef typename Traits::RealScalar RealScalar;
    typedef typename Traits::StorageIndex StorageIndex;
    typedef SparseMatrix<Scalar,RowMajor,StorageIndex> SparseMatrixType;
    typedef Matrix<Scalar,Dynamic,1> VectorType;
    typedef Matrix<StorageIndex, 1, MatrixType::ColsAtCompileTime> IntRowVectorType;
    typedef Matrix<StorageIndex, MatrixType::RowsAtCompileTime, 1> IntColVectorType;
    typedef Array<StorageIndex,64,1,DontAlign> ParameterType;
    enum {
      ScalarIsComplex = NumTraits<Scalar>::IsComplex,
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };

    PardisoImpl()
      : m_analysisIsOk(false), m_factorizationIsOk(false)
    {
      eigen_assert((sizeof(StorageIndex) >= sizeof(_INTEGER_t) && sizeof(StorageIndex) <= 8) && "Non-supported index type");
      m_iparm.setZero();
      m_msglvl = 0; 
      m_isInitialized = false;
    }

    ~PardisoImpl()
    {
      pardisoRelease();
    }

    inline Index cols() const { return m_size; }
    inline Index rows() const { return m_size; }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_info;
    }

    ParameterType& pardisoParameterArray()
    {
      return m_iparm;
    }

    Derived& analyzePattern(const MatrixType& matrix);

    Derived& factorize(const MatrixType& matrix);

    Derived& compute(const MatrixType& matrix);

    template<typename Rhs,typename Dest>
    void _solve_impl(const MatrixBase<Rhs> &b, MatrixBase<Dest> &dest) const;

  protected:
    void pardisoRelease()
    {
      if(m_isInitialized) 
      {
        internal::pardiso_run_selector<StorageIndex>::run(m_pt, 1, 1, m_type, -1, internal::convert_index<StorageIndex>(m_size),0, 0, 0, m_perm.data(), 0,
                                                          m_iparm.data(), m_msglvl, NULL, NULL);
        m_isInitialized = false;
      }
    }

    void pardisoInit(int type)
    {
      m_type = type;
      bool symmetric = std::abs(m_type) < 10;
      m_iparm[0] = 1;   
      m_iparm[1] = 2;   
      m_iparm[2] = 0;   
      m_iparm[3] = 0;   
      m_iparm[4] = 0;   
      m_iparm[5] = 0;   
      m_iparm[6] = 0;   
      m_iparm[7] = 2;   
      m_iparm[8] = 0;   
      m_iparm[9] = 13;  
      m_iparm[10] = symmetric ? 0 : 1; 
      m_iparm[11] = 0;  
      m_iparm[12] = symmetric ? 0 : 1;  
                                        
      m_iparm[13] = 0;  
      m_iparm[14] = 0;  
      m_iparm[15] = 0;  
      m_iparm[16] = 0;  
      m_iparm[17] = -1; 
      m_iparm[18] = -1; 
      m_iparm[19] = 0;  
      
      m_iparm[20] = 0;  
      m_iparm[26] = 0;  
      m_iparm[27] = (sizeof(RealScalar) == 4) ? 1 : 0;
      m_iparm[34] = 1;  
      m_iparm[36] = 0;  
      m_iparm[59] = 0;  
      
      memset(m_pt, 0, sizeof(m_pt));
    }

  protected:

    void manageErrorCode(Index error) const
    {
      switch(error)
      {
        case 0:
          m_info = Success;
          break;
        case -4:
        case -7:
          m_info = NumericalIssue;
          break;
        default:
          m_info = InvalidInput;
      }
    }

    mutable SparseMatrixType m_matrix;
    mutable ComputationInfo m_info;
    bool m_analysisIsOk, m_factorizationIsOk;
    StorageIndex m_type, m_msglvl;
    mutable void *m_pt[64];
    mutable ParameterType m_iparm;
    mutable IntColVectorType m_perm;
    Index m_size;
    
};

template<class Derived>
Derived& PardisoImpl<Derived>::compute(const MatrixType& a)
{
  m_size = a.rows();
  eigen_assert(a.rows() == a.cols());

  pardisoRelease();
  m_perm.setZero(m_size);
  derived().getMatrix(a);
  
  Index error;
  error = internal::pardiso_run_selector<StorageIndex>::run(m_pt, 1, 1, m_type, 12, internal::convert_index<StorageIndex>(m_size),
                                                            m_matrix.valuePtr(), m_matrix.outerIndexPtr(), m_matrix.innerIndexPtr(),
                                                            m_perm.data(), 0, m_iparm.data(), m_msglvl, NULL, NULL);
  manageErrorCode(error);
  m_analysisIsOk = true;
  m_factorizationIsOk = true;
  m_isInitialized = true;
  return derived();
}

template<class Derived>
Derived& PardisoImpl<Derived>::analyzePattern(const MatrixType& a)
{
  m_size = a.rows();
  eigen_assert(m_size == a.cols());

  pardisoRelease();
  m_perm.setZero(m_size);
  derived().getMatrix(a);
  
  Index error;
  error = internal::pardiso_run_selector<StorageIndex>::run(m_pt, 1, 1, m_type, 11, internal::convert_index<StorageIndex>(m_size),
                                                            m_matrix.valuePtr(), m_matrix.outerIndexPtr(), m_matrix.innerIndexPtr(),
                                                            m_perm.data(), 0, m_iparm.data(), m_msglvl, NULL, NULL);
  
  manageErrorCode(error);
  m_analysisIsOk = true;
  m_factorizationIsOk = false;
  m_isInitialized = true;
  return derived();
}

template<class Derived>
Derived& PardisoImpl<Derived>::factorize(const MatrixType& a)
{
  eigen_assert(m_analysisIsOk && "You must first call analyzePattern()");
  eigen_assert(m_size == a.rows() && m_size == a.cols());
  
  derived().getMatrix(a);

  Index error;
  error = internal::pardiso_run_selector<StorageIndex>::run(m_pt, 1, 1, m_type, 22, internal::convert_index<StorageIndex>(m_size),
                                                            m_matrix.valuePtr(), m_matrix.outerIndexPtr(), m_matrix.innerIndexPtr(),
                                                            m_perm.data(), 0, m_iparm.data(), m_msglvl, NULL, NULL);
  
  manageErrorCode(error);
  m_factorizationIsOk = true;
  return derived();
}

template<class Derived>
template<typename BDerived,typename XDerived>
void PardisoImpl<Derived>::_solve_impl(const MatrixBase<BDerived> &b, MatrixBase<XDerived>& x) const
{
  if(m_iparm[0] == 0) 
  {
    m_info = InvalidInput;
    return;
  }

  Index nrhs = Index(b.cols());
  eigen_assert(m_size==b.rows());
  eigen_assert(((MatrixBase<BDerived>::Flags & RowMajorBit) == 0 || nrhs == 1) && "Row-major right hand sides are not supported");
  eigen_assert(((MatrixBase<XDerived>::Flags & RowMajorBit) == 0 || nrhs == 1) && "Row-major matrices of unknowns are not supported");
  eigen_assert(((nrhs == 1) || b.outerStride() == b.rows()));

  Scalar* rhs_ptr = const_cast<Scalar*>(b.derived().data());
  Matrix<Scalar,Dynamic,Dynamic,ColMajor> tmp;

  if(rhs_ptr == x.derived().data())
  {
    tmp = b;
    rhs_ptr = tmp.data();
  }
  
  Index error;
  error = internal::pardiso_run_selector<StorageIndex>::run(m_pt, 1, 1, m_type, 33, internal::convert_index<StorageIndex>(m_size),
                                                            m_matrix.valuePtr(), m_matrix.outerIndexPtr(), m_matrix.innerIndexPtr(),
                                                            m_perm.data(), internal::convert_index<StorageIndex>(nrhs), m_iparm.data(), m_msglvl,
                                                            rhs_ptr, x.derived().data());

  manageErrorCode(error);
}

template<typename MatrixType>
class PardisoLU : public PardisoImpl< PardisoLU<MatrixType> >
{
  protected:
    typedef PardisoImpl<PardisoLU> Base;
    using Base::pardisoInit;
    using Base::m_matrix;
    friend class PardisoImpl< PardisoLU<MatrixType> >;

  public:

    typedef typename Base::Scalar Scalar;
    typedef typename Base::RealScalar RealScalar;

    using Base::compute;
    using Base::solve;

    PardisoLU()
      : Base()
    {
      pardisoInit(Base::ScalarIsComplex ? 13 : 11);
    }

    explicit PardisoLU(const MatrixType& matrix)
      : Base()
    {
      pardisoInit(Base::ScalarIsComplex ? 13 : 11);
      compute(matrix);
    }
  protected:
    void getMatrix(const MatrixType& matrix)
    {
      m_matrix = matrix;
      m_matrix.makeCompressed();
    }
};

template<typename MatrixType, int UpLo_>
class PardisoLLT : public PardisoImpl< PardisoLLT<MatrixType,UpLo_> >
{
  protected:
    typedef PardisoImpl< PardisoLLT<MatrixType,UpLo_> > Base;
    using Base::pardisoInit;
    using Base::m_matrix;
    friend class PardisoImpl< PardisoLLT<MatrixType,UpLo_> >;

  public:

    typedef typename Base::Scalar Scalar;
    typedef typename Base::RealScalar RealScalar;
    typedef typename Base::StorageIndex StorageIndex;
    enum { UpLo = UpLo_ };
    using Base::compute;

    PardisoLLT()
      : Base()
    {
      pardisoInit(Base::ScalarIsComplex ? 4 : 2);
    }

    explicit PardisoLLT(const MatrixType& matrix)
      : Base()
    {
      pardisoInit(Base::ScalarIsComplex ? 4 : 2);
      compute(matrix);
    }
    
  protected:
    
    void getMatrix(const MatrixType& matrix)
    {
      
      PermutationMatrix<Dynamic,Dynamic,StorageIndex> p_null;
      m_matrix.resize(matrix.rows(), matrix.cols());
      m_matrix.template selfadjointView<Upper>() = matrix.template selfadjointView<UpLo>().twistedBy(p_null);
      m_matrix.makeCompressed();
    }
};

template<typename MatrixType, int Options>
class PardisoLDLT : public PardisoImpl< PardisoLDLT<MatrixType,Options> >
{
  protected:
    typedef PardisoImpl< PardisoLDLT<MatrixType,Options> > Base;
    using Base::pardisoInit;
    using Base::m_matrix;
    friend class PardisoImpl< PardisoLDLT<MatrixType,Options> >;

  public:

    typedef typename Base::Scalar Scalar;
    typedef typename Base::RealScalar RealScalar;
    typedef typename Base::StorageIndex StorageIndex;
    using Base::compute;
    enum { UpLo = Options&(Upper|Lower) };

    PardisoLDLT()
      : Base()
    {
      pardisoInit(Base::ScalarIsComplex ? ( bool(Options&Symmetric) ? 6 : -4 ) : -2);
    }

    explicit PardisoLDLT(const MatrixType& matrix)
      : Base()
    {
      pardisoInit(Base::ScalarIsComplex ? ( bool(Options&Symmetric) ? 6 : -4 ) : -2);
      compute(matrix);
    }
    
    void getMatrix(const MatrixType& matrix)
    {
      
      PermutationMatrix<Dynamic,Dynamic,StorageIndex> p_null;
      m_matrix.resize(matrix.rows(), matrix.cols());
      m_matrix.template selfadjointView<Upper>() = matrix.template selfadjointView<UpLo>().twistedBy(p_null);
      m_matrix.makeCompressed();
    }
};

} 

#endif 
