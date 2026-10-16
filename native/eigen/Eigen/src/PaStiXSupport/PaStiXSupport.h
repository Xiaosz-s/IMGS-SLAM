#ifndef EIGEN_PASTIXSUPPORT_H
#define EIGEN_PASTIXSUPPORT_H

namespace Eigen { 

#if defined(DCOMPLEX)
  #define PASTIX_COMPLEX  COMPLEX
  #define PASTIX_DCOMPLEX DCOMPLEX
#else
  #define PASTIX_COMPLEX  std::complex<float>
  #define PASTIX_DCOMPLEX std::complex<double>
#endif

template<typename MatrixType_, bool IsStrSym = false> class PastixLU;
template<typename MatrixType_, int Options> class PastixLLT;
template<typename MatrixType_, int Options> class PastixLDLT;

namespace internal
{
    
  template<class Pastix> struct pastix_traits;

  template<typename MatrixType_>
  struct pastix_traits< PastixLU<MatrixType_> >
  {
    typedef MatrixType_ MatrixType;
    typedef typename MatrixType_::Scalar Scalar;
    typedef typename MatrixType_::RealScalar RealScalar;
    typedef typename MatrixType_::StorageIndex StorageIndex;
  };

  template<typename MatrixType_, int Options>
  struct pastix_traits< PastixLLT<MatrixType_,Options> >
  {
    typedef MatrixType_ MatrixType;
    typedef typename MatrixType_::Scalar Scalar;
    typedef typename MatrixType_::RealScalar RealScalar;
    typedef typename MatrixType_::StorageIndex StorageIndex;
  };

  template<typename MatrixType_, int Options>
  struct pastix_traits< PastixLDLT<MatrixType_,Options> >
  {
    typedef MatrixType_ MatrixType;
    typedef typename MatrixType_::Scalar Scalar;
    typedef typename MatrixType_::RealScalar RealScalar;
    typedef typename MatrixType_::StorageIndex StorageIndex;
  };
  
  inline void eigen_pastix(pastix_data_t **pastix_data, int pastix_comm, int n, int *ptr, int *idx, float *vals, int *perm, int * invp, float *x, int nbrhs, int *iparm, double *dparm)
  {
    if (n == 0) { ptr = NULL; idx = NULL; vals = NULL; }
    if (nbrhs == 0) {x = NULL; nbrhs=1;}
    s_pastix(pastix_data, pastix_comm, n, ptr, idx, vals, perm, invp, x, nbrhs, iparm, dparm); 
  }
  
  inline void eigen_pastix(pastix_data_t **pastix_data, int pastix_comm, int n, int *ptr, int *idx, double *vals, int *perm, int * invp, double *x, int nbrhs, int *iparm, double *dparm)
  {
    if (n == 0) { ptr = NULL; idx = NULL; vals = NULL; }
    if (nbrhs == 0) {x = NULL; nbrhs=1;}
    d_pastix(pastix_data, pastix_comm, n, ptr, idx, vals, perm, invp, x, nbrhs, iparm, dparm); 
  }
  
  inline void eigen_pastix(pastix_data_t **pastix_data, int pastix_comm, int n, int *ptr, int *idx, std::complex<float> *vals, int *perm, int * invp, std::complex<float> *x, int nbrhs, int *iparm, double *dparm)
  {
    if (n == 0) { ptr = NULL; idx = NULL; vals = NULL; }
    if (nbrhs == 0) {x = NULL; nbrhs=1;}
    c_pastix(pastix_data, pastix_comm, n, ptr, idx, reinterpret_cast<PASTIX_COMPLEX*>(vals), perm, invp, reinterpret_cast<PASTIX_COMPLEX*>(x), nbrhs, iparm, dparm); 
  }
  
  inline void eigen_pastix(pastix_data_t **pastix_data, int pastix_comm, int n, int *ptr, int *idx, std::complex<double> *vals, int *perm, int * invp, std::complex<double> *x, int nbrhs, int *iparm, double *dparm)
  {
    if (n == 0) { ptr = NULL; idx = NULL; vals = NULL; }
    if (nbrhs == 0) {x = NULL; nbrhs=1;}
    z_pastix(pastix_data, pastix_comm, n, ptr, idx, reinterpret_cast<PASTIX_DCOMPLEX*>(vals), perm, invp, reinterpret_cast<PASTIX_DCOMPLEX*>(x), nbrhs, iparm, dparm); 
  }

  template <typename MatrixType>
  void c_to_fortran_numbering (MatrixType& mat)
  {
    if ( !(mat.outerIndexPtr()[0]) ) 
    { 
      int i;
      for(i = 0; i <= mat.rows(); ++i)
        ++mat.outerIndexPtr()[i];
      for(i = 0; i < mat.nonZeros(); ++i)
        ++mat.innerIndexPtr()[i];
    }
  }

  template <typename MatrixType>
  void fortran_to_c_numbering (MatrixType& mat)
  {
    
    if ( mat.outerIndexPtr()[0] == 1 ) 
    { 
      int i;
      for(i = 0; i <= mat.rows(); ++i)
        --mat.outerIndexPtr()[i];
      for(i = 0; i < mat.nonZeros(); ++i)
        --mat.innerIndexPtr()[i];
    }
  }
}

template <class Derived>
class PastixBase : public SparseSolverBase<Derived>
{
  protected:
    typedef SparseSolverBase<Derived> Base;
    using Base::derived;
    using Base::m_isInitialized;
  public:
    using Base::_solve_impl;
    
    typedef typename internal::pastix_traits<Derived>::MatrixType MatrixType_;
    typedef MatrixType_ MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename MatrixType::StorageIndex StorageIndex;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    typedef SparseMatrix<Scalar, ColMajor> ColSpMatrix;
    enum {
      ColsAtCompileTime = MatrixType::ColsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };
    
  public:
    
    PastixBase() : m_initisOk(false), m_analysisIsOk(false), m_factorizationIsOk(false), m_pastixdata(0), m_size(0)
    {
      init();
    }
    
    ~PastixBase() 
    {
      clean();
    }
    
    template<typename Rhs,typename Dest>
    bool _solve_impl(const MatrixBase<Rhs> &b, MatrixBase<Dest> &x) const;

    Array<StorageIndex,IPARM_SIZE,1>& iparm()
    {
      return m_iparm; 
    }

    int& iparm(int idxparam)
    {
      return m_iparm(idxparam);
    }

    Array<double,DPARM_SIZE,1>& dparm()
    {
      return m_dparm; 
    }

    double& dparm(int idxparam)
    {
      return m_dparm(idxparam);
    }
    
    inline Index cols() const { return m_size; }
    inline Index rows() const { return m_size; }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_info;
    }
    
  protected:

    void init(); 

    void analyzePattern(ColSpMatrix& mat);

    void factorize(ColSpMatrix& mat);

    void clean()
    {
      eigen_assert(m_initisOk && "The Pastix structure should be allocated first"); 
      m_iparm(IPARM_START_TASK) = API_TASK_CLEAN;
      m_iparm(IPARM_END_TASK) = API_TASK_CLEAN;
      internal::eigen_pastix(&m_pastixdata, MPI_COMM_WORLD, 0, 0, 0, (Scalar*)0,
                             m_perm.data(), m_invp.data(), 0, 0, m_iparm.data(), m_dparm.data());
    }
    
    void compute(ColSpMatrix& mat);
    
    int m_initisOk; 
    int m_analysisIsOk;
    int m_factorizationIsOk;
    mutable ComputationInfo m_info; 
    mutable pastix_data_t *m_pastixdata; 
    mutable int m_comm; 
    mutable Array<int,IPARM_SIZE,1> m_iparm; 
    mutable Array<double,DPARM_SIZE,1> m_dparm; 
    mutable Matrix<StorageIndex,Dynamic,1> m_perm;  
    mutable Matrix<StorageIndex,Dynamic,1> m_invp;  
    mutable int m_size; 
}; 

template <class Derived>
void PastixBase<Derived>::init()
{
  m_size = 0; 
  m_iparm.setZero(IPARM_SIZE);
  m_dparm.setZero(DPARM_SIZE);
  
  m_iparm(IPARM_MODIFY_PARAMETER) = API_NO;
  pastix(&m_pastixdata, MPI_COMM_WORLD,
         0, 0, 0, 0,
         0, 0, 0, 1, m_iparm.data(), m_dparm.data());
  
  m_iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
  m_iparm[IPARM_VERBOSE]             = API_VERBOSE_NOT;
  m_iparm[IPARM_ORDERING]            = API_ORDER_SCOTCH;
  m_iparm[IPARM_INCOMPLETE]          = API_NO;
  m_iparm[IPARM_OOC_LIMIT]           = 2000;
  m_iparm[IPARM_RHS_MAKING]          = API_RHS_B;
  m_iparm(IPARM_MATRIX_VERIFICATION) = API_NO;
  
  m_iparm(IPARM_START_TASK) = API_TASK_INIT;
  m_iparm(IPARM_END_TASK) = API_TASK_INIT;
  internal::eigen_pastix(&m_pastixdata, MPI_COMM_WORLD, 0, 0, 0, (Scalar*)0,
                         0, 0, 0, 0, m_iparm.data(), m_dparm.data());

  if(m_iparm(IPARM_ERROR_NUMBER)) {
    m_info = InvalidInput;
    m_initisOk = false;
  }
  else { 
    m_info = Success;
    m_initisOk = true;
  }
}

template <class Derived>
void PastixBase<Derived>::compute(ColSpMatrix& mat)
{
  eigen_assert(mat.rows() == mat.cols() && "The input matrix should be squared");
  
  analyzePattern(mat);  
  factorize(mat);
  
  m_iparm(IPARM_MATRIX_VERIFICATION) = API_NO;
}

template <class Derived>
void PastixBase<Derived>::analyzePattern(ColSpMatrix& mat)
{                         
  eigen_assert(m_initisOk && "The initialization of PaSTiX failed");

  if(m_size>0)
    clean();
  
  m_size = internal::convert_index<int>(mat.rows());
  m_perm.resize(m_size);
  m_invp.resize(m_size);
  
  m_iparm(IPARM_START_TASK) = API_TASK_ORDERING;
  m_iparm(IPARM_END_TASK) = API_TASK_ANALYSE;
  internal::eigen_pastix(&m_pastixdata, MPI_COMM_WORLD, m_size, mat.outerIndexPtr(), mat.innerIndexPtr(),
               mat.valuePtr(), m_perm.data(), m_invp.data(), 0, 0, m_iparm.data(), m_dparm.data());

  if(m_iparm(IPARM_ERROR_NUMBER))
  {
    m_info = NumericalIssue;
    m_analysisIsOk = false;
  }
  else
  { 
    m_info = Success;
    m_analysisIsOk = true;
  }
}

template <class Derived>
void PastixBase<Derived>::factorize(ColSpMatrix& mat)
{

  eigen_assert(m_analysisIsOk && "The analysis phase should be called before the factorization phase");
  m_iparm(IPARM_START_TASK) = API_TASK_NUMFACT;
  m_iparm(IPARM_END_TASK) = API_TASK_NUMFACT;
  m_size = internal::convert_index<int>(mat.rows());
  
  internal::eigen_pastix(&m_pastixdata, MPI_COMM_WORLD, m_size, mat.outerIndexPtr(), mat.innerIndexPtr(),
               mat.valuePtr(), m_perm.data(), m_invp.data(), 0, 0, m_iparm.data(), m_dparm.data());

  if(m_iparm(IPARM_ERROR_NUMBER))
  {
    m_info = NumericalIssue;
    m_factorizationIsOk = false;
    m_isInitialized = false;
  }
  else
  {
    m_info = Success;
    m_factorizationIsOk = true;
    m_isInitialized = true;
  }
}

template<typename Base>
template<typename Rhs,typename Dest>
bool PastixBase<Base>::_solve_impl(const MatrixBase<Rhs> &b, MatrixBase<Dest> &x) const
{
  eigen_assert(m_isInitialized && "The matrix should be factorized first");
  EIGEN_STATIC_ASSERT((Dest::Flags&RowMajorBit)==0,
                     THIS_METHOD_IS_ONLY_FOR_COLUMN_MAJOR_MATRICES);
  int rhs = 1;
  
  x = b; 
  
  for (int i = 0; i < b.cols(); i++){
    m_iparm[IPARM_START_TASK]          = API_TASK_SOLVE;
    m_iparm[IPARM_END_TASK]            = API_TASK_REFINE;
  
    internal::eigen_pastix(&m_pastixdata, MPI_COMM_WORLD, internal::convert_index<int>(x.rows()), 0, 0, 0,
                           m_perm.data(), m_invp.data(), &x(0, i), rhs, m_iparm.data(), m_dparm.data());
  }

  m_info = m_iparm(IPARM_ERROR_NUMBER)==0 ? Success : NumericalIssue;
  
  return m_iparm(IPARM_ERROR_NUMBER)==0;
}

template<typename MatrixType_, bool IsStrSym>
class PastixLU : public PastixBase< PastixLU<MatrixType_> >
{
  public:
    typedef MatrixType_ MatrixType;
    typedef PastixBase<PastixLU<MatrixType> > Base;
    typedef typename Base::ColSpMatrix ColSpMatrix;
    typedef typename MatrixType::StorageIndex StorageIndex;
    
  public:
    PastixLU() : Base()
    {
      init();
    }
    
    explicit PastixLU(const MatrixType& matrix):Base()
    {
      init();
      compute(matrix);
    }

    void compute (const MatrixType& matrix)
    {
      m_structureIsUptodate = false;
      ColSpMatrix temp;
      grabMatrix(matrix, temp);
      Base::compute(temp);
    }

    void analyzePattern(const MatrixType& matrix)
    {
      m_structureIsUptodate = false;
      ColSpMatrix temp;
      grabMatrix(matrix, temp);
      Base::analyzePattern(temp);
    }

    void factorize(const MatrixType& matrix)
    {
      ColSpMatrix temp;
      grabMatrix(matrix, temp);
      Base::factorize(temp);
    }
  protected:
    
    void init()
    {
      m_structureIsUptodate = false;
      m_iparm(IPARM_SYM) = API_SYM_NO;
      m_iparm(IPARM_FACTORIZATION) = API_FACT_LU;
    }
    
    void grabMatrix(const MatrixType& matrix, ColSpMatrix& out)
    {
      if(IsStrSym)
        out = matrix;
      else
      {
        if(!m_structureIsUptodate)
        {
          
          m_transposedStructure = matrix.transpose();

          for (Index j=0; j<m_transposedStructure.outerSize(); ++j) 
            for(typename ColSpMatrix::InnerIterator it(m_transposedStructure, j); it; ++it)
              it.valueRef() = 0.0;

          m_structureIsUptodate = true;
        }
        
        out = m_transposedStructure + matrix;
      }
      internal::c_to_fortran_numbering(out);
    }
    
    using Base::m_iparm;
    using Base::m_dparm;
    
    ColSpMatrix m_transposedStructure;
    bool m_structureIsUptodate;
};

template<typename MatrixType_, int UpLo_>
class PastixLLT : public PastixBase< PastixLLT<MatrixType_, UpLo_> >
{
  public:
    typedef MatrixType_ MatrixType;
    typedef PastixBase<PastixLLT<MatrixType, UpLo_> > Base;
    typedef typename Base::ColSpMatrix ColSpMatrix;
    
  public:
    enum { UpLo = UpLo_ };
    PastixLLT() : Base()
    {
      init();
    }
    
    explicit PastixLLT(const MatrixType& matrix):Base()
    {
      init();
      compute(matrix);
    }

    void compute (const MatrixType& matrix)
    {
      ColSpMatrix temp;
      grabMatrix(matrix, temp);
      Base::compute(temp);
    }

    void analyzePattern(const MatrixType& matrix)
    {
      ColSpMatrix temp;
      grabMatrix(matrix, temp);
      Base::analyzePattern(temp);
    }

    void factorize(const MatrixType& matrix)
    {
      ColSpMatrix temp;
      grabMatrix(matrix, temp);
      Base::factorize(temp);
    }
  protected:
    using Base::m_iparm;
    
    void init()
    {
      m_iparm(IPARM_SYM) = API_SYM_YES;
      m_iparm(IPARM_FACTORIZATION) = API_FACT_LLT;
    }
    
    void grabMatrix(const MatrixType& matrix, ColSpMatrix& out)
    {
      out.resize(matrix.rows(), matrix.cols());
      
      out.template selfadjointView<Lower>() = matrix.template selfadjointView<UpLo>();
      internal::c_to_fortran_numbering(out);
    }
};

template<typename MatrixType_, int UpLo_>
class PastixLDLT : public PastixBase< PastixLDLT<MatrixType_, UpLo_> >
{
  public:
    typedef MatrixType_ MatrixType;
    typedef PastixBase<PastixLDLT<MatrixType, UpLo_> > Base;
    typedef typename Base::ColSpMatrix ColSpMatrix;
    
  public:
    enum { UpLo = UpLo_ };
    PastixLDLT():Base()
    {
      init();
    }
    
    explicit PastixLDLT(const MatrixType& matrix):Base()
    {
      init();
      compute(matrix);
    }

    void compute (const MatrixType& matrix)
    {
      ColSpMatrix temp;
      grabMatrix(matrix, temp);
      Base::compute(temp);
    }

    void analyzePattern(const MatrixType& matrix)
    { 
      ColSpMatrix temp;
      grabMatrix(matrix, temp);
      Base::analyzePattern(temp);
    }

    void factorize(const MatrixType& matrix)
    {
      ColSpMatrix temp;
      grabMatrix(matrix, temp);
      Base::factorize(temp);
    }

  protected:
    using Base::m_iparm;
    
    void init()
    {
      m_iparm(IPARM_SYM) = API_SYM_YES;
      m_iparm(IPARM_FACTORIZATION) = API_FACT_LDLT;
    }
    
    void grabMatrix(const MatrixType& matrix, ColSpMatrix& out)
    {
      
      out.resize(matrix.rows(), matrix.cols());
      out.template selfadjointView<Lower>() = matrix.template selfadjointView<UpLo>();
      internal::c_to_fortran_numbering(out);
    }
};

} 

#endif
