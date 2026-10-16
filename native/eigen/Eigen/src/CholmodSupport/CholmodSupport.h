#ifndef EIGEN_CHOLMODSUPPORT_H
#define EIGEN_CHOLMODSUPPORT_H

namespace Eigen {

namespace internal {

template<typename Scalar> struct cholmod_configure_matrix;

template<> struct cholmod_configure_matrix<double> {
  template<typename CholmodType>
  static void run(CholmodType& mat) {
    mat.xtype = CHOLMOD_REAL;
    mat.dtype = CHOLMOD_DOUBLE;
  }
};

template<> struct cholmod_configure_matrix<std::complex<double> > {
  template<typename CholmodType>
  static void run(CholmodType& mat) {
    mat.xtype = CHOLMOD_COMPLEX;
    mat.dtype = CHOLMOD_DOUBLE;
  }
};

} 

template<typename Scalar_, int Options_, typename StorageIndex_>
cholmod_sparse viewAsCholmod(Ref<SparseMatrix<Scalar_,Options_,StorageIndex_> > mat)
{
  cholmod_sparse res;
  res.nzmax   = mat.nonZeros();
  res.nrow    = mat.rows();
  res.ncol    = mat.cols();
  res.p       = mat.outerIndexPtr();
  res.i       = mat.innerIndexPtr();
  res.x       = mat.valuePtr();
  res.z       = 0;
  res.sorted  = 1;
  if(mat.isCompressed())
  {
    res.packed  = 1;
    res.nz = 0;
  }
  else
  {
    res.packed  = 0;
    res.nz = mat.innerNonZeroPtr();
  }

  res.dtype   = 0;
  res.stype   = -1;

  if (internal::is_same<StorageIndex_,int>::value)
  {
    res.itype = CHOLMOD_INT;
  }
  else if (internal::is_same<StorageIndex_,SuiteSparse_long>::value)
  {
    res.itype = CHOLMOD_LONG;
  }
  else
  {
    eigen_assert(false && "Index type not supported yet");
  }

  internal::cholmod_configure_matrix<Scalar_>::run(res);

  res.stype = 0;

  return res;
}

template<typename Scalar_, int Options_, typename Index_>
const cholmod_sparse viewAsCholmod(const SparseMatrix<Scalar_,Options_,Index_>& mat)
{
  cholmod_sparse res = viewAsCholmod(Ref<SparseMatrix<Scalar_,Options_,Index_> >(mat.const_cast_derived()));
  return res;
}

template<typename Scalar_, int Options_, typename Index_>
const cholmod_sparse viewAsCholmod(const SparseVector<Scalar_,Options_,Index_>& mat)
{
  cholmod_sparse res = viewAsCholmod(Ref<SparseMatrix<Scalar_,Options_,Index_> >(mat.const_cast_derived()));
  return res;
}

template<typename Scalar_, int Options_, typename Index_, unsigned int UpLo>
cholmod_sparse viewAsCholmod(const SparseSelfAdjointView<const SparseMatrix<Scalar_,Options_,Index_>, UpLo>& mat)
{
  cholmod_sparse res = viewAsCholmod(Ref<SparseMatrix<Scalar_,Options_,Index_> >(mat.matrix().const_cast_derived()));

  if(UpLo==Upper) res.stype =  1;
  if(UpLo==Lower) res.stype = -1;
  
  EIGEN_STATIC_ASSERT((Options_ & RowMajorBit) == 0 || NumTraits<Scalar_>::IsComplex == 0, THIS_METHOD_IS_ONLY_FOR_COLUMN_MAJOR_MATRICES);
  if(Options_ & RowMajorBit) res.stype *=-1;

  return res;
}

template<typename Derived>
cholmod_dense viewAsCholmod(MatrixBase<Derived>& mat)
{
  EIGEN_STATIC_ASSERT((internal::traits<Derived>::Flags&RowMajorBit)==0,THIS_METHOD_IS_ONLY_FOR_COLUMN_MAJOR_MATRICES);
  typedef typename Derived::Scalar Scalar;

  cholmod_dense res;
  res.nrow   = mat.rows();
  res.ncol   = mat.cols();
  res.nzmax  = res.nrow * res.ncol;
  res.d      = Derived::IsVectorAtCompileTime ? mat.derived().size() : mat.derived().outerStride();
  res.x      = (void*)(mat.derived().data());
  res.z      = 0;

  internal::cholmod_configure_matrix<Scalar>::run(res);

  return res;
}

template<typename Scalar, int Flags, typename StorageIndex>
MappedSparseMatrix<Scalar,Flags,StorageIndex> viewAsEigen(cholmod_sparse& cm)
{
  return MappedSparseMatrix<Scalar,Flags,StorageIndex>
         (cm.nrow, cm.ncol, static_cast<StorageIndex*>(cm.p)[cm.ncol],
          static_cast<StorageIndex*>(cm.p), static_cast<StorageIndex*>(cm.i),static_cast<Scalar*>(cm.x) );
}

namespace internal {

#define EIGEN_CHOLMOD_SPECIALIZE0(ret, name) \
    template<typename StorageIndex_> inline ret cm_ ## name       (cholmod_common &Common) { return cholmod_ ## name   (&Common); } \
    template<>                       inline ret cm_ ## name<SuiteSparse_long> (cholmod_common &Common) { return cholmod_l_ ## name (&Common); }

#define EIGEN_CHOLMOD_SPECIALIZE1(ret, name, t1, a1) \
    template<typename StorageIndex_> inline ret cm_ ## name       (t1& a1, cholmod_common &Common) { return cholmod_ ## name   (&a1, &Common); } \
    template<>                       inline ret cm_ ## name<SuiteSparse_long> (t1& a1, cholmod_common &Common) { return cholmod_l_ ## name (&a1, &Common); }

EIGEN_CHOLMOD_SPECIALIZE0(int, start)
EIGEN_CHOLMOD_SPECIALIZE0(int, finish)

EIGEN_CHOLMOD_SPECIALIZE1(int, free_factor, cholmod_factor*, L)
EIGEN_CHOLMOD_SPECIALIZE1(int, free_dense,  cholmod_dense*,  X)
EIGEN_CHOLMOD_SPECIALIZE1(int, free_sparse, cholmod_sparse*, A)

EIGEN_CHOLMOD_SPECIALIZE1(cholmod_factor*, analyze, cholmod_sparse, A)

template<typename StorageIndex_> inline cholmod_dense*  cm_solve         (int sys, cholmod_factor& L, cholmod_dense&  B, cholmod_common &Common) { return cholmod_solve     (sys, &L, &B, &Common); }
template<>                       inline cholmod_dense*  cm_solve<SuiteSparse_long>   (int sys, cholmod_factor& L, cholmod_dense&  B, cholmod_common &Common) { return cholmod_l_solve   (sys, &L, &B, &Common); }

template<typename StorageIndex_> inline cholmod_sparse* cm_spsolve       (int sys, cholmod_factor& L, cholmod_sparse& B, cholmod_common &Common) { return cholmod_spsolve   (sys, &L, &B, &Common); }
template<>                       inline cholmod_sparse* cm_spsolve<SuiteSparse_long> (int sys, cholmod_factor& L, cholmod_sparse& B, cholmod_common &Common) { return cholmod_l_spsolve (sys, &L, &B, &Common); }

template<typename StorageIndex_>
inline int  cm_factorize_p       (cholmod_sparse*  A, double beta[2], StorageIndex_* fset, std::size_t fsize, cholmod_factor* L, cholmod_common &Common) { return cholmod_factorize_p   (A, beta, fset, fsize, L, &Common); }
template<>
inline int  cm_factorize_p<SuiteSparse_long> (cholmod_sparse*  A, double beta[2], SuiteSparse_long* fset,          std::size_t fsize, cholmod_factor* L, cholmod_common &Common) { return cholmod_l_factorize_p (A, beta, fset, fsize, L, &Common); }

#undef EIGEN_CHOLMOD_SPECIALIZE0
#undef EIGEN_CHOLMOD_SPECIALIZE1

}  

enum CholmodMode {
  CholmodAuto, CholmodSimplicialLLt, CholmodSupernodalLLt, CholmodLDLt
};

template<typename MatrixType_, int UpLo_, typename Derived>
class CholmodBase : public SparseSolverBase<Derived>
{
  protected:
    typedef SparseSolverBase<Derived> Base;
    using Base::derived;
    using Base::m_isInitialized;
  public:
    typedef MatrixType_ MatrixType;
    enum { UpLo = UpLo_ };
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef MatrixType CholMatrixType;
    typedef typename MatrixType::StorageIndex StorageIndex;
    enum {
      ColsAtCompileTime = MatrixType::ColsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };

  public:

    CholmodBase()
      : m_cholmodFactor(0), m_info(Success), m_factorizationIsOk(false), m_analysisIsOk(false)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<double,RealScalar>::value), CHOLMOD_SUPPORTS_DOUBLE_PRECISION_ONLY);
      m_shiftOffset[0] = m_shiftOffset[1] = 0.0;
      internal::cm_start<StorageIndex>(m_cholmod);
    }

    explicit CholmodBase(const MatrixType& matrix)
      : m_cholmodFactor(0), m_info(Success), m_factorizationIsOk(false), m_analysisIsOk(false)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<double,RealScalar>::value), CHOLMOD_SUPPORTS_DOUBLE_PRECISION_ONLY);
      m_shiftOffset[0] = m_shiftOffset[1] = 0.0;
      internal::cm_start<StorageIndex>(m_cholmod);
      compute(matrix);
    }

    ~CholmodBase()
    {
      if(m_cholmodFactor)
        internal::cm_free_factor<StorageIndex>(m_cholmodFactor, m_cholmod);
      internal::cm_finish<StorageIndex>(m_cholmod);
    }

    inline StorageIndex cols() const { return internal::convert_index<StorageIndex, Index>(m_cholmodFactor->n); }
    inline StorageIndex rows() const { return internal::convert_index<StorageIndex, Index>(m_cholmodFactor->n); }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_info;
    }

    Derived& compute(const MatrixType& matrix)
    {
      analyzePattern(matrix);
      factorize(matrix);
      return derived();
    }

    void analyzePattern(const MatrixType& matrix)
    {
      if(m_cholmodFactor)
      {
        internal::cm_free_factor<StorageIndex>(m_cholmodFactor, m_cholmod);
        m_cholmodFactor = 0;
      }
      cholmod_sparse A = viewAsCholmod(matrix.template selfadjointView<UpLo>());
      m_cholmodFactor = internal::cm_analyze<StorageIndex>(A, m_cholmod);

      this->m_isInitialized = true;
      this->m_info = Success;
      m_analysisIsOk = true;
      m_factorizationIsOk = false;
    }

    void factorize(const MatrixType& matrix)
    {
      eigen_assert(m_analysisIsOk && "You must first call analyzePattern()");
      cholmod_sparse A = viewAsCholmod(matrix.template selfadjointView<UpLo>());
      internal::cm_factorize_p<StorageIndex>(&A, m_shiftOffset, 0, 0, m_cholmodFactor, m_cholmod);

      this->m_info = (m_cholmodFactor->minor == m_cholmodFactor->n ? Success : NumericalIssue);
      m_factorizationIsOk = true;
    }

    cholmod_common& cholmod() { return m_cholmod; }

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    
    template<typename Rhs,typename Dest>
    void _solve_impl(const MatrixBase<Rhs> &b, MatrixBase<Dest> &dest) const
    {
      eigen_assert(m_factorizationIsOk && "The decomposition is not in a valid state for solving, you must first call either compute() or symbolic()/numeric()");
      const Index size = m_cholmodFactor->n;
      EIGEN_UNUSED_VARIABLE(size);
      eigen_assert(size==b.rows());

      Ref<const Matrix<typename Rhs::Scalar,Dynamic,Dynamic,ColMajor> > b_ref(b.derived());

      cholmod_dense b_cd = viewAsCholmod(b_ref);
      cholmod_dense* x_cd = internal::cm_solve<StorageIndex>(CHOLMOD_A, *m_cholmodFactor, b_cd, m_cholmod);
      if(!x_cd)
      {
        this->m_info = NumericalIssue;
        return;
      }

      dest = Matrix<Scalar,Dest::RowsAtCompileTime,Dest::ColsAtCompileTime>::Map(reinterpret_cast<Scalar*>(x_cd->x),b.rows(),b.cols());
      internal::cm_free_dense<StorageIndex>(x_cd, m_cholmod);
    }

    template<typename RhsDerived, typename DestDerived>
    void _solve_impl(const SparseMatrixBase<RhsDerived> &b, SparseMatrixBase<DestDerived> &dest) const
    {
      eigen_assert(m_factorizationIsOk && "The decomposition is not in a valid state for solving, you must first call either compute() or symbolic()/numeric()");
      const Index size = m_cholmodFactor->n;
      EIGEN_UNUSED_VARIABLE(size);
      eigen_assert(size==b.rows());

      Ref<SparseMatrix<typename RhsDerived::Scalar,ColMajor,typename RhsDerived::StorageIndex> > b_ref(b.const_cast_derived());
      cholmod_sparse b_cs = viewAsCholmod(b_ref);
      cholmod_sparse* x_cs = internal::cm_spsolve<StorageIndex>(CHOLMOD_A, *m_cholmodFactor, b_cs, m_cholmod);
      if(!x_cs)
      {
        this->m_info = NumericalIssue;
        return;
      }

      dest.derived() = viewAsEigen<typename DestDerived::Scalar,ColMajor,typename DestDerived::StorageIndex>(*x_cs);
      internal::cm_free_sparse<StorageIndex>(x_cs, m_cholmod);
    }
    #endif 

    Derived& setShift(const RealScalar& offset)
    {
      m_shiftOffset[0] = double(offset);
      return derived();
    }

    Scalar determinant() const
    {
      using std::exp;
      return exp(logDeterminant());
    }

    Scalar logDeterminant() const
    {
      using std::log;
      using numext::real;
      eigen_assert(m_factorizationIsOk && "The decomposition is not in a valid state for solving, you must first call either compute() or symbolic()/numeric()");

      RealScalar logDet = 0;
      Scalar *x = static_cast<Scalar*>(m_cholmodFactor->x);
      if (m_cholmodFactor->is_super)
      {

        StorageIndex *super = static_cast<StorageIndex*>(m_cholmodFactor->super);
        
        StorageIndex *pi = static_cast<StorageIndex*>(m_cholmodFactor->pi);
        
        StorageIndex *px = static_cast<StorageIndex*>(m_cholmodFactor->px);

        Index nb_super_nodes = m_cholmodFactor->nsuper;
        for (Index k=0; k < nb_super_nodes; ++k)
        {
          StorageIndex ncols = super[k + 1] - super[k];
          StorageIndex nrows = pi[k + 1] - pi[k];

          Map<const Array<Scalar,1,Dynamic>, 0, InnerStride<> > sk(x + px[k], ncols, InnerStride<>(nrows+1));
          logDet += sk.real().log().sum();
        }
      }
      else
      {
        
        StorageIndex *p = static_cast<StorageIndex*>(m_cholmodFactor->p);
        Index size = m_cholmodFactor->n;
        for (Index k=0; k<size; ++k)
          logDet += log(real( x[p[k]] ));
      }
      if (m_cholmodFactor->is_ll)
        logDet *= 2.0;
      return logDet;
    };

    template<typename Stream>
    void dumpMemory(Stream& )
    {}

  protected:
    mutable cholmod_common m_cholmod;
    cholmod_factor* m_cholmodFactor;
    double m_shiftOffset[2];
    mutable ComputationInfo m_info;
    int m_factorizationIsOk;
    int m_analysisIsOk;
};

template<typename MatrixType_, int UpLo_ = Lower>
class CholmodSimplicialLLT : public CholmodBase<MatrixType_, UpLo_, CholmodSimplicialLLT<MatrixType_, UpLo_> >
{
    typedef CholmodBase<MatrixType_, UpLo_, CholmodSimplicialLLT> Base;
    using Base::m_cholmod;

  public:

    typedef MatrixType_ MatrixType;

    CholmodSimplicialLLT() : Base() { init(); }

    CholmodSimplicialLLT(const MatrixType& matrix) : Base()
    {
      init();
      this->compute(matrix);
    }

    ~CholmodSimplicialLLT() {}
  protected:
    void init()
    {
      m_cholmod.final_asis = 0;
      m_cholmod.supernodal = CHOLMOD_SIMPLICIAL;
      m_cholmod.final_ll = 1;
    }
};

template<typename MatrixType_, int UpLo_ = Lower>
class CholmodSimplicialLDLT : public CholmodBase<MatrixType_, UpLo_, CholmodSimplicialLDLT<MatrixType_, UpLo_> >
{
    typedef CholmodBase<MatrixType_, UpLo_, CholmodSimplicialLDLT> Base;
    using Base::m_cholmod;

  public:

    typedef MatrixType_ MatrixType;

    CholmodSimplicialLDLT() : Base() { init(); }

    CholmodSimplicialLDLT(const MatrixType& matrix) : Base()
    {
      init();
      this->compute(matrix);
    }

    ~CholmodSimplicialLDLT() {}
  protected:
    void init()
    {
      m_cholmod.final_asis = 1;
      m_cholmod.supernodal = CHOLMOD_SIMPLICIAL;
    }
};

template<typename MatrixType_, int UpLo_ = Lower>
class CholmodSupernodalLLT : public CholmodBase<MatrixType_, UpLo_, CholmodSupernodalLLT<MatrixType_, UpLo_> >
{
    typedef CholmodBase<MatrixType_, UpLo_, CholmodSupernodalLLT> Base;
    using Base::m_cholmod;

  public:

    typedef MatrixType_ MatrixType;

    CholmodSupernodalLLT() : Base() { init(); }

    CholmodSupernodalLLT(const MatrixType& matrix) : Base()
    {
      init();
      this->compute(matrix);
    }

    ~CholmodSupernodalLLT() {}
  protected:
    void init()
    {
      m_cholmod.final_asis = 1;
      m_cholmod.supernodal = CHOLMOD_SUPERNODAL;
    }
};

template<typename MatrixType_, int UpLo_ = Lower>
class CholmodDecomposition : public CholmodBase<MatrixType_, UpLo_, CholmodDecomposition<MatrixType_, UpLo_> >
{
    typedef CholmodBase<MatrixType_, UpLo_, CholmodDecomposition> Base;
    using Base::m_cholmod;

  public:

    typedef MatrixType_ MatrixType;

    CholmodDecomposition() : Base() { init(); }

    CholmodDecomposition(const MatrixType& matrix) : Base()
    {
      init();
      this->compute(matrix);
    }

    ~CholmodDecomposition() {}

    void setMode(CholmodMode mode)
    {
      switch(mode)
      {
        case CholmodAuto:
          m_cholmod.final_asis = 1;
          m_cholmod.supernodal = CHOLMOD_AUTO;
          break;
        case CholmodSimplicialLLt:
          m_cholmod.final_asis = 0;
          m_cholmod.supernodal = CHOLMOD_SIMPLICIAL;
          m_cholmod.final_ll = 1;
          break;
        case CholmodSupernodalLLt:
          m_cholmod.final_asis = 1;
          m_cholmod.supernodal = CHOLMOD_SUPERNODAL;
          break;
        case CholmodLDLt:
          m_cholmod.final_asis = 1;
          m_cholmod.supernodal = CHOLMOD_SIMPLICIAL;
          break;
        default:
          break;
      }
    }
  protected:
    void init()
    {
      m_cholmod.final_asis = 1;
      m_cholmod.supernodal = CHOLMOD_AUTO;
    }
};

} 

#endif 
