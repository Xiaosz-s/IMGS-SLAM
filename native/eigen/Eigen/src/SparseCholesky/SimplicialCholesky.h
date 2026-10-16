#ifndef EIGEN_SIMPLICIAL_CHOLESKY_H
#define EIGEN_SIMPLICIAL_CHOLESKY_H

namespace Eigen { 

enum SimplicialCholeskyMode {
  SimplicialCholeskyLLT,
  SimplicialCholeskyLDLT
};

namespace internal {
  template<typename CholMatrixType, typename InputMatrixType>
  struct simplicial_cholesky_grab_input {
    typedef CholMatrixType const * ConstCholMatrixPtr;
    static void run(const InputMatrixType& input, ConstCholMatrixPtr &pmat, CholMatrixType &tmp)
    {
      tmp = input;
      pmat = &tmp;
    }
  };
  
  template<typename MatrixType>
  struct simplicial_cholesky_grab_input<MatrixType,MatrixType> {
    typedef MatrixType const * ConstMatrixPtr;
    static void run(const MatrixType& input, ConstMatrixPtr &pmat, MatrixType &)
    {
      pmat = &input;
    }
  };
} 

template<typename Derived>
class SimplicialCholeskyBase : public SparseSolverBase<Derived>
{
    typedef SparseSolverBase<Derived> Base;
    using Base::m_isInitialized;
    
  public:
    typedef typename internal::traits<Derived>::MatrixType MatrixType;
    typedef typename internal::traits<Derived>::OrderingType OrderingType;
    enum { UpLo = internal::traits<Derived>::UpLo };
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename MatrixType::StorageIndex StorageIndex;
    typedef SparseMatrix<Scalar,ColMajor,StorageIndex> CholMatrixType;
    typedef CholMatrixType const * ConstCholMatrixPtr;
    typedef Matrix<Scalar,Dynamic,1> VectorType;
    typedef Matrix<StorageIndex,Dynamic,1> VectorI;

    enum {
      ColsAtCompileTime = MatrixType::ColsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };

  public:
    
    using Base::derived;

    SimplicialCholeskyBase()
      : m_info(Success),
        m_factorizationIsOk(false),
        m_analysisIsOk(false),
        m_shiftOffset(0),
        m_shiftScale(1)
    {}

    explicit SimplicialCholeskyBase(const MatrixType& matrix)
      : m_info(Success),
        m_factorizationIsOk(false),
        m_analysisIsOk(false),
        m_shiftOffset(0),
        m_shiftScale(1)
    {
      derived().compute(matrix);
    }

    ~SimplicialCholeskyBase()
    {
    }

    Derived& derived() { return *static_cast<Derived*>(this); }
    const Derived& derived() const { return *static_cast<const Derived*>(this); }
    
    inline Index cols() const { return m_matrix.cols(); }
    inline Index rows() const { return m_matrix.rows(); }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_info;
    }

    const PermutationMatrix<Dynamic,Dynamic,StorageIndex>& permutationP() const
    { return m_P; }

    const PermutationMatrix<Dynamic,Dynamic,StorageIndex>& permutationPinv() const
    { return m_Pinv; }

    Derived& setShift(const RealScalar& offset, const RealScalar& scale = 1)
    {
      m_shiftOffset = offset;
      m_shiftScale = scale;
      return derived();
    }

#ifndef EIGEN_PARSED_BY_DOXYGEN
    
    template<typename Stream>
    void dumpMemory(Stream& s)
    {
      int total = 0;
      s << "  L:        " << ((total+=(m_matrix.cols()+1) * sizeof(int) + m_matrix.nonZeros()*(sizeof(int)+sizeof(Scalar))) >> 20) << "Mb" << "\n";
      s << "  diag:     " << ((total+=m_diag.size() * sizeof(Scalar)) >> 20) << "Mb" << "\n";
      s << "  tree:     " << ((total+=m_parent.size() * sizeof(int)) >> 20) << "Mb" << "\n";
      s << "  nonzeros: " << ((total+=m_nonZerosPerCol.size() * sizeof(int)) >> 20) << "Mb" << "\n";
      s << "  perm:     " << ((total+=m_P.size() * sizeof(int)) >> 20) << "Mb" << "\n";
      s << "  perm^-1:  " << ((total+=m_Pinv.size() * sizeof(int)) >> 20) << "Mb" << "\n";
      s << "  TOTAL:    " << (total>> 20) << "Mb" << "\n";
    }

    template<typename Rhs,typename Dest>
    void _solve_impl(const MatrixBase<Rhs> &b, MatrixBase<Dest> &dest) const
    {
      eigen_assert(m_factorizationIsOk && "The decomposition is not in a valid state for solving, you must first call either compute() or symbolic()/numeric()");
      eigen_assert(m_matrix.rows()==b.rows());

      if(m_info!=Success)
        return;

      if(m_P.size()>0)
        dest = m_P * b;
      else
        dest = b;

      if(m_matrix.nonZeros()>0) 
        derived().matrixL().solveInPlace(dest);

      if(m_diag.size()>0)
        dest = m_diag.asDiagonal().inverse() * dest;

      if (m_matrix.nonZeros()>0) 
        derived().matrixU().solveInPlace(dest);

      if(m_P.size()>0)
        dest = m_Pinv * dest;
    }
    
    template<typename Rhs,typename Dest>
    void _solve_impl(const SparseMatrixBase<Rhs> &b, SparseMatrixBase<Dest> &dest) const
    {
      internal::solve_sparse_through_dense_panels(derived(), b, dest);
    }

#endif 

  protected:

    template<bool DoLDLT>
    void compute(const MatrixType& matrix)
    {
      eigen_assert(matrix.rows()==matrix.cols());
      Index size = matrix.cols();
      CholMatrixType tmp(size,size);
      ConstCholMatrixPtr pmat;
      ordering(matrix, pmat, tmp);
      analyzePattern_preordered(*pmat, DoLDLT);
      factorize_preordered<DoLDLT>(*pmat);
    }
    
    template<bool DoLDLT>
    void factorize(const MatrixType& a)
    {
      eigen_assert(a.rows()==a.cols());
      Index size = a.cols();
      CholMatrixType tmp(size,size);
      ConstCholMatrixPtr pmat;
      
      if(m_P.size() == 0 && (int(UpLo) & int(Upper)) == Upper)
      {
        
        internal::simplicial_cholesky_grab_input<CholMatrixType,MatrixType>::run(a, pmat, tmp);
      }
      else
      {
        tmp.template selfadjointView<Upper>() = a.template selfadjointView<UpLo>().twistedBy(m_P);
        pmat = &tmp;
      }
      
      factorize_preordered<DoLDLT>(*pmat);
    }

    template<bool DoLDLT>
    void factorize_preordered(const CholMatrixType& a);

    void analyzePattern(const MatrixType& a, bool doLDLT)
    {
      eigen_assert(a.rows()==a.cols());
      Index size = a.cols();
      CholMatrixType tmp(size,size);
      ConstCholMatrixPtr pmat;
      ordering(a, pmat, tmp);
      analyzePattern_preordered(*pmat,doLDLT);
    }
    void analyzePattern_preordered(const CholMatrixType& a, bool doLDLT);
    
    void ordering(const MatrixType& a, ConstCholMatrixPtr &pmat, CholMatrixType& ap);

    struct keep_diag {
      inline bool operator() (const Index& row, const Index& col, const Scalar&) const
      {
        return row!=col;
      }
    };

    mutable ComputationInfo m_info;
    bool m_factorizationIsOk;
    bool m_analysisIsOk;
    
    CholMatrixType m_matrix;
    VectorType m_diag;                                
    VectorI m_parent;                                 
    VectorI m_nonZerosPerCol;
    PermutationMatrix<Dynamic,Dynamic,StorageIndex> m_P;     
    PermutationMatrix<Dynamic,Dynamic,StorageIndex> m_Pinv;  

    RealScalar m_shiftOffset;
    RealScalar m_shiftScale;
};

template<typename MatrixType_, int UpLo_ = Lower, typename Ordering_ = AMDOrdering<typename MatrixType_::StorageIndex> > class SimplicialLLT;
template<typename MatrixType_, int UpLo_ = Lower, typename Ordering_ = AMDOrdering<typename MatrixType_::StorageIndex> > class SimplicialLDLT;
template<typename MatrixType_, int UpLo_ = Lower, typename Ordering_ = AMDOrdering<typename MatrixType_::StorageIndex> > class SimplicialCholesky;

namespace internal {

template<typename MatrixType_, int UpLo_, typename Ordering_> struct traits<SimplicialLLT<MatrixType_,UpLo_,Ordering_> >
{
  typedef MatrixType_ MatrixType;
  typedef Ordering_ OrderingType;
  enum { UpLo = UpLo_ };
  typedef typename MatrixType::Scalar                         Scalar;
  typedef typename MatrixType::StorageIndex                   StorageIndex;
  typedef SparseMatrix<Scalar, ColMajor, StorageIndex>        CholMatrixType;
  typedef TriangularView<const CholMatrixType, Eigen::Lower>  MatrixL;
  typedef TriangularView<const typename CholMatrixType::AdjointReturnType, Eigen::Upper>   MatrixU;
  static inline MatrixL getL(const CholMatrixType& m) { return MatrixL(m); }
  static inline MatrixU getU(const CholMatrixType& m) { return MatrixU(m.adjoint()); }
};

template<typename MatrixType_,int UpLo_, typename Ordering_> struct traits<SimplicialLDLT<MatrixType_,UpLo_,Ordering_> >
{
  typedef MatrixType_ MatrixType;
  typedef Ordering_ OrderingType;
  enum { UpLo = UpLo_ };
  typedef typename MatrixType::Scalar                             Scalar;
  typedef typename MatrixType::StorageIndex                       StorageIndex;
  typedef SparseMatrix<Scalar, ColMajor, StorageIndex>            CholMatrixType;
  typedef TriangularView<const CholMatrixType, Eigen::UnitLower>  MatrixL;
  typedef TriangularView<const typename CholMatrixType::AdjointReturnType, Eigen::UnitUpper> MatrixU;
  static inline MatrixL getL(const CholMatrixType& m) { return MatrixL(m); }
  static inline MatrixU getU(const CholMatrixType& m) { return MatrixU(m.adjoint()); }
};

template<typename MatrixType_, int UpLo_, typename Ordering_> struct traits<SimplicialCholesky<MatrixType_,UpLo_,Ordering_> >
{
  typedef MatrixType_ MatrixType;
  typedef Ordering_ OrderingType;
  enum { UpLo = UpLo_ };
};

}

template<typename MatrixType_, int UpLo_, typename Ordering_>
    class SimplicialLLT : public SimplicialCholeskyBase<SimplicialLLT<MatrixType_,UpLo_,Ordering_> >
{
public:
    typedef MatrixType_ MatrixType;
    enum { UpLo = UpLo_ };
    typedef SimplicialCholeskyBase<SimplicialLLT> Base;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename MatrixType::StorageIndex StorageIndex;
    typedef SparseMatrix<Scalar,ColMajor,Index> CholMatrixType;
    typedef Matrix<Scalar,Dynamic,1> VectorType;
    typedef internal::traits<SimplicialLLT> Traits;
    typedef typename Traits::MatrixL  MatrixL;
    typedef typename Traits::MatrixU  MatrixU;
public:
    
    SimplicialLLT() : Base() {}
    
    explicit SimplicialLLT(const MatrixType& matrix)
        : Base(matrix) {}

    inline const MatrixL matrixL() const {
        eigen_assert(Base::m_factorizationIsOk && "Simplicial LLT not factorized");
        return Traits::getL(Base::m_matrix);
    }

    inline const MatrixU matrixU() const {
        eigen_assert(Base::m_factorizationIsOk && "Simplicial LLT not factorized");
        return Traits::getU(Base::m_matrix);
    }

    SimplicialLLT& compute(const MatrixType& matrix)
    {
      Base::template compute<false>(matrix);
      return *this;
    }

    void analyzePattern(const MatrixType& a)
    {
      Base::analyzePattern(a, false);
    }

    void factorize(const MatrixType& a)
    {
      Base::template factorize<false>(a);
    }

    Scalar determinant() const
    {
      Scalar detL = Base::m_matrix.diagonal().prod();
      return numext::abs2(detL);
    }
};

template<typename MatrixType_, int UpLo_, typename Ordering_>
    class SimplicialLDLT : public SimplicialCholeskyBase<SimplicialLDLT<MatrixType_,UpLo_,Ordering_> >
{
public:
    typedef MatrixType_ MatrixType;
    enum { UpLo = UpLo_ };
    typedef SimplicialCholeskyBase<SimplicialLDLT> Base;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename MatrixType::StorageIndex StorageIndex;
    typedef SparseMatrix<Scalar,ColMajor,StorageIndex> CholMatrixType;
    typedef Matrix<Scalar,Dynamic,1> VectorType;
    typedef internal::traits<SimplicialLDLT> Traits;
    typedef typename Traits::MatrixL  MatrixL;
    typedef typename Traits::MatrixU  MatrixU;
public:
    
    SimplicialLDLT() : Base() {}

    explicit SimplicialLDLT(const MatrixType& matrix)
        : Base(matrix) {}

    inline const VectorType vectorD() const {
        eigen_assert(Base::m_factorizationIsOk && "Simplicial LDLT not factorized");
        return Base::m_diag;
    }
    
    inline const MatrixL matrixL() const {
        eigen_assert(Base::m_factorizationIsOk && "Simplicial LDLT not factorized");
        return Traits::getL(Base::m_matrix);
    }

    inline const MatrixU matrixU() const {
        eigen_assert(Base::m_factorizationIsOk && "Simplicial LDLT not factorized");
        return Traits::getU(Base::m_matrix);
    }

    SimplicialLDLT& compute(const MatrixType& matrix)
    {
      Base::template compute<true>(matrix);
      return *this;
    }

    void analyzePattern(const MatrixType& a)
    {
      Base::analyzePattern(a, true);
    }

    void factorize(const MatrixType& a)
    {
      Base::template factorize<true>(a);
    }

    Scalar determinant() const
    {
      return Base::m_diag.prod();
    }
};

template<typename MatrixType_, int UpLo_, typename Ordering_>
    class SimplicialCholesky : public SimplicialCholeskyBase<SimplicialCholesky<MatrixType_,UpLo_,Ordering_> >
{
public:
    typedef MatrixType_ MatrixType;
    enum { UpLo = UpLo_ };
    typedef SimplicialCholeskyBase<SimplicialCholesky> Base;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename MatrixType::StorageIndex StorageIndex;
    typedef SparseMatrix<Scalar,ColMajor,StorageIndex> CholMatrixType;
    typedef Matrix<Scalar,Dynamic,1> VectorType;
    typedef internal::traits<SimplicialCholesky> Traits;
    typedef internal::traits<SimplicialLDLT<MatrixType,UpLo> > LDLTTraits;
    typedef internal::traits<SimplicialLLT<MatrixType,UpLo>  > LLTTraits;
  public:
    SimplicialCholesky() : Base(), m_LDLT(true) {}

    explicit SimplicialCholesky(const MatrixType& matrix)
      : Base(), m_LDLT(true)
    {
      compute(matrix);
    }

    SimplicialCholesky& setMode(SimplicialCholeskyMode mode)
    {
      switch(mode)
      {
      case SimplicialCholeskyLLT:
        m_LDLT = false;
        break;
      case SimplicialCholeskyLDLT:
        m_LDLT = true;
        break;
      default:
        break;
      }

      return *this;
    }

    inline const VectorType vectorD() const {
        eigen_assert(Base::m_factorizationIsOk && "Simplicial Cholesky not factorized");
        return Base::m_diag;
    }
    inline const CholMatrixType rawMatrix() const {
        eigen_assert(Base::m_factorizationIsOk && "Simplicial Cholesky not factorized");
        return Base::m_matrix;
    }

    SimplicialCholesky& compute(const MatrixType& matrix)
    {
      if(m_LDLT)
        Base::template compute<true>(matrix);
      else
        Base::template compute<false>(matrix);
      return *this;
    }

    void analyzePattern(const MatrixType& a)
    {
      Base::analyzePattern(a, m_LDLT);
    }

    void factorize(const MatrixType& a)
    {
      if(m_LDLT)
        Base::template factorize<true>(a);
      else
        Base::template factorize<false>(a);
    }

    template<typename Rhs,typename Dest>
    void _solve_impl(const MatrixBase<Rhs> &b, MatrixBase<Dest> &dest) const
    {
      eigen_assert(Base::m_factorizationIsOk && "The decomposition is not in a valid state for solving, you must first call either compute() or symbolic()/numeric()");
      eigen_assert(Base::m_matrix.rows()==b.rows());

      if(Base::m_info!=Success)
        return;

      if(Base::m_P.size()>0)
        dest = Base::m_P * b;
      else
        dest = b;

      if(Base::m_matrix.nonZeros()>0) 
      {
        if(m_LDLT)
          LDLTTraits::getL(Base::m_matrix).solveInPlace(dest);
        else
          LLTTraits::getL(Base::m_matrix).solveInPlace(dest);
      }

      if(Base::m_diag.size()>0)
        dest = Base::m_diag.real().asDiagonal().inverse() * dest;

      if (Base::m_matrix.nonZeros()>0) 
      {
        if(m_LDLT)
          LDLTTraits::getU(Base::m_matrix).solveInPlace(dest);
        else
          LLTTraits::getU(Base::m_matrix).solveInPlace(dest);
      }

      if(Base::m_P.size()>0)
        dest = Base::m_Pinv * dest;
    }

    template<typename Rhs,typename Dest>
    void _solve_impl(const SparseMatrixBase<Rhs> &b, SparseMatrixBase<Dest> &dest) const
    {
      internal::solve_sparse_through_dense_panels(*this, b, dest);
    }
    
    Scalar determinant() const
    {
      if(m_LDLT)
      {
        return Base::m_diag.prod();
      }
      else
      {
        Scalar detL = Diagonal<const CholMatrixType>(Base::m_matrix).prod();
        return numext::abs2(detL);
      }
    }
    
  protected:
    bool m_LDLT;
};

template<typename Derived>
void SimplicialCholeskyBase<Derived>::ordering(const MatrixType& a, ConstCholMatrixPtr &pmat, CholMatrixType& ap)
{
  eigen_assert(a.rows()==a.cols());
  const Index size = a.rows();
  pmat = &ap;
  
  if(!internal::is_same<OrderingType,NaturalOrdering<Index> >::value)
  {
    {
      CholMatrixType C;
      C = a.template selfadjointView<UpLo>();
      
      OrderingType ordering;
      ordering(C,m_Pinv);
    }

    if(m_Pinv.size()>0) m_P = m_Pinv.inverse();
    else                m_P.resize(0);
    
    ap.resize(size,size);
    ap.template selfadjointView<Upper>() = a.template selfadjointView<UpLo>().twistedBy(m_P);
  }
  else
  {
    m_Pinv.resize(0);
    m_P.resize(0);
    if(int(UpLo)==int(Lower) || MatrixType::IsRowMajor)
    {
      
      ap.resize(size,size);
      ap.template selfadjointView<Upper>() = a.template selfadjointView<UpLo>();
    }
    else
      internal::simplicial_cholesky_grab_input<CholMatrixType,MatrixType>::run(a, pmat, ap);
  }  
}

} 

#endif 
