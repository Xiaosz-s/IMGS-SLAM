#ifndef EIGEN_KLUSUPPORT_H
#define EIGEN_KLUSUPPORT_H

namespace Eigen {

inline int klu_solve(klu_symbolic *Symbolic, klu_numeric *Numeric, Index ldim, Index nrhs, double B [ ], klu_common *Common, double) {
   return klu_solve(Symbolic, Numeric, internal::convert_index<int>(ldim), internal::convert_index<int>(nrhs), B, Common);
}

inline int klu_solve(klu_symbolic *Symbolic, klu_numeric *Numeric, Index ldim, Index nrhs, std::complex<double>B[], klu_common *Common, std::complex<double>) {
   return klu_z_solve(Symbolic, Numeric, internal::convert_index<int>(ldim), internal::convert_index<int>(nrhs), &numext::real_ref(B[0]), Common);
}

inline int klu_tsolve(klu_symbolic *Symbolic, klu_numeric *Numeric, Index ldim, Index nrhs, double B[], klu_common *Common, double) {
   return klu_tsolve(Symbolic, Numeric, internal::convert_index<int>(ldim), internal::convert_index<int>(nrhs), B, Common);
}

inline int klu_tsolve(klu_symbolic *Symbolic, klu_numeric *Numeric, Index ldim, Index nrhs, std::complex<double>B[], klu_common *Common, std::complex<double>) {
   return klu_z_tsolve(Symbolic, Numeric, internal::convert_index<int>(ldim), internal::convert_index<int>(nrhs), &numext::real_ref(B[0]), 0, Common);
}

inline klu_numeric* klu_factor(int Ap [ ], int Ai [ ], double Ax [ ], klu_symbolic *Symbolic, klu_common *Common, double) {
   return klu_factor(Ap, Ai, Ax, Symbolic, Common);
}

inline klu_numeric* klu_factor(int Ap[], int Ai[], std::complex<double> Ax[], klu_symbolic *Symbolic, klu_common *Common, std::complex<double>) {
   return klu_z_factor(Ap, Ai, &numext::real_ref(Ax[0]), Symbolic, Common);
}

template<typename MatrixType_>
class KLU : public SparseSolverBase<KLU<MatrixType_> >
{
  protected:
    typedef SparseSolverBase<KLU<MatrixType_> > Base;
    using Base::m_isInitialized;
  public:
    using Base::_solve_impl;
    typedef MatrixType_ MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename MatrixType::StorageIndex StorageIndex;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    typedef Matrix<int, 1, MatrixType::ColsAtCompileTime> IntRowVectorType;
    typedef Matrix<int, MatrixType::RowsAtCompileTime, 1> IntColVectorType;
    typedef SparseMatrix<Scalar> LUMatrixType;
    typedef SparseMatrix<Scalar,ColMajor,int> KLUMatrixType;
    typedef Ref<const KLUMatrixType, StandardCompressedFormat> KLUMatrixRef;
    enum {
      ColsAtCompileTime = MatrixType::ColsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };

  public:

    KLU()
      : m_dummy(0,0), mp_matrix(m_dummy)
    {
      init();
    }

    template<typename InputMatrixType>
    explicit KLU(const InputMatrixType& matrix)
      : mp_matrix(matrix)
    {
      init();
      compute(matrix);
    }

    ~KLU()
    {
      if(m_symbolic) klu_free_symbolic(&m_symbolic,&m_common);
      if(m_numeric)  klu_free_numeric(&m_numeric,&m_common);
    }

    EIGEN_CONSTEXPR inline Index rows() const EIGEN_NOEXCEPT { return mp_matrix.rows(); }
    EIGEN_CONSTEXPR inline Index cols() const EIGEN_NOEXCEPT { return mp_matrix.cols(); }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_info;
    }
#if 0 
    inline const LUMatrixType& matrixL() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_l;
    }

    inline const LUMatrixType& matrixU() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_u;
    }

    inline const IntColVectorType& permutationP() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_p;
    }

    inline const IntRowVectorType& permutationQ() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_q;
    }
#endif

    template<typename InputMatrixType>
    void compute(const InputMatrixType& matrix)
    {
      if(m_symbolic) klu_free_symbolic(&m_symbolic, &m_common);
      if(m_numeric)  klu_free_numeric(&m_numeric, &m_common);
      grab(matrix.derived());
      analyzePattern_impl();
      factorize_impl();
    }

    template<typename InputMatrixType>
    void analyzePattern(const InputMatrixType& matrix)
    {
      if(m_symbolic) klu_free_symbolic(&m_symbolic, &m_common);
      if(m_numeric)  klu_free_numeric(&m_numeric, &m_common);

      grab(matrix.derived());

      analyzePattern_impl();
    }

    inline const klu_common& kluCommon() const
    {
      return m_common;
    }

    inline klu_common& kluCommon()
    {
      return m_common;
    }

    template<typename InputMatrixType>
    void factorize(const InputMatrixType& matrix)
    {
      eigen_assert(m_analysisIsOk && "KLU: you must first call analyzePattern()");
      if(m_numeric)
        klu_free_numeric(&m_numeric,&m_common);

      grab(matrix.derived());

      factorize_impl();
    }

    template<typename BDerived,typename XDerived>
    bool _solve_impl(const MatrixBase<BDerived> &b, MatrixBase<XDerived> &x) const;

#if 0 
    Scalar determinant() const;

    void extractData() const;
#endif

  protected:

    void init()
    {
      m_info                  = InvalidInput;
      m_isInitialized         = false;
      m_numeric               = 0;
      m_symbolic              = 0;
      m_extractedDataAreDirty = true;

      klu_defaults(&m_common);
    }

    void analyzePattern_impl()
    {
      m_info = InvalidInput;
      m_analysisIsOk = false;
      m_factorizationIsOk = false;
      m_symbolic = klu_analyze(internal::convert_index<int>(mp_matrix.rows()),
                                     const_cast<StorageIndex*>(mp_matrix.outerIndexPtr()), const_cast<StorageIndex*>(mp_matrix.innerIndexPtr()),
                                     &m_common);
      if (m_symbolic) {
         m_isInitialized = true;
         m_info = Success;
         m_analysisIsOk = true;
         m_extractedDataAreDirty = true;
      }
    }

    void factorize_impl()
    {

      m_numeric = klu_factor(const_cast<StorageIndex*>(mp_matrix.outerIndexPtr()), const_cast<StorageIndex*>(mp_matrix.innerIndexPtr()), const_cast<Scalar*>(mp_matrix.valuePtr()),
                                    m_symbolic, &m_common, Scalar());

      m_info = m_numeric ? Success : NumericalIssue;
      m_factorizationIsOk = m_numeric ? 1 : 0;
      m_extractedDataAreDirty = true;
    }

    template<typename MatrixDerived>
    void grab(const EigenBase<MatrixDerived> &A)
    {
      mp_matrix.~KLUMatrixRef();
      ::new (&mp_matrix) KLUMatrixRef(A.derived());
    }

    void grab(const KLUMatrixRef &A)
    {
      if(&(A.derived()) != &mp_matrix)
      {
        mp_matrix.~KLUMatrixRef();
        ::new (&mp_matrix) KLUMatrixRef(A);
      }
    }

#if 0 
    mutable LUMatrixType m_l;
    mutable LUMatrixType m_u;
    mutable IntColVectorType m_p;
    mutable IntRowVectorType m_q;
#endif

    KLUMatrixType m_dummy;
    KLUMatrixRef mp_matrix;

    klu_numeric* m_numeric;
    klu_symbolic* m_symbolic;
    klu_common m_common;
    mutable ComputationInfo m_info;
    int m_factorizationIsOk;
    int m_analysisIsOk;
    mutable bool m_extractedDataAreDirty;

  private:
    KLU(const KLU& ) { }
};

#if 0 
template<typename MatrixType>
void KLU<MatrixType>::extractData() const
{
  if (m_extractedDataAreDirty)
  {
     eigen_assert(false && "KLU: extractData Not Yet Implemented");

    int lnz, unz, rows, cols, nz_udiag;
    umfpack_get_lunz(&lnz, &unz, &rows, &cols, &nz_udiag, m_numeric, Scalar());

    m_l.resize(rows,(std::min)(rows,cols));
    m_l.resizeNonZeros(lnz);

    m_u.resize((std::min)(rows,cols),cols);
    m_u.resizeNonZeros(unz);

    m_p.resize(rows);
    m_q.resize(cols);

    umfpack_get_numeric(m_l.outerIndexPtr(), m_l.innerIndexPtr(), m_l.valuePtr(),
                        m_u.outerIndexPtr(), m_u.innerIndexPtr(), m_u.valuePtr(),
                        m_p.data(), m_q.data(), 0, 0, 0, m_numeric);

    m_extractedDataAreDirty = false;
  }
}

template<typename MatrixType>
typename KLU<MatrixType>::Scalar KLU<MatrixType>::determinant() const
{
  eigen_assert(false && "KLU: extractData Not Yet Implemented");
  return Scalar();
}
#endif

template<typename MatrixType>
template<typename BDerived,typename XDerived>
bool KLU<MatrixType>::_solve_impl(const MatrixBase<BDerived> &b, MatrixBase<XDerived> &x) const
{
  Index rhsCols = b.cols();
  EIGEN_STATIC_ASSERT((XDerived::Flags&RowMajorBit)==0, THIS_METHOD_IS_ONLY_FOR_COLUMN_MAJOR_MATRICES);
  eigen_assert(m_factorizationIsOk && "The decomposition is not in a valid state for solving, you must first call either compute() or analyzePattern()/factorize()");

  x = b;
  int info = klu_solve(m_symbolic, m_numeric, b.rows(), rhsCols, x.const_cast_derived().data(), const_cast<klu_common*>(&m_common), Scalar());

  m_info = info!=0 ? Success : NumericalIssue;
  return true;
}

} 

#endif 
