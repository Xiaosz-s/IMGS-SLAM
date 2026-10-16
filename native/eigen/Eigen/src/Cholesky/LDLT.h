#ifndef EIGEN_LDLT_H
#define EIGEN_LDLT_H

namespace Eigen {

namespace internal {
  template<typename MatrixType_, int UpLo_> struct traits<LDLT<MatrixType_, UpLo_> >
   : traits<MatrixType_>
  {
    typedef MatrixXpr XprKind;
    typedef SolverStorage StorageKind;
    typedef int StorageIndex;
    enum { Flags = 0 };
  };

  template<typename MatrixType, int UpLo> struct LDLT_Traits;

  enum SignMatrix { PositiveSemiDef, NegativeSemiDef, ZeroSign, Indefinite };
}

template<typename MatrixType_, int UpLo_> class LDLT
        : public SolverBase<LDLT<MatrixType_, UpLo_> >
{
  public:
    typedef MatrixType_ MatrixType;
    typedef SolverBase<LDLT> Base;
    friend class SolverBase<LDLT>;

    EIGEN_GENERIC_PUBLIC_INTERFACE(LDLT)
    enum {
      MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime,
      UpLo = UpLo_
    };
    typedef Matrix<Scalar, RowsAtCompileTime, 1, 0, MaxRowsAtCompileTime, 1> TmpMatrixType;

    typedef Transpositions<RowsAtCompileTime, MaxRowsAtCompileTime> TranspositionType;
    typedef PermutationMatrix<RowsAtCompileTime, MaxRowsAtCompileTime> PermutationType;

    typedef internal::LDLT_Traits<MatrixType,UpLo> Traits;

    LDLT()
      : m_matrix(),
        m_transpositions(),
        m_sign(internal::ZeroSign),
        m_isInitialized(false)
    {}

    explicit LDLT(Index size)
      : m_matrix(size, size),
        m_transpositions(size),
        m_temporary(size),
        m_sign(internal::ZeroSign),
        m_isInitialized(false)
    {}

    template<typename InputType>
    explicit LDLT(const EigenBase<InputType>& matrix)
      : m_matrix(matrix.rows(), matrix.cols()),
        m_transpositions(matrix.rows()),
        m_temporary(matrix.rows()),
        m_sign(internal::ZeroSign),
        m_isInitialized(false)
    {
      compute(matrix.derived());
    }

    template<typename InputType>
    explicit LDLT(EigenBase<InputType>& matrix)
      : m_matrix(matrix.derived()),
        m_transpositions(matrix.rows()),
        m_temporary(matrix.rows()),
        m_sign(internal::ZeroSign),
        m_isInitialized(false)
    {
      compute(matrix.derived());
    }

    void setZero()
    {
      m_isInitialized = false;
    }

    inline typename Traits::MatrixU matrixU() const
    {
      eigen_assert(m_isInitialized && "LDLT is not initialized.");
      return Traits::getU(m_matrix);
    }

    inline typename Traits::MatrixL matrixL() const
    {
      eigen_assert(m_isInitialized && "LDLT is not initialized.");
      return Traits::getL(m_matrix);
    }

    inline const TranspositionType& transpositionsP() const
    {
      eigen_assert(m_isInitialized && "LDLT is not initialized.");
      return m_transpositions;
    }

    inline Diagonal<const MatrixType> vectorD() const
    {
      eigen_assert(m_isInitialized && "LDLT is not initialized.");
      return m_matrix.diagonal();
    }

    inline bool isPositive() const
    {
      eigen_assert(m_isInitialized && "LDLT is not initialized.");
      return m_sign == internal::PositiveSemiDef || m_sign == internal::ZeroSign;
    }

    inline bool isNegative(void) const
    {
      eigen_assert(m_isInitialized && "LDLT is not initialized.");
      return m_sign == internal::NegativeSemiDef || m_sign == internal::ZeroSign;
    }

    #ifdef EIGEN_PARSED_BY_DOXYGEN

    template<typename Rhs>
    inline const Solve<LDLT, Rhs>
    solve(const MatrixBase<Rhs>& b) const;
    #endif

    template<typename Derived>
    bool solveInPlace(MatrixBase<Derived> &bAndX) const;

    template<typename InputType>
    LDLT& compute(const EigenBase<InputType>& matrix);

    RealScalar rcond() const
    {
      eigen_assert(m_isInitialized && "LDLT is not initialized.");
      return internal::rcond_estimate_helper(m_l1_norm, *this);
    }

    template <typename Derived>
    LDLT& rankUpdate(const MatrixBase<Derived>& w, const RealScalar& alpha=1);

    inline const MatrixType& matrixLDLT() const
    {
      eigen_assert(m_isInitialized && "LDLT is not initialized.");
      return m_matrix;
    }

    MatrixType reconstructedMatrix() const;

    const LDLT& adjoint() const { return *this; };

    EIGEN_DEVICE_FUNC inline EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return m_matrix.rows(); }
    EIGEN_DEVICE_FUNC inline EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return m_matrix.cols(); }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "LDLT is not initialized.");
      return m_info;
    }

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    template<typename RhsType, typename DstType>
    void _solve_impl(const RhsType &rhs, DstType &dst) const;

    template<bool Conjugate, typename RhsType, typename DstType>
    void _solve_impl_transposed(const RhsType &rhs, DstType &dst) const;
    #endif

  protected:

    static void check_template_parameters()
    {
      EIGEN_STATIC_ASSERT_NON_INTEGER(Scalar);
    }

    MatrixType m_matrix;
    RealScalar m_l1_norm;
    TranspositionType m_transpositions;
    TmpMatrixType m_temporary;
    internal::SignMatrix m_sign;
    bool m_isInitialized;
    ComputationInfo m_info;
};

namespace internal {

template<int UpLo> struct ldlt_inplace;

template<> struct ldlt_inplace<Lower>
{
  template<typename MatrixType, typename TranspositionType, typename Workspace>
  static bool unblocked(MatrixType& mat, TranspositionType& transpositions, Workspace& temp, SignMatrix& sign)
  {
    using std::abs;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename TranspositionType::StorageIndex IndexType;
    eigen_assert(mat.rows()==mat.cols());
    const Index size = mat.rows();
    bool found_zero_pivot = false;
    bool ret = true;

    if (size <= 1)
    {
      transpositions.setIdentity();
      if(size==0) sign = ZeroSign;
      else if (numext::real(mat.coeff(0,0)) > static_cast<RealScalar>(0) ) sign = PositiveSemiDef;
      else if (numext::real(mat.coeff(0,0)) < static_cast<RealScalar>(0)) sign = NegativeSemiDef;
      else sign = ZeroSign;
      return true;
    }

    for (Index k = 0; k < size; ++k)
    {
      
      Index index_of_biggest_in_corner;
      mat.diagonal().tail(size-k).cwiseAbs().maxCoeff(&index_of_biggest_in_corner);
      index_of_biggest_in_corner += k;

      transpositions.coeffRef(k) = IndexType(index_of_biggest_in_corner);
      if(k != index_of_biggest_in_corner)
      {

        Index s = size-index_of_biggest_in_corner-1; 
        mat.row(k).head(k).swap(mat.row(index_of_biggest_in_corner).head(k));
        mat.col(k).tail(s).swap(mat.col(index_of_biggest_in_corner).tail(s));
        std::swap(mat.coeffRef(k,k),mat.coeffRef(index_of_biggest_in_corner,index_of_biggest_in_corner));
        for(Index i=k+1;i<index_of_biggest_in_corner;++i)
        {
          Scalar tmp = mat.coeffRef(i,k);
          mat.coeffRef(i,k) = numext::conj(mat.coeffRef(index_of_biggest_in_corner,i));
          mat.coeffRef(index_of_biggest_in_corner,i) = numext::conj(tmp);
        }
        if(NumTraits<Scalar>::IsComplex)
          mat.coeffRef(index_of_biggest_in_corner,k) = numext::conj(mat.coeff(index_of_biggest_in_corner,k));
      }

      Index rs = size - k - 1;
      Block<MatrixType,Dynamic,1> A21(mat,k+1,k,rs,1);
      Block<MatrixType,1,Dynamic> A10(mat,k,0,1,k);
      Block<MatrixType,Dynamic,Dynamic> A20(mat,k+1,0,rs,k);

      if(k>0)
      {
        temp.head(k) = mat.diagonal().real().head(k).asDiagonal() * A10.adjoint();
        mat.coeffRef(k,k) -= (A10 * temp.head(k)).value();
        if(rs>0)
          A21.noalias() -= A20 * temp.head(k);
      }

      RealScalar realAkk = numext::real(mat.coeffRef(k,k));
      bool pivot_is_valid = (abs(realAkk) > RealScalar(0));

      if(k==0 && !pivot_is_valid)
      {

        sign = ZeroSign;
        for(Index j = 0; j<size; ++j)
        {
          transpositions.coeffRef(j) = IndexType(j);
          ret = ret && (mat.col(j).tail(size-j-1).array()==Scalar(0)).all();
        }
        return ret;
      }

      if((rs>0) && pivot_is_valid)
        A21 /= realAkk;
      else if(rs>0)
        ret = ret && (A21.array()==Scalar(0)).all();

      if(found_zero_pivot && pivot_is_valid) ret = false; 
      else if(!pivot_is_valid) found_zero_pivot = true;

      if (sign == PositiveSemiDef) {
        if (realAkk < static_cast<RealScalar>(0)) sign = Indefinite;
      } else if (sign == NegativeSemiDef) {
        if (realAkk > static_cast<RealScalar>(0)) sign = Indefinite;
      } else if (sign == ZeroSign) {
        if (realAkk > static_cast<RealScalar>(0)) sign = PositiveSemiDef;
        else if (realAkk < static_cast<RealScalar>(0)) sign = NegativeSemiDef;
      }
    }

    return ret;
  }

  template<typename MatrixType, typename WDerived>
  static bool updateInPlace(MatrixType& mat, MatrixBase<WDerived>& w, const typename MatrixType::RealScalar& sigma=1)
  {
    using numext::isfinite;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;

    const Index size = mat.rows();
    eigen_assert(mat.cols() == size && w.size()==size);

    RealScalar alpha = 1;

    for (Index j = 0; j < size; j++)
    {
      
      if (!(isfinite)(alpha))
        break;

      RealScalar dj = numext::real(mat.coeff(j,j));
      Scalar wj = w.coeff(j);
      RealScalar swj2 = sigma*numext::abs2(wj);
      RealScalar gamma = dj*alpha + swj2;

      mat.coeffRef(j,j) += swj2/alpha;
      alpha += swj2/dj;

      Index rs = size-j-1;
      w.tail(rs) -= wj * mat.col(j).tail(rs);
      if(gamma != 0)
        mat.col(j).tail(rs) += (sigma*numext::conj(wj)/gamma)*w.tail(rs);
    }
    return true;
  }

  template<typename MatrixType, typename TranspositionType, typename Workspace, typename WType>
  static bool update(MatrixType& mat, const TranspositionType& transpositions, Workspace& tmp, const WType& w, const typename MatrixType::RealScalar& sigma=1)
  {
    
    tmp = transpositions * w;

    return ldlt_inplace<Lower>::updateInPlace(mat,tmp,sigma);
  }
};

template<> struct ldlt_inplace<Upper>
{
  template<typename MatrixType, typename TranspositionType, typename Workspace>
  static EIGEN_STRONG_INLINE bool unblocked(MatrixType& mat, TranspositionType& transpositions, Workspace& temp, SignMatrix& sign)
  {
    Transpose<MatrixType> matt(mat);
    return ldlt_inplace<Lower>::unblocked(matt, transpositions, temp, sign);
  }

  template<typename MatrixType, typename TranspositionType, typename Workspace, typename WType>
  static EIGEN_STRONG_INLINE bool update(MatrixType& mat, TranspositionType& transpositions, Workspace& tmp, WType& w, const typename MatrixType::RealScalar& sigma=1)
  {
    Transpose<MatrixType> matt(mat);
    return ldlt_inplace<Lower>::update(matt, transpositions, tmp, w.conjugate(), sigma);
  }
};

template<typename MatrixType> struct LDLT_Traits<MatrixType,Lower>
{
  typedef const TriangularView<const MatrixType, UnitLower> MatrixL;
  typedef const TriangularView<const typename MatrixType::AdjointReturnType, UnitUpper> MatrixU;
  static inline MatrixL getL(const MatrixType& m) { return MatrixL(m); }
  static inline MatrixU getU(const MatrixType& m) { return MatrixU(m.adjoint()); }
};

template<typename MatrixType> struct LDLT_Traits<MatrixType,Upper>
{
  typedef const TriangularView<const typename MatrixType::AdjointReturnType, UnitLower> MatrixL;
  typedef const TriangularView<const MatrixType, UnitUpper> MatrixU;
  static inline MatrixL getL(const MatrixType& m) { return MatrixL(m.adjoint()); }
  static inline MatrixU getU(const MatrixType& m) { return MatrixU(m); }
};

} 

template<typename MatrixType, int UpLo_>
template<typename InputType>
LDLT<MatrixType,UpLo_>& LDLT<MatrixType,UpLo_>::compute(const EigenBase<InputType>& a)
{
  check_template_parameters();

  eigen_assert(a.rows()==a.cols());
  const Index size = a.rows();

  m_matrix = a.derived();

  m_l1_norm = RealScalar(0);
  
  for (Index col = 0; col < size; ++col) {
    RealScalar abs_col_sum;
    if (UpLo_ == Lower)
      abs_col_sum = m_matrix.col(col).tail(size - col).template lpNorm<1>() + m_matrix.row(col).head(col).template lpNorm<1>();
    else
      abs_col_sum = m_matrix.col(col).head(col).template lpNorm<1>() + m_matrix.row(col).tail(size - col).template lpNorm<1>();
    if (abs_col_sum > m_l1_norm)
      m_l1_norm = abs_col_sum;
  }

  m_transpositions.resize(size);
  m_isInitialized = false;
  m_temporary.resize(size);
  m_sign = internal::ZeroSign;

  m_info = internal::ldlt_inplace<UpLo>::unblocked(m_matrix, m_transpositions, m_temporary, m_sign) ? Success : NumericalIssue;

  m_isInitialized = true;
  return *this;
}

template<typename MatrixType, int UpLo_>
template<typename Derived>
LDLT<MatrixType,UpLo_>& LDLT<MatrixType,UpLo_>::rankUpdate(const MatrixBase<Derived>& w, const typename LDLT<MatrixType,UpLo_>::RealScalar& sigma)
{
  typedef typename TranspositionType::StorageIndex IndexType;
  const Index size = w.rows();
  if (m_isInitialized)
  {
    eigen_assert(m_matrix.rows()==size);
  }
  else
  {
    m_matrix.resize(size,size);
    m_matrix.setZero();
    m_transpositions.resize(size);
    for (Index i = 0; i < size; i++)
      m_transpositions.coeffRef(i) = IndexType(i);
    m_temporary.resize(size);
    m_sign = sigma>=0 ? internal::PositiveSemiDef : internal::NegativeSemiDef;
    m_isInitialized = true;
  }

  internal::ldlt_inplace<UpLo>::update(m_matrix, m_transpositions, m_temporary, w, sigma);

  return *this;
}

#ifndef EIGEN_PARSED_BY_DOXYGEN
template<typename MatrixType_, int UpLo_>
template<typename RhsType, typename DstType>
void LDLT<MatrixType_,UpLo_>::_solve_impl(const RhsType &rhs, DstType &dst) const
{
  _solve_impl_transposed<true>(rhs, dst);
}

template<typename MatrixType_,int UpLo_>
template<bool Conjugate, typename RhsType, typename DstType>
void LDLT<MatrixType_,UpLo_>::_solve_impl_transposed(const RhsType &rhs, DstType &dst) const
{
  
  dst = m_transpositions * rhs;

  matrixL().template conjugateIf<!Conjugate>().solveInPlace(dst);

  using std::abs;
  const typename Diagonal<const MatrixType>::RealReturnType vecD(vectorD());

  RealScalar tolerance = (std::numeric_limits<RealScalar>::min)();
  for (Index i = 0; i < vecD.size(); ++i)
  {
    if(abs(vecD(i)) > tolerance)
      dst.row(i) /= vecD(i);
    else
      dst.row(i).setZero();
  }

  matrixL().transpose().template conjugateIf<Conjugate>().solveInPlace(dst);

  dst = m_transpositions.transpose() * dst;
}
#endif

template<typename MatrixType,int UpLo_>
template<typename Derived>
bool LDLT<MatrixType,UpLo_>::solveInPlace(MatrixBase<Derived> &bAndX) const
{
  eigen_assert(m_isInitialized && "LDLT is not initialized.");
  eigen_assert(m_matrix.rows() == bAndX.rows());

  bAndX = this->solve(bAndX);

  return true;
}

template<typename MatrixType, int UpLo_>
MatrixType LDLT<MatrixType,UpLo_>::reconstructedMatrix() const
{
  eigen_assert(m_isInitialized && "LDLT is not initialized.");
  const Index size = m_matrix.rows();
  MatrixType res(size,size);

  res.setIdentity();
  res = transpositionsP() * res;
  
  res = matrixU() * res;
  
  res = vectorD().real().asDiagonal() * res;
  
  res = matrixL() * res;
  
  res = transpositionsP().transpose() * res;

  return res;
}

template<typename MatrixType, unsigned int UpLo>
inline const LDLT<typename SelfAdjointView<MatrixType, UpLo>::PlainObject, UpLo>
SelfAdjointView<MatrixType, UpLo>::ldlt() const
{
  return LDLT<PlainObject,UpLo>(m_matrix);
}

template<typename Derived>
inline const LDLT<typename MatrixBase<Derived>::PlainObject>
MatrixBase<Derived>::ldlt() const
{
  return LDLT<PlainObject>(derived());
}

} 

#endif 
