#ifndef EIGEN_LLT_H
#define EIGEN_LLT_H

namespace Eigen {

namespace internal{

template<typename MatrixType_, int UpLo_> struct traits<LLT<MatrixType_, UpLo_> >
 : traits<MatrixType_>
{
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef int StorageIndex;
  enum { Flags = 0 };
};

template<typename MatrixType, int UpLo> struct LLT_Traits;
}

template<typename MatrixType_, int UpLo_> class LLT
        : public SolverBase<LLT<MatrixType_, UpLo_> >
{
  public:
    typedef MatrixType_ MatrixType;
    typedef SolverBase<LLT> Base;
    friend class SolverBase<LLT>;

    EIGEN_GENERIC_PUBLIC_INTERFACE(LLT)
    enum {
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };

    enum {
      PacketSize = internal::packet_traits<Scalar>::size,
      AlignmentMask = int(PacketSize)-1,
      UpLo = UpLo_
    };

    typedef internal::LLT_Traits<MatrixType,UpLo> Traits;

    LLT() : m_matrix(), m_isInitialized(false) {}

    explicit LLT(Index size) : m_matrix(size, size),
                    m_isInitialized(false) {}

    template<typename InputType>
    explicit LLT(const EigenBase<InputType>& matrix)
      : m_matrix(matrix.rows(), matrix.cols()),
        m_isInitialized(false)
    {
      compute(matrix.derived());
    }

    template<typename InputType>
    explicit LLT(EigenBase<InputType>& matrix)
      : m_matrix(matrix.derived()),
        m_isInitialized(false)
    {
      compute(matrix.derived());
    }

    inline typename Traits::MatrixU matrixU() const
    {
      eigen_assert(m_isInitialized && "LLT is not initialized.");
      return Traits::getU(m_matrix);
    }

    inline typename Traits::MatrixL matrixL() const
    {
      eigen_assert(m_isInitialized && "LLT is not initialized.");
      return Traits::getL(m_matrix);
    }

    #ifdef EIGEN_PARSED_BY_DOXYGEN

    template<typename Rhs>
    inline const Solve<LLT, Rhs>
    solve(const MatrixBase<Rhs>& b) const;
    #endif

    template<typename Derived>
    void solveInPlace(const MatrixBase<Derived> &bAndX) const;

    template<typename InputType>
    LLT& compute(const EigenBase<InputType>& matrix);

    RealScalar rcond() const
    {
      eigen_assert(m_isInitialized && "LLT is not initialized.");
      eigen_assert(m_info == Success && "LLT failed because matrix appears to be negative");
      return internal::rcond_estimate_helper(m_l1_norm, *this);
    }

    inline const MatrixType& matrixLLT() const
    {
      eigen_assert(m_isInitialized && "LLT is not initialized.");
      return m_matrix;
    }

    MatrixType reconstructedMatrix() const;

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "LLT is not initialized.");
      return m_info;
    }

    const LLT& adjoint() const EIGEN_NOEXCEPT { return *this; };

    inline EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return m_matrix.rows(); }
    inline EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return m_matrix.cols(); }

    template<typename VectorType>
    LLT & rankUpdate(const VectorType& vec, const RealScalar& sigma = 1);

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
    bool m_isInitialized;
    ComputationInfo m_info;
};

namespace internal {

template<typename Scalar, int UpLo> struct llt_inplace;

template<typename MatrixType, typename VectorType>
static Index llt_rank_update_lower(MatrixType& mat, const VectorType& vec, const typename MatrixType::RealScalar& sigma)
{
  using std::sqrt;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;
  typedef typename MatrixType::ColXpr ColXpr;
  typedef typename internal::remove_all<ColXpr>::type ColXprCleaned;
  typedef typename ColXprCleaned::SegmentReturnType ColXprSegment;
  typedef Matrix<Scalar,Dynamic,1> TempVectorType;
  typedef typename TempVectorType::SegmentReturnType TempVecSegment;

  Index n = mat.cols();
  eigen_assert(mat.rows()==n && vec.size()==n);

  TempVectorType temp;

  if(sigma>0)
  {

    temp = sqrt(sigma) * vec;

    for(Index i=0; i<n; ++i)
    {
      JacobiRotation<Scalar> g;
      g.makeGivens(mat(i,i), -temp(i), &mat(i,i));

      Index rs = n-i-1;
      if(rs>0)
      {
        ColXprSegment x(mat.col(i).tail(rs));
        TempVecSegment y(temp.tail(rs));
        apply_rotation_in_the_plane(x, y, g);
      }
    }
  }
  else
  {
    temp = vec;
    RealScalar beta = 1;
    for(Index j=0; j<n; ++j)
    {
      RealScalar Ljj = numext::real(mat.coeff(j,j));
      RealScalar dj = numext::abs2(Ljj);
      Scalar wj = temp.coeff(j);
      RealScalar swj2 = sigma*numext::abs2(wj);
      RealScalar gamma = dj*beta + swj2;

      RealScalar x = dj + swj2/beta;
      if (x<=RealScalar(0))
        return j;
      RealScalar nLjj = sqrt(x);
      mat.coeffRef(j,j) = nLjj;
      beta += swj2/dj;

      Index rs = n-j-1;
      if(rs)
      {
        temp.tail(rs) -= (wj/Ljj) * mat.col(j).tail(rs);
        if(gamma != 0)
          mat.col(j).tail(rs) = (nLjj/Ljj) * mat.col(j).tail(rs) + (nLjj * sigma*numext::conj(wj)/gamma)*temp.tail(rs);
      }
    }
  }
  return -1;
}

template<typename Scalar> struct llt_inplace<Scalar, Lower>
{
  typedef typename NumTraits<Scalar>::Real RealScalar;
  template<typename MatrixType>
  static Index unblocked(MatrixType& mat)
  {
    using std::sqrt;

    eigen_assert(mat.rows()==mat.cols());
    const Index size = mat.rows();
    for(Index k = 0; k < size; ++k)
    {
      Index rs = size-k-1; 

      Block<MatrixType,Dynamic,1> A21(mat,k+1,k,rs,1);
      Block<MatrixType,1,Dynamic> A10(mat,k,0,1,k);
      Block<MatrixType,Dynamic,Dynamic> A20(mat,k+1,0,rs,k);

      RealScalar x = numext::real(mat.coeff(k,k));
      if (k>0) x -= A10.squaredNorm();
      if (x<=RealScalar(0))
        return k;
      mat.coeffRef(k,k) = x = sqrt(x);
      if (k>0 && rs>0) A21.noalias() -= A20 * A10.adjoint();
      if (rs>0) A21 /= x;
    }
    return -1;
  }

  template<typename MatrixType>
  static Index blocked(MatrixType& m)
  {
    eigen_assert(m.rows()==m.cols());
    Index size = m.rows();
    if(size<32)
      return unblocked(m);

    Index blockSize = size/8;
    blockSize = (blockSize/16)*16;
    blockSize = (std::min)((std::max)(blockSize,Index(8)), Index(128));

    for (Index k=0; k<size; k+=blockSize)
    {

      Index bs = (std::min)(blockSize, size-k);
      Index rs = size - k - bs;
      Block<MatrixType,Dynamic,Dynamic> A11(m,k,   k,   bs,bs);
      Block<MatrixType,Dynamic,Dynamic> A21(m,k+bs,k,   rs,bs);
      Block<MatrixType,Dynamic,Dynamic> A22(m,k+bs,k+bs,rs,rs);

      Index ret;
      if((ret=unblocked(A11))>=0) return k+ret;
      if(rs>0) A11.adjoint().template triangularView<Upper>().template solveInPlace<OnTheRight>(A21);
      if(rs>0) A22.template selfadjointView<Lower>().rankUpdate(A21,typename NumTraits<RealScalar>::Literal(-1)); 
    }
    return -1;
  }

  template<typename MatrixType, typename VectorType>
  static Index rankUpdate(MatrixType& mat, const VectorType& vec, const RealScalar& sigma)
  {
    return Eigen::internal::llt_rank_update_lower(mat, vec, sigma);
  }
};

template<typename Scalar> struct llt_inplace<Scalar, Upper>
{
  typedef typename NumTraits<Scalar>::Real RealScalar;

  template<typename MatrixType>
  static EIGEN_STRONG_INLINE Index unblocked(MatrixType& mat)
  {
    Transpose<MatrixType> matt(mat);
    return llt_inplace<Scalar, Lower>::unblocked(matt);
  }
  template<typename MatrixType>
  static EIGEN_STRONG_INLINE Index blocked(MatrixType& mat)
  {
    Transpose<MatrixType> matt(mat);
    return llt_inplace<Scalar, Lower>::blocked(matt);
  }
  template<typename MatrixType, typename VectorType>
  static Index rankUpdate(MatrixType& mat, const VectorType& vec, const RealScalar& sigma)
  {
    Transpose<MatrixType> matt(mat);
    return llt_inplace<Scalar, Lower>::rankUpdate(matt, vec.conjugate(), sigma);
  }
};

template<typename MatrixType> struct LLT_Traits<MatrixType,Lower>
{
  typedef const TriangularView<const MatrixType, Lower> MatrixL;
  typedef const TriangularView<const typename MatrixType::AdjointReturnType, Upper> MatrixU;
  static inline MatrixL getL(const MatrixType& m) { return MatrixL(m); }
  static inline MatrixU getU(const MatrixType& m) { return MatrixU(m.adjoint()); }
  static bool inplace_decomposition(MatrixType& m)
  { return llt_inplace<typename MatrixType::Scalar, Lower>::blocked(m)==-1; }
};

template<typename MatrixType> struct LLT_Traits<MatrixType,Upper>
{
  typedef const TriangularView<const typename MatrixType::AdjointReturnType, Lower> MatrixL;
  typedef const TriangularView<const MatrixType, Upper> MatrixU;
  static inline MatrixL getL(const MatrixType& m) { return MatrixL(m.adjoint()); }
  static inline MatrixU getU(const MatrixType& m) { return MatrixU(m); }
  static bool inplace_decomposition(MatrixType& m)
  { return llt_inplace<typename MatrixType::Scalar, Upper>::blocked(m)==-1; }
};

} 

template<typename MatrixType, int UpLo_>
template<typename InputType>
LLT<MatrixType,UpLo_>& LLT<MatrixType,UpLo_>::compute(const EigenBase<InputType>& a)
{
  check_template_parameters();

  eigen_assert(a.rows()==a.cols());
  const Index size = a.rows();
  m_matrix.resize(size, size);
  if (!internal::is_same_dense(m_matrix, a.derived()))
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

  m_isInitialized = true;
  bool ok = Traits::inplace_decomposition(m_matrix);
  m_info = ok ? Success : NumericalIssue;

  return *this;
}

template<typename MatrixType_, int UpLo_>
template<typename VectorType>
LLT<MatrixType_,UpLo_> & LLT<MatrixType_,UpLo_>::rankUpdate(const VectorType& v, const RealScalar& sigma)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorType);
  eigen_assert(v.size()==m_matrix.cols());
  eigen_assert(m_isInitialized);
  if(internal::llt_inplace<typename MatrixType::Scalar, UpLo>::rankUpdate(m_matrix,v,sigma)>=0)
    m_info = NumericalIssue;
  else
    m_info = Success;

  return *this;
}

#ifndef EIGEN_PARSED_BY_DOXYGEN
template<typename MatrixType_,int UpLo_>
template<typename RhsType, typename DstType>
void LLT<MatrixType_,UpLo_>::_solve_impl(const RhsType &rhs, DstType &dst) const
{
  _solve_impl_transposed<true>(rhs, dst);
}

template<typename MatrixType_,int UpLo_>
template<bool Conjugate, typename RhsType, typename DstType>
void LLT<MatrixType_,UpLo_>::_solve_impl_transposed(const RhsType &rhs, DstType &dst) const
{
    dst = rhs;

    matrixL().template conjugateIf<!Conjugate>().solveInPlace(dst);
    matrixU().template conjugateIf<!Conjugate>().solveInPlace(dst);
}
#endif

template<typename MatrixType, int UpLo_>
template<typename Derived>
void LLT<MatrixType,UpLo_>::solveInPlace(const MatrixBase<Derived> &bAndX) const
{
  eigen_assert(m_isInitialized && "LLT is not initialized.");
  eigen_assert(m_matrix.rows()==bAndX.rows());
  matrixL().solveInPlace(bAndX);
  matrixU().solveInPlace(bAndX);
}

template<typename MatrixType, int UpLo_>
MatrixType LLT<MatrixType,UpLo_>::reconstructedMatrix() const
{
  eigen_assert(m_isInitialized && "LLT is not initialized.");
  return matrixL() * matrixL().adjoint().toDenseMatrix();
}

template<typename Derived>
inline const LLT<typename MatrixBase<Derived>::PlainObject>
MatrixBase<Derived>::llt() const
{
  return LLT<PlainObject>(derived());
}

template<typename MatrixType, unsigned int UpLo>
inline const LLT<typename SelfAdjointView<MatrixType, UpLo>::PlainObject, UpLo>
SelfAdjointView<MatrixType, UpLo>::llt() const
{
  return LLT<PlainObject,UpLo>(m_matrix);
}

} 

#endif 
