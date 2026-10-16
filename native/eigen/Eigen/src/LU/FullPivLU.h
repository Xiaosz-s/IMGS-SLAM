#ifndef EIGEN_LU_H
#define EIGEN_LU_H

namespace Eigen {

namespace internal {
template<typename MatrixType_> struct traits<FullPivLU<MatrixType_> >
 : traits<MatrixType_>
{
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef int StorageIndex;
  enum { Flags = 0 };
};

} 

template<typename MatrixType_> class FullPivLU
  : public SolverBase<FullPivLU<MatrixType_> >
{
  public:
    typedef MatrixType_ MatrixType;
    typedef SolverBase<FullPivLU> Base;
    friend class SolverBase<FullPivLU>;

    EIGEN_GENERIC_PUBLIC_INTERFACE(FullPivLU)
    enum {
      MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };
    typedef typename internal::plain_row_type<MatrixType, StorageIndex>::type IntRowVectorType;
    typedef typename internal::plain_col_type<MatrixType, StorageIndex>::type IntColVectorType;
    typedef PermutationMatrix<ColsAtCompileTime, MaxColsAtCompileTime> PermutationQType;
    typedef PermutationMatrix<RowsAtCompileTime, MaxRowsAtCompileTime> PermutationPType;
    typedef typename MatrixType::PlainObject PlainObject;

    FullPivLU();

    FullPivLU(Index rows, Index cols);

    template<typename InputType>
    explicit FullPivLU(const EigenBase<InputType>& matrix);

    template<typename InputType>
    explicit FullPivLU(EigenBase<InputType>& matrix);

    template<typename InputType>
    FullPivLU& compute(const EigenBase<InputType>& matrix) {
      m_lu = matrix.derived();
      computeInPlace();
      return *this;
    }

    inline const MatrixType& matrixLU() const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      return m_lu;
    }

    inline Index nonzeroPivots() const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      return m_nonzero_pivots;
    }

    RealScalar maxPivot() const { return m_maxpivot; }

    EIGEN_DEVICE_FUNC inline const PermutationPType& permutationP() const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      return m_p;
    }

    inline const PermutationQType& permutationQ() const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      return m_q;
    }

    inline const internal::kernel_retval<FullPivLU> kernel() const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      return internal::kernel_retval<FullPivLU>(*this);
    }

    inline const internal::image_retval<FullPivLU>
      image(const MatrixType& originalMatrix) const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      return internal::image_retval<FullPivLU>(*this, originalMatrix);
    }

    #ifdef EIGEN_PARSED_BY_DOXYGEN

    template<typename Rhs>
    inline const Solve<FullPivLU, Rhs>
    solve(const MatrixBase<Rhs>& b) const;
    #endif

    inline RealScalar rcond() const
    {
      eigen_assert(m_isInitialized && "PartialPivLU is not initialized.");
      return internal::rcond_estimate_helper(m_l1_norm, *this);
    }

    typename internal::traits<MatrixType>::Scalar determinant() const;

    FullPivLU& setThreshold(const RealScalar& threshold)
    {
      m_usePrescribedThreshold = true;
      m_prescribedThreshold = threshold;
      return *this;
    }

    FullPivLU& setThreshold(Default_t)
    {
      m_usePrescribedThreshold = false;
      return *this;
    }

    RealScalar threshold() const
    {
      eigen_assert(m_isInitialized || m_usePrescribedThreshold);
      return m_usePrescribedThreshold ? m_prescribedThreshold

          : NumTraits<Scalar>::epsilon() * RealScalar(m_lu.diagonalSize());
    }

    inline Index rank() const
    {
      using std::abs;
      eigen_assert(m_isInitialized && "LU is not initialized.");
      RealScalar premultiplied_threshold = abs(m_maxpivot) * threshold();
      Index result = 0;
      for(Index i = 0; i < m_nonzero_pivots; ++i)
        result += (abs(m_lu.coeff(i,i)) > premultiplied_threshold);
      return result;
    }

    inline Index dimensionOfKernel() const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      return cols() - rank();
    }

    inline bool isInjective() const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      return rank() == cols();
    }

    inline bool isSurjective() const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      return rank() == rows();
    }

    inline bool isInvertible() const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      return isInjective() && (m_lu.rows() == m_lu.cols());
    }

    inline const Inverse<FullPivLU> inverse() const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      eigen_assert(m_lu.rows() == m_lu.cols() && "You can't take the inverse of a non-square matrix!");
      return Inverse<FullPivLU>(*this);
    }

    MatrixType reconstructedMatrix() const;

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index rows() const EIGEN_NOEXCEPT { return m_lu.rows(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index cols() const EIGEN_NOEXCEPT { return m_lu.cols(); }

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

    void computeInPlace();

    MatrixType m_lu;
    PermutationPType m_p;
    PermutationQType m_q;
    IntColVectorType m_rowsTranspositions;
    IntRowVectorType m_colsTranspositions;
    Index m_nonzero_pivots;
    RealScalar m_l1_norm;
    RealScalar m_maxpivot, m_prescribedThreshold;
    signed char m_det_pq;
    bool m_isInitialized, m_usePrescribedThreshold;
};

template<typename MatrixType>
FullPivLU<MatrixType>::FullPivLU()
  : m_isInitialized(false), m_usePrescribedThreshold(false)
{
}

template<typename MatrixType>
FullPivLU<MatrixType>::FullPivLU(Index rows, Index cols)
  : m_lu(rows, cols),
    m_p(rows),
    m_q(cols),
    m_rowsTranspositions(rows),
    m_colsTranspositions(cols),
    m_isInitialized(false),
    m_usePrescribedThreshold(false)
{
}

template<typename MatrixType>
template<typename InputType>
FullPivLU<MatrixType>::FullPivLU(const EigenBase<InputType>& matrix)
  : m_lu(matrix.rows(), matrix.cols()),
    m_p(matrix.rows()),
    m_q(matrix.cols()),
    m_rowsTranspositions(matrix.rows()),
    m_colsTranspositions(matrix.cols()),
    m_isInitialized(false),
    m_usePrescribedThreshold(false)
{
  compute(matrix.derived());
}

template<typename MatrixType>
template<typename InputType>
FullPivLU<MatrixType>::FullPivLU(EigenBase<InputType>& matrix)
  : m_lu(matrix.derived()),
    m_p(matrix.rows()),
    m_q(matrix.cols()),
    m_rowsTranspositions(matrix.rows()),
    m_colsTranspositions(matrix.cols()),
    m_isInitialized(false),
    m_usePrescribedThreshold(false)
{
  computeInPlace();
}

template<typename MatrixType>
void FullPivLU<MatrixType>::computeInPlace()
{
  check_template_parameters();

  eigen_assert(m_lu.rows()<=NumTraits<int>::highest() && m_lu.cols()<=NumTraits<int>::highest());

  m_l1_norm = m_lu.cwiseAbs().colwise().sum().maxCoeff();

  const Index size = m_lu.diagonalSize();
  const Index rows = m_lu.rows();
  const Index cols = m_lu.cols();

  m_rowsTranspositions.resize(m_lu.rows());
  m_colsTranspositions.resize(m_lu.cols());
  Index number_of_transpositions = 0; 

  m_nonzero_pivots = size; 
  m_maxpivot = RealScalar(0);

  for(Index k = 0; k < size; ++k)
  {

    Index row_of_biggest_in_corner, col_of_biggest_in_corner;
    typedef internal::scalar_score_coeff_op<Scalar> Scoring;
    typedef typename Scoring::result_type Score;
    Score biggest_in_corner;
    biggest_in_corner = m_lu.bottomRightCorner(rows-k, cols-k)
                        .unaryExpr(Scoring())
                        .maxCoeff(&row_of_biggest_in_corner, &col_of_biggest_in_corner);
    row_of_biggest_in_corner += k; 
    col_of_biggest_in_corner += k; 

    if(biggest_in_corner==Score(0))
    {

      m_nonzero_pivots = k;
      for(Index i = k; i < size; ++i)
      {
        m_rowsTranspositions.coeffRef(i) = internal::convert_index<StorageIndex>(i);
        m_colsTranspositions.coeffRef(i) = internal::convert_index<StorageIndex>(i);
      }
      break;
    }

    RealScalar abs_pivot = internal::abs_knowing_score<Scalar>()(m_lu(row_of_biggest_in_corner, col_of_biggest_in_corner), biggest_in_corner);
    if(abs_pivot > m_maxpivot) m_maxpivot = abs_pivot;

    m_rowsTranspositions.coeffRef(k) = internal::convert_index<StorageIndex>(row_of_biggest_in_corner);
    m_colsTranspositions.coeffRef(k) = internal::convert_index<StorageIndex>(col_of_biggest_in_corner);
    if(k != row_of_biggest_in_corner) {
      m_lu.row(k).swap(m_lu.row(row_of_biggest_in_corner));
      ++number_of_transpositions;
    }
    if(k != col_of_biggest_in_corner) {
      m_lu.col(k).swap(m_lu.col(col_of_biggest_in_corner));
      ++number_of_transpositions;
    }

    if(k<rows-1)
      m_lu.col(k).tail(rows-k-1) /= m_lu.coeff(k,k);
    if(k<size-1)
      m_lu.block(k+1,k+1,rows-k-1,cols-k-1).noalias() -= m_lu.col(k).tail(rows-k-1) * m_lu.row(k).tail(cols-k-1);
  }

  m_p.setIdentity(rows);
  for(Index k = size-1; k >= 0; --k)
    m_p.applyTranspositionOnTheRight(k, m_rowsTranspositions.coeff(k));

  m_q.setIdentity(cols);
  for(Index k = 0; k < size; ++k)
    m_q.applyTranspositionOnTheRight(k, m_colsTranspositions.coeff(k));

  m_det_pq = (number_of_transpositions%2) ? -1 : 1;

  m_isInitialized = true;
}

template<typename MatrixType>
typename internal::traits<MatrixType>::Scalar FullPivLU<MatrixType>::determinant() const
{
  eigen_assert(m_isInitialized && "LU is not initialized.");
  eigen_assert(m_lu.rows() == m_lu.cols() && "You can't take the determinant of a non-square matrix!");
  return Scalar(m_det_pq) * Scalar(m_lu.diagonal().prod());
}

template<typename MatrixType>
MatrixType FullPivLU<MatrixType>::reconstructedMatrix() const
{
  eigen_assert(m_isInitialized && "LU is not initialized.");
  const Index smalldim = (std::min)(m_lu.rows(), m_lu.cols());
  
  MatrixType res(m_lu.rows(),m_lu.cols());
  
  res = m_lu.leftCols(smalldim)
            .template triangularView<UnitLower>().toDenseMatrix()
      * m_lu.topRows(smalldim)
            .template triangularView<Upper>().toDenseMatrix();

  res = m_p.inverse() * res;

  res = res * m_q.inverse();

  return res;
}

namespace internal {
template<typename MatrixType_>
struct kernel_retval<FullPivLU<MatrixType_> >
  : kernel_retval_base<FullPivLU<MatrixType_> >
{
  EIGEN_MAKE_KERNEL_HELPERS(FullPivLU<MatrixType_>)

  enum { MaxSmallDimAtCompileTime = EIGEN_SIZE_MIN_PREFER_FIXED(
            MatrixType::MaxColsAtCompileTime,
            MatrixType::MaxRowsAtCompileTime)
  };

  template<typename Dest> void evalTo(Dest& dst) const
  {
    using std::abs;
    const Index cols = dec().matrixLU().cols(), dimker = cols - rank();
    if(dimker == 0)
    {

      dst.setZero();
      return;
    }

    Matrix<Index, Dynamic, 1, 0, MaxSmallDimAtCompileTime, 1> pivots(rank());
    RealScalar premultiplied_threshold = dec().maxPivot() * dec().threshold();
    Index p = 0;
    for(Index i = 0; i < dec().nonzeroPivots(); ++i)
      if(abs(dec().matrixLU().coeff(i,i)) > premultiplied_threshold)
        pivots.coeffRef(p++) = i;
    eigen_internal_assert(p == rank());

    Matrix<typename MatrixType::Scalar, Dynamic, Dynamic, MatrixType::Options,
           MaxSmallDimAtCompileTime, MatrixType::MaxColsAtCompileTime>
      m(dec().matrixLU().block(0, 0, rank(), cols));
    for(Index i = 0; i < rank(); ++i)
    {
      if(i) m.row(i).head(i).setZero();
      m.row(i).tail(cols-i) = dec().matrixLU().row(pivots.coeff(i)).tail(cols-i);
    }
    m.block(0, 0, rank(), rank());
    m.block(0, 0, rank(), rank()).template triangularView<StrictlyLower>().setZero();
    for(Index i = 0; i < rank(); ++i)
      m.col(i).swap(m.col(pivots.coeff(i)));

    m.topLeftCorner(rank(), rank())
     .template triangularView<Upper>().solveInPlace(
        m.topRightCorner(rank(), dimker)
      );

    for(Index i = rank()-1; i >= 0; --i)
      m.col(i).swap(m.col(pivots.coeff(i)));

    for(Index i = 0; i < rank(); ++i) dst.row(dec().permutationQ().indices().coeff(i)) = -m.row(i).tail(dimker);
    for(Index i = rank(); i < cols; ++i) dst.row(dec().permutationQ().indices().coeff(i)).setZero();
    for(Index k = 0; k < dimker; ++k) dst.coeffRef(dec().permutationQ().indices().coeff(rank()+k), k) = Scalar(1);
  }
};

template<typename MatrixType_>
struct image_retval<FullPivLU<MatrixType_> >
  : image_retval_base<FullPivLU<MatrixType_> >
{
  EIGEN_MAKE_IMAGE_HELPERS(FullPivLU<MatrixType_>)

  enum { MaxSmallDimAtCompileTime = EIGEN_SIZE_MIN_PREFER_FIXED(
            MatrixType::MaxColsAtCompileTime,
            MatrixType::MaxRowsAtCompileTime)
  };

  template<typename Dest> void evalTo(Dest& dst) const
  {
    using std::abs;
    if(rank() == 0)
    {

      dst.setZero();
      return;
    }

    Matrix<Index, Dynamic, 1, 0, MaxSmallDimAtCompileTime, 1> pivots(rank());
    RealScalar premultiplied_threshold = dec().maxPivot() * dec().threshold();
    Index p = 0;
    for(Index i = 0; i < dec().nonzeroPivots(); ++i)
      if(abs(dec().matrixLU().coeff(i,i)) > premultiplied_threshold)
        pivots.coeffRef(p++) = i;
    eigen_internal_assert(p == rank());

    for(Index i = 0; i < rank(); ++i)
      dst.col(i) = originalMatrix().col(dec().permutationQ().indices().coeff(pivots.coeff(i)));
  }
};

} 

#ifndef EIGEN_PARSED_BY_DOXYGEN
template<typename MatrixType_>
template<typename RhsType, typename DstType>
void FullPivLU<MatrixType_>::_solve_impl(const RhsType &rhs, DstType &dst) const
{

  const Index rows = this->rows(),
              cols = this->cols(),
              nonzero_pivots = this->rank();
  const Index smalldim = (std::min)(rows, cols);

  if(nonzero_pivots == 0)
  {
    dst.setZero();
    return;
  }

  typename RhsType::PlainObject c(rhs.rows(), rhs.cols());

  c = permutationP() * rhs;

  m_lu.topLeftCorner(smalldim,smalldim)
      .template triangularView<UnitLower>()
      .solveInPlace(c.topRows(smalldim));
  if(rows>cols)
    c.bottomRows(rows-cols) -= m_lu.bottomRows(rows-cols) * c.topRows(cols);

  m_lu.topLeftCorner(nonzero_pivots, nonzero_pivots)
      .template triangularView<Upper>()
      .solveInPlace(c.topRows(nonzero_pivots));

  for(Index i = 0; i < nonzero_pivots; ++i)
    dst.row(permutationQ().indices().coeff(i)) = c.row(i);
  for(Index i = nonzero_pivots; i < m_lu.cols(); ++i)
    dst.row(permutationQ().indices().coeff(i)).setZero();
}

template<typename MatrixType_>
template<bool Conjugate, typename RhsType, typename DstType>
void FullPivLU<MatrixType_>::_solve_impl_transposed(const RhsType &rhs, DstType &dst) const
{

  const Index rows = this->rows(), cols = this->cols(),
    nonzero_pivots = this->rank();
  const Index smalldim = (std::min)(rows, cols);

  if(nonzero_pivots == 0)
  {
    dst.setZero();
    return;
  }

  typename RhsType::PlainObject c(rhs.rows(), rhs.cols());

  c = permutationQ().inverse() * rhs;

  m_lu.topLeftCorner(nonzero_pivots, nonzero_pivots)
      .template triangularView<Upper>()
      .transpose()
      .template conjugateIf<Conjugate>()
      .solveInPlace(c.topRows(nonzero_pivots));

  m_lu.topLeftCorner(smalldim, smalldim)
      .template triangularView<UnitLower>()
      .transpose()
      .template conjugateIf<Conjugate>()
      .solveInPlace(c.topRows(smalldim));

  PermutationPType invp = permutationP().inverse().eval();
  for(Index i = 0; i < smalldim; ++i)
    dst.row(invp.indices().coeff(i)) = c.row(i);
  for(Index i = smalldim; i < rows; ++i)
    dst.row(invp.indices().coeff(i)).setZero();
}

#endif

namespace internal {

template<typename DstXprType, typename MatrixType>
struct Assignment<DstXprType, Inverse<FullPivLU<MatrixType> >, internal::assign_op<typename DstXprType::Scalar,typename FullPivLU<MatrixType>::Scalar>, Dense2Dense>
{
  typedef FullPivLU<MatrixType> LuType;
  typedef Inverse<LuType> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<typename DstXprType::Scalar,typename MatrixType::Scalar> &)
  {
    dst = src.nestedExpression().solve(MatrixType::Identity(src.rows(), src.cols()));
  }
};
} 

template<typename Derived>
inline const FullPivLU<typename MatrixBase<Derived>::PlainObject>
MatrixBase<Derived>::fullPivLu() const
{
  return FullPivLU<PlainObject>(eval());
}

} 

#endif 
