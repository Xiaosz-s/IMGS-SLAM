#ifndef EIGEN_FULLPIVOTINGHOUSEHOLDERQR_H
#define EIGEN_FULLPIVOTINGHOUSEHOLDERQR_H

namespace Eigen { 

namespace internal {

template<typename MatrixType_> struct traits<FullPivHouseholderQR<MatrixType_> >
 : traits<MatrixType_>
{
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef int StorageIndex;
  enum { Flags = 0 };
};

template<typename MatrixType> struct FullPivHouseholderQRMatrixQReturnType;

template<typename MatrixType>
struct traits<FullPivHouseholderQRMatrixQReturnType<MatrixType> >
{
  typedef typename MatrixType::PlainObject ReturnType;
};

} 

template<typename MatrixType_> class FullPivHouseholderQR
        : public SolverBase<FullPivHouseholderQR<MatrixType_> >
{
  public:

    typedef MatrixType_ MatrixType;
    typedef SolverBase<FullPivHouseholderQR> Base;
    friend class SolverBase<FullPivHouseholderQR>;

    EIGEN_GENERIC_PUBLIC_INTERFACE(FullPivHouseholderQR)
    enum {
      MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };
    typedef internal::FullPivHouseholderQRMatrixQReturnType<MatrixType> MatrixQReturnType;
    typedef typename internal::plain_diag_type<MatrixType>::type HCoeffsType;
    typedef Matrix<StorageIndex, 1,
                   EIGEN_SIZE_MIN_PREFER_DYNAMIC(ColsAtCompileTime,RowsAtCompileTime), RowMajor, 1,
                   EIGEN_SIZE_MIN_PREFER_FIXED(MaxColsAtCompileTime,MaxRowsAtCompileTime)> IntDiagSizeVectorType;
    typedef PermutationMatrix<ColsAtCompileTime, MaxColsAtCompileTime> PermutationType;
    typedef typename internal::plain_row_type<MatrixType>::type RowVectorType;
    typedef typename internal::plain_col_type<MatrixType>::type ColVectorType;
    typedef typename MatrixType::PlainObject PlainObject;

    FullPivHouseholderQR()
      : m_qr(),
        m_hCoeffs(),
        m_rows_transpositions(),
        m_cols_transpositions(),
        m_cols_permutation(),
        m_temp(),
        m_isInitialized(false),
        m_usePrescribedThreshold(false) {}

    FullPivHouseholderQR(Index rows, Index cols)
      : m_qr(rows, cols),
        m_hCoeffs((std::min)(rows,cols)),
        m_rows_transpositions((std::min)(rows,cols)),
        m_cols_transpositions((std::min)(rows,cols)),
        m_cols_permutation(cols),
        m_temp(cols),
        m_isInitialized(false),
        m_usePrescribedThreshold(false) {}

    template<typename InputType>
    explicit FullPivHouseholderQR(const EigenBase<InputType>& matrix)
      : m_qr(matrix.rows(), matrix.cols()),
        m_hCoeffs((std::min)(matrix.rows(), matrix.cols())),
        m_rows_transpositions((std::min)(matrix.rows(), matrix.cols())),
        m_cols_transpositions((std::min)(matrix.rows(), matrix.cols())),
        m_cols_permutation(matrix.cols()),
        m_temp(matrix.cols()),
        m_isInitialized(false),
        m_usePrescribedThreshold(false)
    {
      compute(matrix.derived());
    }

    template<typename InputType>
    explicit FullPivHouseholderQR(EigenBase<InputType>& matrix)
      : m_qr(matrix.derived()),
        m_hCoeffs((std::min)(matrix.rows(), matrix.cols())),
        m_rows_transpositions((std::min)(matrix.rows(), matrix.cols())),
        m_cols_transpositions((std::min)(matrix.rows(), matrix.cols())),
        m_cols_permutation(matrix.cols()),
        m_temp(matrix.cols()),
        m_isInitialized(false),
        m_usePrescribedThreshold(false)
    {
      computeInPlace();
    }

    #ifdef EIGEN_PARSED_BY_DOXYGEN

    template<typename Rhs>
    inline const Solve<FullPivHouseholderQR, Rhs>
    solve(const MatrixBase<Rhs>& b) const;
    #endif

    MatrixQReturnType matrixQ(void) const;

    const MatrixType& matrixQR() const
    {
      eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
      return m_qr;
    }

    template<typename InputType>
    FullPivHouseholderQR& compute(const EigenBase<InputType>& matrix);

    const PermutationType& colsPermutation() const
    {
      eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
      return m_cols_permutation;
    }

    const IntDiagSizeVectorType& rowsTranspositions() const
    {
      eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
      return m_rows_transpositions;
    }

    typename MatrixType::RealScalar absDeterminant() const;

    typename MatrixType::RealScalar logAbsDeterminant() const;

    inline Index rank() const
    {
      using std::abs;
      eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
      RealScalar premultiplied_threshold = abs(m_maxpivot) * threshold();
      Index result = 0;
      for(Index i = 0; i < m_nonzero_pivots; ++i)
        result += (abs(m_qr.coeff(i,i)) > premultiplied_threshold);
      return result;
    }

    inline Index dimensionOfKernel() const
    {
      eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
      return cols() - rank();
    }

    inline bool isInjective() const
    {
      eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
      return rank() == cols();
    }

    inline bool isSurjective() const
    {
      eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
      return rank() == rows();
    }

    inline bool isInvertible() const
    {
      eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
      return isInjective() && isSurjective();
    }

    inline const Inverse<FullPivHouseholderQR> inverse() const
    {
      eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
      return Inverse<FullPivHouseholderQR>(*this);
    }

    inline Index rows() const { return m_qr.rows(); }
    inline Index cols() const { return m_qr.cols(); }

    const HCoeffsType& hCoeffs() const { return m_hCoeffs; }

    FullPivHouseholderQR& setThreshold(const RealScalar& threshold)
    {
      m_usePrescribedThreshold = true;
      m_prescribedThreshold = threshold;
      return *this;
    }

    FullPivHouseholderQR& setThreshold(Default_t)
    {
      m_usePrescribedThreshold = false;
      return *this;
    }

    RealScalar threshold() const
    {
      eigen_assert(m_isInitialized || m_usePrescribedThreshold);
      return m_usePrescribedThreshold ? m_prescribedThreshold

                                      : NumTraits<Scalar>::epsilon() * RealScalar(m_qr.diagonalSize());
    }

    inline Index nonzeroPivots() const
    {
      eigen_assert(m_isInitialized && "LU is not initialized.");
      return m_nonzero_pivots;
    }

    RealScalar maxPivot() const { return m_maxpivot; }

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

    MatrixType m_qr;
    HCoeffsType m_hCoeffs;
    IntDiagSizeVectorType m_rows_transpositions;
    IntDiagSizeVectorType m_cols_transpositions;
    PermutationType m_cols_permutation;
    RowVectorType m_temp;
    bool m_isInitialized, m_usePrescribedThreshold;
    RealScalar m_prescribedThreshold, m_maxpivot;
    Index m_nonzero_pivots;
    RealScalar m_precision;
    Index m_det_pq;
};

template<typename MatrixType>
typename MatrixType::RealScalar FullPivHouseholderQR<MatrixType>::absDeterminant() const
{
  using std::abs;
  eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
  eigen_assert(m_qr.rows() == m_qr.cols() && "You can't take the determinant of a non-square matrix!");
  return abs(m_qr.diagonal().prod());
}

template<typename MatrixType>
typename MatrixType::RealScalar FullPivHouseholderQR<MatrixType>::logAbsDeterminant() const
{
  eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
  eigen_assert(m_qr.rows() == m_qr.cols() && "You can't take the determinant of a non-square matrix!");
  return m_qr.diagonal().cwiseAbs().array().log().sum();
}

template<typename MatrixType>
template<typename InputType>
FullPivHouseholderQR<MatrixType>& FullPivHouseholderQR<MatrixType>::compute(const EigenBase<InputType>& matrix)
{
  m_qr = matrix.derived();
  computeInPlace();
  return *this;
}

template<typename MatrixType>
void FullPivHouseholderQR<MatrixType>::computeInPlace()
{
  check_template_parameters();

  using std::abs;
  Index rows = m_qr.rows();
  Index cols = m_qr.cols();
  Index size = (std::min)(rows,cols);

  m_hCoeffs.resize(size);

  m_temp.resize(cols);

  m_precision = NumTraits<Scalar>::epsilon() * RealScalar(size);

  m_rows_transpositions.resize(size);
  m_cols_transpositions.resize(size);
  Index number_of_transpositions = 0;

  RealScalar biggest(0);

  m_nonzero_pivots = size; 
  m_maxpivot = RealScalar(0);

  for (Index k = 0; k < size; ++k)
  {
    Index row_of_biggest_in_corner, col_of_biggest_in_corner;
    typedef internal::scalar_score_coeff_op<Scalar> Scoring;
    typedef typename Scoring::result_type Score;

    Score score = m_qr.bottomRightCorner(rows-k, cols-k)
                      .unaryExpr(Scoring())
                      .maxCoeff(&row_of_biggest_in_corner, &col_of_biggest_in_corner);
    row_of_biggest_in_corner += k;
    col_of_biggest_in_corner += k;
    RealScalar biggest_in_corner = internal::abs_knowing_score<Scalar>()(m_qr(row_of_biggest_in_corner, col_of_biggest_in_corner), score);
    if(k==0) biggest = biggest_in_corner;

    if(internal::isMuchSmallerThan(biggest_in_corner, biggest, m_precision))
    {
      m_nonzero_pivots = k;
      for(Index i = k; i < size; i++)
      {
        m_rows_transpositions.coeffRef(i) = internal::convert_index<StorageIndex>(i);
        m_cols_transpositions.coeffRef(i) = internal::convert_index<StorageIndex>(i);
        m_hCoeffs.coeffRef(i) = Scalar(0);
      }
      break;
    }

    m_rows_transpositions.coeffRef(k) = internal::convert_index<StorageIndex>(row_of_biggest_in_corner);
    m_cols_transpositions.coeffRef(k) = internal::convert_index<StorageIndex>(col_of_biggest_in_corner);
    if(k != row_of_biggest_in_corner) {
      m_qr.row(k).tail(cols-k).swap(m_qr.row(row_of_biggest_in_corner).tail(cols-k));
      ++number_of_transpositions;
    }
    if(k != col_of_biggest_in_corner) {
      m_qr.col(k).swap(m_qr.col(col_of_biggest_in_corner));
      ++number_of_transpositions;
    }

    RealScalar beta;
    m_qr.col(k).tail(rows-k).makeHouseholderInPlace(m_hCoeffs.coeffRef(k), beta);
    m_qr.coeffRef(k,k) = beta;

    if(abs(beta) > m_maxpivot) m_maxpivot = abs(beta);

    m_qr.bottomRightCorner(rows-k, cols-k-1)
        .applyHouseholderOnTheLeft(m_qr.col(k).tail(rows-k-1), m_hCoeffs.coeffRef(k), &m_temp.coeffRef(k+1));
  }

  m_cols_permutation.setIdentity(cols);
  for(Index k = 0; k < size; ++k)
    m_cols_permutation.applyTranspositionOnTheRight(k, m_cols_transpositions.coeff(k));

  m_det_pq = (number_of_transpositions%2) ? -1 : 1;
  m_isInitialized = true;
}

#ifndef EIGEN_PARSED_BY_DOXYGEN
template<typename MatrixType_>
template<typename RhsType, typename DstType>
void FullPivHouseholderQR<MatrixType_>::_solve_impl(const RhsType &rhs, DstType &dst) const
{
  const Index l_rank = rank();

  if(l_rank==0)
  {
    dst.setZero();
    return;
  }

  typename RhsType::PlainObject c(rhs);

  Matrix<typename RhsType::Scalar,1,RhsType::ColsAtCompileTime> temp(rhs.cols());
  for (Index k = 0; k < l_rank; ++k)
  {
    Index remainingSize = rows()-k;
    c.row(k).swap(c.row(m_rows_transpositions.coeff(k)));
    c.bottomRightCorner(remainingSize, rhs.cols())
      .applyHouseholderOnTheLeft(m_qr.col(k).tail(remainingSize-1),
                               m_hCoeffs.coeff(k), &temp.coeffRef(0));
  }

  m_qr.topLeftCorner(l_rank, l_rank)
      .template triangularView<Upper>()
      .solveInPlace(c.topRows(l_rank));

  for(Index i = 0; i < l_rank; ++i) dst.row(m_cols_permutation.indices().coeff(i)) = c.row(i);
  for(Index i = l_rank; i < cols(); ++i) dst.row(m_cols_permutation.indices().coeff(i)).setZero();
}

template<typename MatrixType_>
template<bool Conjugate, typename RhsType, typename DstType>
void FullPivHouseholderQR<MatrixType_>::_solve_impl_transposed(const RhsType &rhs, DstType &dst) const
{
  const Index l_rank = rank();

  if(l_rank == 0)
  {
    dst.setZero();
    return;
  }

  typename RhsType::PlainObject c(m_cols_permutation.transpose()*rhs);

  m_qr.topLeftCorner(l_rank, l_rank)
         .template triangularView<Upper>()
         .transpose().template conjugateIf<Conjugate>()
         .solveInPlace(c.topRows(l_rank));

  dst.topRows(l_rank) = c.topRows(l_rank);
  dst.bottomRows(rows()-l_rank).setZero();

  Matrix<Scalar, 1, DstType::ColsAtCompileTime> temp(dst.cols());
  const Index size = (std::min)(rows(), cols());
  for (Index k = size-1; k >= 0; --k)
  {
    Index remainingSize = rows()-k;

    dst.bottomRightCorner(remainingSize, dst.cols())
       .applyHouseholderOnTheLeft(m_qr.col(k).tail(remainingSize-1).template conjugateIf<!Conjugate>(),
                                  m_hCoeffs.template conjugateIf<Conjugate>().coeff(k), &temp.coeffRef(0));

    dst.row(k).swap(dst.row(m_rows_transpositions.coeff(k)));
  }
}
#endif

namespace internal {
  
template<typename DstXprType, typename MatrixType>
struct Assignment<DstXprType, Inverse<FullPivHouseholderQR<MatrixType> >, internal::assign_op<typename DstXprType::Scalar,typename FullPivHouseholderQR<MatrixType>::Scalar>, Dense2Dense>
{
  typedef FullPivHouseholderQR<MatrixType> QrType;
  typedef Inverse<QrType> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<typename DstXprType::Scalar,typename QrType::Scalar> &)
  {    
    dst = src.nestedExpression().solve(MatrixType::Identity(src.rows(), src.cols()));
  }
};

template<typename MatrixType> struct FullPivHouseholderQRMatrixQReturnType
  : public ReturnByValue<FullPivHouseholderQRMatrixQReturnType<MatrixType> >
{
public:
  typedef typename FullPivHouseholderQR<MatrixType>::IntDiagSizeVectorType IntDiagSizeVectorType;
  typedef typename internal::plain_diag_type<MatrixType>::type HCoeffsType;
  typedef Matrix<typename MatrixType::Scalar, 1, MatrixType::RowsAtCompileTime, RowMajor, 1,
                 MatrixType::MaxRowsAtCompileTime> WorkVectorType;

  FullPivHouseholderQRMatrixQReturnType(const MatrixType&       qr,
                                        const HCoeffsType&      hCoeffs,
                                        const IntDiagSizeVectorType& rowsTranspositions)
    : m_qr(qr),
      m_hCoeffs(hCoeffs),
      m_rowsTranspositions(rowsTranspositions)
  {}

  template <typename ResultType>
  void evalTo(ResultType& result) const
  {
    const Index rows = m_qr.rows();
    WorkVectorType workspace(rows);
    evalTo(result, workspace);
  }

  template <typename ResultType>
  void evalTo(ResultType& result, WorkVectorType& workspace) const
  {
    using numext::conj;

    const Index rows = m_qr.rows();
    const Index cols = m_qr.cols();
    const Index size = (std::min)(rows, cols);
    workspace.resize(rows);
    result.setIdentity(rows, rows);
    for (Index k = size-1; k >= 0; k--)
    {
      result.block(k, k, rows-k, rows-k)
            .applyHouseholderOnTheLeft(m_qr.col(k).tail(rows-k-1), conj(m_hCoeffs.coeff(k)), &workspace.coeffRef(k));
      result.row(k).swap(result.row(m_rowsTranspositions.coeff(k)));
    }
  }

  Index rows() const { return m_qr.rows(); }
  Index cols() const { return m_qr.rows(); }

protected:
  typename MatrixType::Nested m_qr;
  typename HCoeffsType::Nested m_hCoeffs;
  typename IntDiagSizeVectorType::Nested m_rowsTranspositions;
};

} 

template<typename MatrixType>
inline typename FullPivHouseholderQR<MatrixType>::MatrixQReturnType FullPivHouseholderQR<MatrixType>::matrixQ() const
{
  eigen_assert(m_isInitialized && "FullPivHouseholderQR is not initialized.");
  return MatrixQReturnType(m_qr, m_hCoeffs, m_rows_transpositions);
}

template<typename Derived>
const FullPivHouseholderQR<typename MatrixBase<Derived>::PlainObject>
MatrixBase<Derived>::fullPivHouseholderQr() const
{
  return FullPivHouseholderQR<PlainObject>(eval());
}

} 

#endif 
