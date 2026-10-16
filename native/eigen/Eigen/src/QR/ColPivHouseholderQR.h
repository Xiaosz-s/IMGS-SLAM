#ifndef EIGEN_COLPIVOTINGHOUSEHOLDERQR_H
#define EIGEN_COLPIVOTINGHOUSEHOLDERQR_H

namespace Eigen {

namespace internal {
template<typename MatrixType_> struct traits<ColPivHouseholderQR<MatrixType_> >
 : traits<MatrixType_>
{
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef int StorageIndex;
  enum { Flags = 0 };
};

} 

template<typename MatrixType_> class ColPivHouseholderQR
        : public SolverBase<ColPivHouseholderQR<MatrixType_> >
{
  public:

    typedef MatrixType_ MatrixType;
    typedef SolverBase<ColPivHouseholderQR> Base;
    friend class SolverBase<ColPivHouseholderQR>;

    EIGEN_GENERIC_PUBLIC_INTERFACE(ColPivHouseholderQR)
    enum {
      MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };
    typedef typename internal::plain_diag_type<MatrixType>::type HCoeffsType;
    typedef PermutationMatrix<ColsAtCompileTime, MaxColsAtCompileTime> PermutationType;
    typedef typename internal::plain_row_type<MatrixType, Index>::type IntRowVectorType;
    typedef typename internal::plain_row_type<MatrixType>::type RowVectorType;
    typedef typename internal::plain_row_type<MatrixType, RealScalar>::type RealRowVectorType;
    typedef HouseholderSequence<MatrixType,typename internal::remove_all<typename HCoeffsType::ConjugateReturnType>::type> HouseholderSequenceType;
    typedef typename MatrixType::PlainObject PlainObject;

  private:

    typedef typename PermutationType::StorageIndex PermIndexType;

  public:

    ColPivHouseholderQR()
      : m_qr(),
        m_hCoeffs(),
        m_colsPermutation(),
        m_colsTranspositions(),
        m_temp(),
        m_colNormsUpdated(),
        m_colNormsDirect(),
        m_isInitialized(false),
        m_usePrescribedThreshold(false) {}

    ColPivHouseholderQR(Index rows, Index cols)
      : m_qr(rows, cols),
        m_hCoeffs((std::min)(rows,cols)),
        m_colsPermutation(PermIndexType(cols)),
        m_colsTranspositions(cols),
        m_temp(cols),
        m_colNormsUpdated(cols),
        m_colNormsDirect(cols),
        m_isInitialized(false),
        m_usePrescribedThreshold(false) {}

    template<typename InputType>
    explicit ColPivHouseholderQR(const EigenBase<InputType>& matrix)
      : m_qr(matrix.rows(), matrix.cols()),
        m_hCoeffs((std::min)(matrix.rows(),matrix.cols())),
        m_colsPermutation(PermIndexType(matrix.cols())),
        m_colsTranspositions(matrix.cols()),
        m_temp(matrix.cols()),
        m_colNormsUpdated(matrix.cols()),
        m_colNormsDirect(matrix.cols()),
        m_isInitialized(false),
        m_usePrescribedThreshold(false)
    {
      compute(matrix.derived());
    }

    template<typename InputType>
    explicit ColPivHouseholderQR(EigenBase<InputType>& matrix)
      : m_qr(matrix.derived()),
        m_hCoeffs((std::min)(matrix.rows(),matrix.cols())),
        m_colsPermutation(PermIndexType(matrix.cols())),
        m_colsTranspositions(matrix.cols()),
        m_temp(matrix.cols()),
        m_colNormsUpdated(matrix.cols()),
        m_colNormsDirect(matrix.cols()),
        m_isInitialized(false),
        m_usePrescribedThreshold(false)
    {
      computeInPlace();
    }

    #ifdef EIGEN_PARSED_BY_DOXYGEN

    template<typename Rhs>
    inline const Solve<ColPivHouseholderQR, Rhs>
    solve(const MatrixBase<Rhs>& b) const;
    #endif

    HouseholderSequenceType householderQ() const;
    HouseholderSequenceType matrixQ() const
    {
      return householderQ();
    }

    const MatrixType& matrixQR() const
    {
      eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
      return m_qr;
    }

    const MatrixType& matrixR() const
    {
      eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
      return m_qr;
    }

    template<typename InputType>
    ColPivHouseholderQR& compute(const EigenBase<InputType>& matrix);

    const PermutationType& colsPermutation() const
    {
      eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
      return m_colsPermutation;
    }

    typename MatrixType::RealScalar absDeterminant() const;

    typename MatrixType::RealScalar logAbsDeterminant() const;

    inline Index rank() const
    {
      using std::abs;
      eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
      RealScalar premultiplied_threshold = abs(m_maxpivot) * threshold();
      Index result = 0;
      for(Index i = 0; i < m_nonzero_pivots; ++i)
        result += (abs(m_qr.coeff(i,i)) > premultiplied_threshold);
      return result;
    }

    inline Index dimensionOfKernel() const
    {
      eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
      return cols() - rank();
    }

    inline bool isInjective() const
    {
      eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
      return rank() == cols();
    }

    inline bool isSurjective() const
    {
      eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
      return rank() == rows();
    }

    inline bool isInvertible() const
    {
      eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
      return isInjective() && isSurjective();
    }

    inline const Inverse<ColPivHouseholderQR> inverse() const
    {
      eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
      return Inverse<ColPivHouseholderQR>(*this);
    }

    inline Index rows() const { return m_qr.rows(); }
    inline Index cols() const { return m_qr.cols(); }

    const HCoeffsType& hCoeffs() const { return m_hCoeffs; }

    ColPivHouseholderQR& setThreshold(const RealScalar& threshold)
    {
      m_usePrescribedThreshold = true;
      m_prescribedThreshold = threshold;
      return *this;
    }

    ColPivHouseholderQR& setThreshold(Default_t)
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
      eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
      return m_nonzero_pivots;
    }

    RealScalar maxPivot() const { return m_maxpivot; }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return Success;
    }

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    template<typename RhsType, typename DstType>
    void _solve_impl(const RhsType &rhs, DstType &dst) const;

    template<bool Conjugate, typename RhsType, typename DstType>
    void _solve_impl_transposed(const RhsType &rhs, DstType &dst) const;
    #endif

  protected:

    friend class CompleteOrthogonalDecomposition<MatrixType>;

    static void check_template_parameters()
    {
      EIGEN_STATIC_ASSERT_NON_INTEGER(Scalar);
    }

    void computeInPlace();

    MatrixType m_qr;
    HCoeffsType m_hCoeffs;
    PermutationType m_colsPermutation;
    IntRowVectorType m_colsTranspositions;
    RowVectorType m_temp;
    RealRowVectorType m_colNormsUpdated;
    RealRowVectorType m_colNormsDirect;
    bool m_isInitialized, m_usePrescribedThreshold;
    RealScalar m_prescribedThreshold, m_maxpivot;
    Index m_nonzero_pivots;
    Index m_det_pq;
};

template<typename MatrixType>
typename MatrixType::RealScalar ColPivHouseholderQR<MatrixType>::absDeterminant() const
{
  using std::abs;
  eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
  eigen_assert(m_qr.rows() == m_qr.cols() && "You can't take the determinant of a non-square matrix!");
  return abs(m_qr.diagonal().prod());
}

template<typename MatrixType>
typename MatrixType::RealScalar ColPivHouseholderQR<MatrixType>::logAbsDeterminant() const
{
  eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
  eigen_assert(m_qr.rows() == m_qr.cols() && "You can't take the determinant of a non-square matrix!");
  return m_qr.diagonal().cwiseAbs().array().log().sum();
}

template<typename MatrixType>
template<typename InputType>
ColPivHouseholderQR<MatrixType>& ColPivHouseholderQR<MatrixType>::compute(const EigenBase<InputType>& matrix)
{
  m_qr = matrix.derived();
  computeInPlace();
  return *this;
}

template<typename MatrixType>
void ColPivHouseholderQR<MatrixType>::computeInPlace()
{
  check_template_parameters();

  eigen_assert(m_qr.cols()<=NumTraits<int>::highest());

  using std::abs;

  Index rows = m_qr.rows();
  Index cols = m_qr.cols();
  Index size = m_qr.diagonalSize();

  m_hCoeffs.resize(size);

  m_temp.resize(cols);

  m_colsTranspositions.resize(m_qr.cols());
  Index number_of_transpositions = 0;

  m_colNormsUpdated.resize(cols);
  m_colNormsDirect.resize(cols);
  for (Index k = 0; k < cols; ++k) {

    m_colNormsDirect.coeffRef(k) = m_qr.col(k).norm();
    m_colNormsUpdated.coeffRef(k) = m_colNormsDirect.coeffRef(k);
  }

  RealScalar threshold_helper =  numext::abs2<RealScalar>(m_colNormsUpdated.maxCoeff() * NumTraits<RealScalar>::epsilon()) / RealScalar(rows);
  RealScalar norm_downdate_threshold = numext::sqrt(NumTraits<RealScalar>::epsilon());

  m_nonzero_pivots = size; 
  m_maxpivot = RealScalar(0);

  for(Index k = 0; k < size; ++k)
  {
    
    Index biggest_col_index;
    RealScalar biggest_col_sq_norm = numext::abs2(m_colNormsUpdated.tail(cols-k).maxCoeff(&biggest_col_index));
    biggest_col_index += k;

    if(m_nonzero_pivots==size && biggest_col_sq_norm < threshold_helper * RealScalar(rows-k))
      m_nonzero_pivots = k;

    m_colsTranspositions.coeffRef(k) = biggest_col_index;
    if(k != biggest_col_index) {
      m_qr.col(k).swap(m_qr.col(biggest_col_index));
      std::swap(m_colNormsUpdated.coeffRef(k), m_colNormsUpdated.coeffRef(biggest_col_index));
      std::swap(m_colNormsDirect.coeffRef(k), m_colNormsDirect.coeffRef(biggest_col_index));
      ++number_of_transpositions;
    }

    RealScalar beta;
    m_qr.col(k).tail(rows-k).makeHouseholderInPlace(m_hCoeffs.coeffRef(k), beta);

    m_qr.coeffRef(k,k) = beta;

    if(abs(beta) > m_maxpivot) m_maxpivot = abs(beta);

    m_qr.bottomRightCorner(rows-k, cols-k-1)
        .applyHouseholderOnTheLeft(m_qr.col(k).tail(rows-k-1), m_hCoeffs.coeffRef(k), &m_temp.coeffRef(k+1));

    for (Index j = k + 1; j < cols; ++j) {

      if (m_colNormsUpdated.coeffRef(j) != RealScalar(0)) {
        RealScalar temp = abs(m_qr.coeffRef(k, j)) / m_colNormsUpdated.coeffRef(j);
        temp = (RealScalar(1) + temp) * (RealScalar(1) - temp);
        temp = temp <  RealScalar(0) ? RealScalar(0) : temp;
        RealScalar temp2 = temp * numext::abs2<RealScalar>(m_colNormsUpdated.coeffRef(j) /
                                                           m_colNormsDirect.coeffRef(j));
        if (temp2 <= norm_downdate_threshold) {

          m_colNormsDirect.coeffRef(j) = m_qr.col(j).tail(rows - k - 1).norm();
          m_colNormsUpdated.coeffRef(j) = m_colNormsDirect.coeffRef(j);
        } else {
          m_colNormsUpdated.coeffRef(j) *= numext::sqrt(temp);
        }
      }
    }
  }

  m_colsPermutation.setIdentity(PermIndexType(cols));
  for(PermIndexType k = 0; k < size; ++k)
    m_colsPermutation.applyTranspositionOnTheRight(k, PermIndexType(m_colsTranspositions.coeff(k)));

  m_det_pq = (number_of_transpositions%2) ? -1 : 1;
  m_isInitialized = true;
}

#ifndef EIGEN_PARSED_BY_DOXYGEN
template<typename MatrixType_>
template<typename RhsType, typename DstType>
void ColPivHouseholderQR<MatrixType_>::_solve_impl(const RhsType &rhs, DstType &dst) const
{
  const Index nonzero_pivots = nonzeroPivots();

  if(nonzero_pivots == 0)
  {
    dst.setZero();
    return;
  }

  typename RhsType::PlainObject c(rhs);

  c.applyOnTheLeft(householderQ().setLength(nonzero_pivots).adjoint() );

  m_qr.topLeftCorner(nonzero_pivots, nonzero_pivots)
      .template triangularView<Upper>()
      .solveInPlace(c.topRows(nonzero_pivots));

  for(Index i = 0; i < nonzero_pivots; ++i) dst.row(m_colsPermutation.indices().coeff(i)) = c.row(i);
  for(Index i = nonzero_pivots; i < cols(); ++i) dst.row(m_colsPermutation.indices().coeff(i)).setZero();
}

template<typename MatrixType_>
template<bool Conjugate, typename RhsType, typename DstType>
void ColPivHouseholderQR<MatrixType_>::_solve_impl_transposed(const RhsType &rhs, DstType &dst) const
{
  const Index nonzero_pivots = nonzeroPivots();

  if(nonzero_pivots == 0)
  {
    dst.setZero();
    return;
  }

  typename RhsType::PlainObject c(m_colsPermutation.transpose()*rhs);

  m_qr.topLeftCorner(nonzero_pivots, nonzero_pivots)
        .template triangularView<Upper>()
        .transpose().template conjugateIf<Conjugate>()
        .solveInPlace(c.topRows(nonzero_pivots));

  dst.topRows(nonzero_pivots) = c.topRows(nonzero_pivots);
  dst.bottomRows(rows()-nonzero_pivots).setZero();

  dst.applyOnTheLeft(householderQ().setLength(nonzero_pivots).template conjugateIf<!Conjugate>() );
}
#endif

namespace internal {

template<typename DstXprType, typename MatrixType>
struct Assignment<DstXprType, Inverse<ColPivHouseholderQR<MatrixType> >, internal::assign_op<typename DstXprType::Scalar,typename ColPivHouseholderQR<MatrixType>::Scalar>, Dense2Dense>
{
  typedef ColPivHouseholderQR<MatrixType> QrType;
  typedef Inverse<QrType> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<typename DstXprType::Scalar,typename QrType::Scalar> &)
  {
    dst = src.nestedExpression().solve(MatrixType::Identity(src.rows(), src.cols()));
  }
};

} 

template<typename MatrixType>
typename ColPivHouseholderQR<MatrixType>::HouseholderSequenceType ColPivHouseholderQR<MatrixType>
  ::householderQ() const
{
  eigen_assert(m_isInitialized && "ColPivHouseholderQR is not initialized.");
  return HouseholderSequenceType(m_qr, m_hCoeffs.conjugate());
}

template<typename Derived>
const ColPivHouseholderQR<typename MatrixBase<Derived>::PlainObject>
MatrixBase<Derived>::colPivHouseholderQr() const
{
  return ColPivHouseholderQR<PlainObject>(eval());
}

} 

#endif 
