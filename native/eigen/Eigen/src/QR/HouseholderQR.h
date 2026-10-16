#ifndef EIGEN_QR_H
#define EIGEN_QR_H

namespace Eigen { 

namespace internal {
template<typename MatrixType_> struct traits<HouseholderQR<MatrixType_> >
 : traits<MatrixType_>
{
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef int StorageIndex;
  enum { Flags = 0 };
};

} 

template<typename MatrixType_> class HouseholderQR
        : public SolverBase<HouseholderQR<MatrixType_> >
{
  public:

    typedef MatrixType_ MatrixType;
    typedef SolverBase<HouseholderQR> Base;
    friend class SolverBase<HouseholderQR>;

    EIGEN_GENERIC_PUBLIC_INTERFACE(HouseholderQR)
    enum {
      MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };
    typedef Matrix<Scalar, RowsAtCompileTime, RowsAtCompileTime, (MatrixType::Flags&RowMajorBit) ? RowMajor : ColMajor, MaxRowsAtCompileTime, MaxRowsAtCompileTime> MatrixQType;
    typedef typename internal::plain_diag_type<MatrixType>::type HCoeffsType;
    typedef typename internal::plain_row_type<MatrixType>::type RowVectorType;
    typedef HouseholderSequence<MatrixType,typename internal::remove_all<typename HCoeffsType::ConjugateReturnType>::type> HouseholderSequenceType;

    HouseholderQR() : m_qr(), m_hCoeffs(), m_temp(), m_isInitialized(false) {}

    HouseholderQR(Index rows, Index cols)
      : m_qr(rows, cols),
        m_hCoeffs((std::min)(rows,cols)),
        m_temp(cols),
        m_isInitialized(false) {}

    template<typename InputType>
    explicit HouseholderQR(const EigenBase<InputType>& matrix)
      : m_qr(matrix.rows(), matrix.cols()),
        m_hCoeffs((std::min)(matrix.rows(),matrix.cols())),
        m_temp(matrix.cols()),
        m_isInitialized(false)
    {
      compute(matrix.derived());
    }

    template<typename InputType>
    explicit HouseholderQR(EigenBase<InputType>& matrix)
      : m_qr(matrix.derived()),
        m_hCoeffs((std::min)(matrix.rows(),matrix.cols())),
        m_temp(matrix.cols()),
        m_isInitialized(false)
    {
      computeInPlace();
    }

    #ifdef EIGEN_PARSED_BY_DOXYGEN

    template<typename Rhs>
    inline const Solve<HouseholderQR, Rhs>
    solve(const MatrixBase<Rhs>& b) const;
    #endif

    HouseholderSequenceType householderQ() const
    {
      eigen_assert(m_isInitialized && "HouseholderQR is not initialized.");
      return HouseholderSequenceType(m_qr, m_hCoeffs.conjugate());
    }

    const MatrixType& matrixQR() const
    {
        eigen_assert(m_isInitialized && "HouseholderQR is not initialized.");
        return m_qr;
    }

    template<typename InputType>
    HouseholderQR& compute(const EigenBase<InputType>& matrix) {
      m_qr = matrix.derived();
      computeInPlace();
      return *this;
    }

    typename MatrixType::RealScalar absDeterminant() const;

    typename MatrixType::RealScalar logAbsDeterminant() const;

    inline Index rows() const { return m_qr.rows(); }
    inline Index cols() const { return m_qr.cols(); }

    const HCoeffsType& hCoeffs() const { return m_hCoeffs; }

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
    RowVectorType m_temp;
    bool m_isInitialized;
};

template<typename MatrixType>
typename MatrixType::RealScalar HouseholderQR<MatrixType>::absDeterminant() const
{
  using std::abs;
  eigen_assert(m_isInitialized && "HouseholderQR is not initialized.");
  eigen_assert(m_qr.rows() == m_qr.cols() && "You can't take the determinant of a non-square matrix!");
  return abs(m_qr.diagonal().prod());
}

template<typename MatrixType>
typename MatrixType::RealScalar HouseholderQR<MatrixType>::logAbsDeterminant() const
{
  eigen_assert(m_isInitialized && "HouseholderQR is not initialized.");
  eigen_assert(m_qr.rows() == m_qr.cols() && "You can't take the determinant of a non-square matrix!");
  return m_qr.diagonal().cwiseAbs().array().log().sum();
}

namespace internal {

template<typename MatrixQR, typename HCoeffs>
void householder_qr_inplace_unblocked(MatrixQR& mat, HCoeffs& hCoeffs, typename MatrixQR::Scalar* tempData = 0)
{
  typedef typename MatrixQR::Scalar Scalar;
  typedef typename MatrixQR::RealScalar RealScalar;
  Index rows = mat.rows();
  Index cols = mat.cols();
  Index size = (std::min)(rows,cols);

  eigen_assert(hCoeffs.size() == size);

  typedef Matrix<Scalar,MatrixQR::ColsAtCompileTime,1> TempType;
  TempType tempVector;
  if(tempData==0)
  {
    tempVector.resize(cols);
    tempData = tempVector.data();
  }

  for(Index k = 0; k < size; ++k)
  {
    Index remainingRows = rows - k;
    Index remainingCols = cols - k - 1;

    RealScalar beta;
    mat.col(k).tail(remainingRows).makeHouseholderInPlace(hCoeffs.coeffRef(k), beta);
    mat.coeffRef(k,k) = beta;

    mat.bottomRightCorner(remainingRows, remainingCols)
        .applyHouseholderOnTheLeft(mat.col(k).tail(remainingRows-1), hCoeffs.coeffRef(k), tempData+k+1);
  }
}

template<typename MatrixQR, typename HCoeffs,
  typename MatrixQRScalar = typename MatrixQR::Scalar,
  bool InnerStrideIsOne = (MatrixQR::InnerStrideAtCompileTime == 1 && HCoeffs::InnerStrideAtCompileTime == 1)>
struct householder_qr_inplace_blocked
{
  
  static void run(MatrixQR& mat, HCoeffs& hCoeffs, Index maxBlockSize=32,
      typename MatrixQR::Scalar* tempData = 0)
  {
    typedef typename MatrixQR::Scalar Scalar;
    typedef Block<MatrixQR,Dynamic,Dynamic> BlockType;

    Index rows = mat.rows();
    Index cols = mat.cols();
    Index size = (std::min)(rows, cols);

    typedef Matrix<Scalar,Dynamic,1,ColMajor,MatrixQR::MaxColsAtCompileTime,1> TempType;
    TempType tempVector;
    if(tempData==0)
    {
      tempVector.resize(cols);
      tempData = tempVector.data();
    }

    Index blockSize = (std::min)(maxBlockSize,size);

    Index k = 0;
    for (k = 0; k < size; k += blockSize)
    {
      Index bs = (std::min)(size-k,blockSize);  
      Index tcols = cols - k - bs;              
      Index brows = rows-k;                     

      BlockType A11_21 = mat.block(k,k,brows,bs);
      Block<HCoeffs,Dynamic,1> hCoeffsSegment = hCoeffs.segment(k,bs);

      householder_qr_inplace_unblocked(A11_21, hCoeffsSegment, tempData);

      if(tcols)
      {
        BlockType A21_22 = mat.block(k,k+bs,brows,tcols);
        apply_block_householder_on_the_left(A21_22,A11_21,hCoeffsSegment, false); 
      }
    }
  }
};

} 

#ifndef EIGEN_PARSED_BY_DOXYGEN
template<typename MatrixType_>
template<typename RhsType, typename DstType>
void HouseholderQR<MatrixType_>::_solve_impl(const RhsType &rhs, DstType &dst) const
{
  const Index rank = (std::min)(rows(), cols());

  typename RhsType::PlainObject c(rhs);

  c.applyOnTheLeft(householderQ().setLength(rank).adjoint() );

  m_qr.topLeftCorner(rank, rank)
      .template triangularView<Upper>()
      .solveInPlace(c.topRows(rank));

  dst.topRows(rank) = c.topRows(rank);
  dst.bottomRows(cols()-rank).setZero();
}

template<typename MatrixType_>
template<bool Conjugate, typename RhsType, typename DstType>
void HouseholderQR<MatrixType_>::_solve_impl_transposed(const RhsType &rhs, DstType &dst) const
{
  const Index rank = (std::min)(rows(), cols());

  typename RhsType::PlainObject c(rhs);

  m_qr.topLeftCorner(rank, rank)
      .template triangularView<Upper>()
      .transpose().template conjugateIf<Conjugate>()
      .solveInPlace(c.topRows(rank));

  dst.topRows(rank) = c.topRows(rank);
  dst.bottomRows(rows()-rank).setZero();

  dst.applyOnTheLeft(householderQ().setLength(rank).template conjugateIf<!Conjugate>() );
}
#endif

template<typename MatrixType>
void HouseholderQR<MatrixType>::computeInPlace()
{
  check_template_parameters();
  
  Index rows = m_qr.rows();
  Index cols = m_qr.cols();
  Index size = (std::min)(rows,cols);

  m_hCoeffs.resize(size);

  m_temp.resize(cols);

  internal::householder_qr_inplace_blocked<MatrixType, HCoeffsType>::run(m_qr, m_hCoeffs, 48, m_temp.data());

  m_isInitialized = true;
}

template<typename Derived>
const HouseholderQR<typename MatrixBase<Derived>::PlainObject>
MatrixBase<Derived>::householderQr() const
{
  return HouseholderQR<PlainObject>(eval());
}

} 

#endif 
