#ifndef EIGEN_COMPLETEORTHOGONALDECOMPOSITION_H
#define EIGEN_COMPLETEORTHOGONALDECOMPOSITION_H

namespace Eigen {

namespace internal {
template <typename MatrixType_>
struct traits<CompleteOrthogonalDecomposition<MatrixType_> >
    : traits<MatrixType_> {
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef int StorageIndex;
  enum { Flags = 0 };
};

}  

template <typename MatrixType_> class CompleteOrthogonalDecomposition
          : public SolverBase<CompleteOrthogonalDecomposition<MatrixType_> >
{
 public:
  typedef MatrixType_ MatrixType;
  typedef SolverBase<CompleteOrthogonalDecomposition> Base;

  template<typename Derived>
  friend struct internal::solve_assertion;

  EIGEN_GENERIC_PUBLIC_INTERFACE(CompleteOrthogonalDecomposition)
  enum {
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
  };
  typedef typename internal::plain_diag_type<MatrixType>::type HCoeffsType;
  typedef PermutationMatrix<ColsAtCompileTime, MaxColsAtCompileTime>
      PermutationType;
  typedef typename internal::plain_row_type<MatrixType, Index>::type
      IntRowVectorType;
  typedef typename internal::plain_row_type<MatrixType>::type RowVectorType;
  typedef typename internal::plain_row_type<MatrixType, RealScalar>::type
      RealRowVectorType;
  typedef HouseholderSequence<
      MatrixType, typename internal::remove_all<
                      typename HCoeffsType::ConjugateReturnType>::type>
      HouseholderSequenceType;
  typedef typename MatrixType::PlainObject PlainObject;

 private:
  typedef typename PermutationType::Index PermIndexType;

 public:

  CompleteOrthogonalDecomposition() : m_cpqr(), m_zCoeffs(), m_temp() {}

  CompleteOrthogonalDecomposition(Index rows, Index cols)
      : m_cpqr(rows, cols), m_zCoeffs((std::min)(rows, cols)), m_temp(cols) {}

  template <typename InputType>
  explicit CompleteOrthogonalDecomposition(const EigenBase<InputType>& matrix)
      : m_cpqr(matrix.rows(), matrix.cols()),
        m_zCoeffs((std::min)(matrix.rows(), matrix.cols())),
        m_temp(matrix.cols())
  {
    compute(matrix.derived());
  }

  template<typename InputType>
  explicit CompleteOrthogonalDecomposition(EigenBase<InputType>& matrix)
    : m_cpqr(matrix.derived()),
      m_zCoeffs((std::min)(matrix.rows(), matrix.cols())),
      m_temp(matrix.cols())
  {
    computeInPlace();
  } 

  #ifdef EIGEN_PARSED_BY_DOXYGEN

  template <typename Rhs>
  inline const Solve<CompleteOrthogonalDecomposition, Rhs> solve(
      const MatrixBase<Rhs>& b) const;
  #endif

  HouseholderSequenceType householderQ(void) const;
  HouseholderSequenceType matrixQ(void) const { return m_cpqr.householderQ(); }

  MatrixType matrixZ() const {
    MatrixType Z = MatrixType::Identity(m_cpqr.cols(), m_cpqr.cols());
    applyZOnTheLeftInPlace<false>(Z);
    return Z;
  }

  const MatrixType& matrixQTZ() const { return m_cpqr.matrixQR(); }

  const MatrixType& matrixT() const { return m_cpqr.matrixQR(); }

  template <typename InputType>
  CompleteOrthogonalDecomposition& compute(const EigenBase<InputType>& matrix) {
    
    m_cpqr.compute(matrix);
    computeInPlace();
    return *this;
  }

  const PermutationType& colsPermutation() const {
    return m_cpqr.colsPermutation();
  }

  typename MatrixType::RealScalar absDeterminant() const;

  typename MatrixType::RealScalar logAbsDeterminant() const;

  inline Index rank() const { return m_cpqr.rank(); }

  inline Index dimensionOfKernel() const { return m_cpqr.dimensionOfKernel(); }

  inline bool isInjective() const { return m_cpqr.isInjective(); }

  inline bool isSurjective() const { return m_cpqr.isSurjective(); }

  inline bool isInvertible() const { return m_cpqr.isInvertible(); }

  inline const Inverse<CompleteOrthogonalDecomposition> pseudoInverse() const
  {
    eigen_assert(m_cpqr.m_isInitialized && "CompleteOrthogonalDecomposition is not initialized.");
    return Inverse<CompleteOrthogonalDecomposition>(*this);
  }

  inline Index rows() const { return m_cpqr.rows(); }
  inline Index cols() const { return m_cpqr.cols(); }

  inline const HCoeffsType& hCoeffs() const { return m_cpqr.hCoeffs(); }

  const HCoeffsType& zCoeffs() const { return m_zCoeffs; }

  CompleteOrthogonalDecomposition& setThreshold(const RealScalar& threshold) {
    m_cpqr.setThreshold(threshold);
    return *this;
  }

  CompleteOrthogonalDecomposition& setThreshold(Default_t) {
    m_cpqr.setThreshold(Default);
    return *this;
  }

  RealScalar threshold() const { return m_cpqr.threshold(); }

  inline Index nonzeroPivots() const { return m_cpqr.nonzeroPivots(); }

  inline RealScalar maxPivot() const { return m_cpqr.maxPivot(); }

  ComputationInfo info() const {
    eigen_assert(m_cpqr.m_isInitialized && "Decomposition is not initialized.");
    return Success;
  }

#ifndef EIGEN_PARSED_BY_DOXYGEN
  template <typename RhsType, typename DstType>
  void _solve_impl(const RhsType& rhs, DstType& dst) const;

  template<bool Conjugate, typename RhsType, typename DstType>
  void _solve_impl_transposed(const RhsType &rhs, DstType &dst) const;
#endif

 protected:
  static void check_template_parameters() {
    EIGEN_STATIC_ASSERT_NON_INTEGER(Scalar);
  }

  template<bool Transpose_, typename Rhs>
  void _check_solve_assertion(const Rhs& b) const {
      EIGEN_ONLY_USED_FOR_DEBUG(b);
      eigen_assert(m_cpqr.m_isInitialized && "CompleteOrthogonalDecomposition is not initialized.");
      eigen_assert((Transpose_?derived().cols():derived().rows())==b.rows() && "CompleteOrthogonalDecomposition::solve(): invalid number of rows of the right hand side matrix b");
  }

  void computeInPlace();

  template <bool Conjugate, typename Rhs>
  void applyZOnTheLeftInPlace(Rhs& rhs) const;

  template <typename Rhs>
  void applyZAdjointOnTheLeftInPlace(Rhs& rhs) const;

  ColPivHouseholderQR<MatrixType> m_cpqr;
  HCoeffsType m_zCoeffs;
  RowVectorType m_temp;
};

template <typename MatrixType>
typename MatrixType::RealScalar
CompleteOrthogonalDecomposition<MatrixType>::absDeterminant() const {
  return m_cpqr.absDeterminant();
}

template <typename MatrixType>
typename MatrixType::RealScalar
CompleteOrthogonalDecomposition<MatrixType>::logAbsDeterminant() const {
  return m_cpqr.logAbsDeterminant();
}

template <typename MatrixType>
void CompleteOrthogonalDecomposition<MatrixType>::computeInPlace()
{
  check_template_parameters();

  eigen_assert(m_cpqr.cols() <= NumTraits<int>::highest());

  const Index rank = m_cpqr.rank();
  const Index cols = m_cpqr.cols();
  const Index rows = m_cpqr.rows();
  m_zCoeffs.resize((std::min)(rows, cols));
  m_temp.resize(cols);

  if (rank < cols) {

    for (Index k = rank - 1; k >= 0; --k) {
      if (k != rank - 1) {

        m_cpqr.m_qr.col(k).head(k + 1).swap(
            m_cpqr.m_qr.col(rank - 1).head(k + 1));
      }

      RealScalar beta;
      m_cpqr.m_qr.row(k)
          .tail(cols - rank + 1)
          .makeHouseholderInPlace(m_zCoeffs(k), beta);
      m_cpqr.m_qr(k, rank - 1) = beta;
      if (k > 0) {
        
        m_cpqr.m_qr.topRightCorner(k, cols - rank + 1)
            .applyHouseholderOnTheRight(
                m_cpqr.m_qr.row(k).tail(cols - rank).adjoint(), m_zCoeffs(k),
                &m_temp(0));
      }
      if (k != rank - 1) {
        
        m_cpqr.m_qr.col(k).head(k + 1).swap(
            m_cpqr.m_qr.col(rank - 1).head(k + 1));
      }
    }
  }
}

template <typename MatrixType>
template <bool Conjugate, typename Rhs>
void CompleteOrthogonalDecomposition<MatrixType>::applyZOnTheLeftInPlace(
    Rhs& rhs) const {
  const Index cols = this->cols();
  const Index nrhs = rhs.cols();
  const Index rank = this->rank();
  Matrix<typename Rhs::Scalar, Dynamic, 1> temp((std::max)(cols, nrhs));
  for (Index k = rank-1; k >= 0; --k) {
    if (k != rank - 1) {
      rhs.row(k).swap(rhs.row(rank - 1));
    }
    rhs.middleRows(rank - 1, cols - rank + 1)
        .applyHouseholderOnTheLeft(
            matrixQTZ().row(k).tail(cols - rank).transpose().template conjugateIf<!Conjugate>(), zCoeffs().template conjugateIf<Conjugate>()(k),
            &temp(0));
    if (k != rank - 1) {
      rhs.row(k).swap(rhs.row(rank - 1));
    }
  }
}

template <typename MatrixType>
template <typename Rhs>
void CompleteOrthogonalDecomposition<MatrixType>::applyZAdjointOnTheLeftInPlace(
    Rhs& rhs) const {
  const Index cols = this->cols();
  const Index nrhs = rhs.cols();
  const Index rank = this->rank();
  Matrix<typename Rhs::Scalar, Dynamic, 1> temp((std::max)(cols, nrhs));
  for (Index k = 0; k < rank; ++k) {
    if (k != rank - 1) {
      rhs.row(k).swap(rhs.row(rank - 1));
    }
    rhs.middleRows(rank - 1, cols - rank + 1)
        .applyHouseholderOnTheLeft(
            matrixQTZ().row(k).tail(cols - rank).adjoint(), zCoeffs()(k),
            &temp(0));
    if (k != rank - 1) {
      rhs.row(k).swap(rhs.row(rank - 1));
    }
  }
}

#ifndef EIGEN_PARSED_BY_DOXYGEN
template <typename MatrixType_>
template <typename RhsType, typename DstType>
void CompleteOrthogonalDecomposition<MatrixType_>::_solve_impl(
    const RhsType& rhs, DstType& dst) const {
  const Index rank = this->rank();
  if (rank == 0) {
    dst.setZero();
    return;
  }

  typename RhsType::PlainObject c(rhs);
  c.applyOnTheLeft(matrixQ().setLength(rank).adjoint());

  dst.topRows(rank) = matrixT()
                          .topLeftCorner(rank, rank)
                          .template triangularView<Upper>()
                          .solve(c.topRows(rank));

  const Index cols = this->cols();
  if (rank < cols) {

    dst.bottomRows(cols - rank).setZero();
    applyZAdjointOnTheLeftInPlace(dst);
  }

  dst = colsPermutation() * dst;
}

template<typename MatrixType_>
template<bool Conjugate, typename RhsType, typename DstType>
void CompleteOrthogonalDecomposition<MatrixType_>::_solve_impl_transposed(const RhsType &rhs, DstType &dst) const
{
  const Index rank = this->rank();

  if (rank == 0) {
    dst.setZero();
    return;
  }

  typename RhsType::PlainObject c(colsPermutation().transpose()*rhs);

  if (rank < cols()) {
    applyZOnTheLeftInPlace<!Conjugate>(c);
  }

  matrixT().topLeftCorner(rank, rank)
           .template triangularView<Upper>()
           .transpose().template conjugateIf<Conjugate>()
           .solveInPlace(c.topRows(rank));

  dst.topRows(rank) = c.topRows(rank);
  dst.bottomRows(rows()-rank).setZero();

  dst.applyOnTheLeft(householderQ().setLength(rank).template conjugateIf<!Conjugate>() );
}
#endif

namespace internal {

template<typename MatrixType>
struct traits<Inverse<CompleteOrthogonalDecomposition<MatrixType> > >
  : traits<typename Transpose<typename MatrixType::PlainObject>::PlainObject>
{
  enum { Flags = 0 };
};

template<typename DstXprType, typename MatrixType>
struct Assignment<DstXprType, Inverse<CompleteOrthogonalDecomposition<MatrixType> >, internal::assign_op<typename DstXprType::Scalar,typename CompleteOrthogonalDecomposition<MatrixType>::Scalar>, Dense2Dense>
{
  typedef CompleteOrthogonalDecomposition<MatrixType> CodType;
  typedef Inverse<CodType> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<typename DstXprType::Scalar,typename CodType::Scalar> &)
  {
    typedef Matrix<typename CodType::Scalar, CodType::RowsAtCompileTime, CodType::RowsAtCompileTime, 0, CodType::MaxRowsAtCompileTime, CodType::MaxRowsAtCompileTime> IdentityMatrixType;
    dst = src.nestedExpression().solve(IdentityMatrixType::Identity(src.cols(), src.cols()));
  }
};

} 

template <typename MatrixType>
typename CompleteOrthogonalDecomposition<MatrixType>::HouseholderSequenceType
CompleteOrthogonalDecomposition<MatrixType>::householderQ() const {
  return m_cpqr.householderQ();
}

template <typename Derived>
const CompleteOrthogonalDecomposition<typename MatrixBase<Derived>::PlainObject>
MatrixBase<Derived>::completeOrthogonalDecomposition() const {
  return CompleteOrthogonalDecomposition<PlainObject>(eval());
}

}  

#endif  
