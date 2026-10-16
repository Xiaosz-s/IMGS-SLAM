#ifndef EIGEN_SVDBASE_H
#define EIGEN_SVDBASE_H

namespace Eigen {

namespace internal {
template<typename Derived> struct traits<SVDBase<Derived> >
 : traits<Derived>
{
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef int StorageIndex;
  enum { Flags = 0 };
};
}

template<typename Derived> class SVDBase
 : public SolverBase<SVDBase<Derived> >
{
public: 
   
  template<typename Derived_>
  friend struct internal::solve_assertion;

  typedef typename internal::traits<Derived>::MatrixType MatrixType;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<typename MatrixType::Scalar>::Real RealScalar;
  typedef typename Eigen::internal::traits<SVDBase>::StorageIndex StorageIndex;
  typedef Eigen::Index Index; 
  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime,
    DiagSizeAtCompileTime = EIGEN_SIZE_MIN_PREFER_DYNAMIC(RowsAtCompileTime,ColsAtCompileTime),
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime,
    MaxDiagSizeAtCompileTime = EIGEN_SIZE_MIN_PREFER_FIXED(MaxRowsAtCompileTime,MaxColsAtCompileTime),
    MatrixOptions = MatrixType::Options
  };

  typedef Matrix<Scalar, RowsAtCompileTime, RowsAtCompileTime, MatrixOptions, MaxRowsAtCompileTime, MaxRowsAtCompileTime> MatrixUType;
  typedef Matrix<Scalar, ColsAtCompileTime, ColsAtCompileTime, MatrixOptions, MaxColsAtCompileTime, MaxColsAtCompileTime> MatrixVType;
  typedef typename internal::plain_diag_type<MatrixType, RealScalar>::type SingularValuesType;
  
  Derived& derived() { return *static_cast<Derived*>(this); }
  const Derived& derived() const { return *static_cast<const Derived*>(this); }

  const MatrixUType& matrixU() const
  {
    _check_compute_assertions();
    eigen_assert(computeU() && "This SVD decomposition didn't compute U. Did you ask for it?");
    return m_matrixU;
  }

  const MatrixVType& matrixV() const
  {
    _check_compute_assertions();
    eigen_assert(computeV() && "This SVD decomposition didn't compute V. Did you ask for it?");
    return m_matrixV;
  }

  const SingularValuesType& singularValues() const
  {
    _check_compute_assertions();
    return m_singularValues;
  }

  Index nonzeroSingularValues() const
  {
    _check_compute_assertions();
    return m_nonzeroSingularValues;
  }

  inline Index rank() const
  {
    using std::abs;
    _check_compute_assertions();
    if(m_singularValues.size()==0) return 0;
    RealScalar premultiplied_threshold = numext::maxi<RealScalar>(m_singularValues.coeff(0) * threshold(), (std::numeric_limits<RealScalar>::min)());
    Index i = m_nonzeroSingularValues-1;
    while(i>=0 && m_singularValues.coeff(i) < premultiplied_threshold) --i;
    return i+1;
  }

  Derived& setThreshold(const RealScalar& threshold)
  {
    m_usePrescribedThreshold = true;
    m_prescribedThreshold = threshold;
    return derived();
  }

  Derived& setThreshold(Default_t)
  {
    m_usePrescribedThreshold = false;
    return derived();
  }

  RealScalar threshold() const
  {
    eigen_assert(m_isInitialized || m_usePrescribedThreshold);
    
    Index diagSize = (std::max<Index>)(1,m_diagSize);
    return m_usePrescribedThreshold ? m_prescribedThreshold
                                    : RealScalar(diagSize)*NumTraits<Scalar>::epsilon();
  }

  inline bool computeU() const { return m_computeFullU || m_computeThinU; }
  
  inline bool computeV() const { return m_computeFullV || m_computeThinV; }

  inline Index rows() const { return m_rows; }
  inline Index cols() const { return m_cols; }
  
  #ifdef EIGEN_PARSED_BY_DOXYGEN

  template<typename Rhs>
  inline const Solve<Derived, Rhs>
  solve(const MatrixBase<Rhs>& b) const;
  #endif

  EIGEN_DEVICE_FUNC
  ComputationInfo info() const
  {
    eigen_assert(m_isInitialized && "SVD is not initialized.");
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

  void _check_compute_assertions() const {
    eigen_assert(m_isInitialized && "SVD is not initialized.");
  }

  template<bool Transpose_, typename Rhs>
  void _check_solve_assertion(const Rhs& b) const {
      EIGEN_ONLY_USED_FOR_DEBUG(b);
      _check_compute_assertions();
      eigen_assert(computeU() && computeV() && "SVDBase::solve(): Both unitaries U and V are required to be computed (thin unitaries suffice).");
      eigen_assert((Transpose_?cols():rows())==b.rows() && "SVDBase::solve(): invalid number of rows of the right hand side matrix b");
  }

  bool allocate(Index rows, Index cols, unsigned int computationOptions) ;

  MatrixUType m_matrixU;
  MatrixVType m_matrixV;
  SingularValuesType m_singularValues;
  ComputationInfo m_info;
  bool m_isInitialized, m_isAllocated, m_usePrescribedThreshold;
  bool m_computeFullU, m_computeThinU;
  bool m_computeFullV, m_computeThinV;
  unsigned int m_computationOptions;
  Index m_nonzeroSingularValues, m_rows, m_cols, m_diagSize;
  RealScalar m_prescribedThreshold;

  SVDBase()
    : m_info(Success),
      m_isInitialized(false),
      m_isAllocated(false),
      m_usePrescribedThreshold(false),
      m_computeFullU(false),
      m_computeThinU(false),
      m_computeFullV(false),
      m_computeThinV(false),
      m_computationOptions(0),
      m_rows(-1), m_cols(-1), m_diagSize(0)
  {
    check_template_parameters();
  }

};

#ifndef EIGEN_PARSED_BY_DOXYGEN
template<typename Derived>
template<typename RhsType, typename DstType>
void SVDBase<Derived>::_solve_impl(const RhsType &rhs, DstType &dst) const
{

  Matrix<typename RhsType::Scalar, Dynamic, RhsType::ColsAtCompileTime, 0, MatrixType::MaxRowsAtCompileTime, RhsType::MaxColsAtCompileTime> tmp;
  Index l_rank = rank();
  tmp.noalias() =  m_matrixU.leftCols(l_rank).adjoint() * rhs;
  tmp = m_singularValues.head(l_rank).asDiagonal().inverse() * tmp;
  dst = m_matrixV.leftCols(l_rank) * tmp;
}

template<typename Derived>
template<bool Conjugate, typename RhsType, typename DstType>
void SVDBase<Derived>::_solve_impl_transposed(const RhsType &rhs, DstType &dst) const
{

  Matrix<typename RhsType::Scalar, Dynamic, RhsType::ColsAtCompileTime, 0, MatrixType::MaxRowsAtCompileTime, RhsType::MaxColsAtCompileTime> tmp;
  Index l_rank = rank();

  tmp.noalias() =  m_matrixV.leftCols(l_rank).transpose().template conjugateIf<Conjugate>() * rhs;
  tmp = m_singularValues.head(l_rank).asDiagonal().inverse() * tmp;
  dst = m_matrixU.template conjugateIf<!Conjugate>().leftCols(l_rank) * tmp;
}
#endif

template<typename MatrixType>
bool SVDBase<MatrixType>::allocate(Index rows, Index cols, unsigned int computationOptions)
{
  eigen_assert(rows >= 0 && cols >= 0);

  if (m_isAllocated &&
      rows == m_rows &&
      cols == m_cols &&
      computationOptions == m_computationOptions)
  {
    return true;
  }

  m_rows = rows;
  m_cols = cols;
  m_info = Success;
  m_isInitialized = false;
  m_isAllocated = true;
  m_computationOptions = computationOptions;
  m_computeFullU = (computationOptions & ComputeFullU) != 0;
  m_computeThinU = (computationOptions & ComputeThinU) != 0;
  m_computeFullV = (computationOptions & ComputeFullV) != 0;
  m_computeThinV = (computationOptions & ComputeThinV) != 0;
  eigen_assert(!(m_computeFullU && m_computeThinU) && "SVDBase: you can't ask for both full and thin U");
  eigen_assert(!(m_computeFullV && m_computeThinV) && "SVDBase: you can't ask for both full and thin V");
  eigen_assert(EIGEN_IMPLIES(m_computeThinU || m_computeThinV, MatrixType::ColsAtCompileTime==Dynamic) &&
	       "SVDBase: thin U and V are only available when your matrix has a dynamic number of columns.");

  m_diagSize = (std::min)(m_rows, m_cols);
  m_singularValues.resize(m_diagSize);
  if(RowsAtCompileTime==Dynamic)
    m_matrixU.resize(m_rows, m_computeFullU ? m_rows : m_computeThinU ? m_diagSize : 0);
  if(ColsAtCompileTime==Dynamic)
    m_matrixV.resize(m_cols, m_computeFullV ? m_cols : m_computeThinV ? m_diagSize : 0);

  return false;
}

}

#endif 
