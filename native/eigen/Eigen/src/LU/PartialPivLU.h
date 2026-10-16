#ifndef EIGEN_PARTIALLU_H
#define EIGEN_PARTIALLU_H

namespace Eigen {

namespace internal {
template<typename MatrixType_> struct traits<PartialPivLU<MatrixType_> >
 : traits<MatrixType_>
{
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef int StorageIndex;
  typedef traits<MatrixType_> BaseTraits;
  enum {
    Flags = BaseTraits::Flags & RowMajorBit,
    CoeffReadCost = Dynamic
  };
};

template<typename T,typename Derived>
struct enable_if_ref;

template<typename T,typename Derived>
struct enable_if_ref<Ref<T>,Derived> {
  typedef Derived type;
};

} 

template<typename MatrixType_> class PartialPivLU
  : public SolverBase<PartialPivLU<MatrixType_> >
{
  public:

    typedef MatrixType_ MatrixType;
    typedef SolverBase<PartialPivLU> Base;
    friend class SolverBase<PartialPivLU>;

    EIGEN_GENERIC_PUBLIC_INTERFACE(PartialPivLU)
    enum {
      MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };
    typedef PermutationMatrix<RowsAtCompileTime, MaxRowsAtCompileTime> PermutationType;
    typedef Transpositions<RowsAtCompileTime, MaxRowsAtCompileTime> TranspositionType;
    typedef typename MatrixType::PlainObject PlainObject;

    PartialPivLU();

    explicit PartialPivLU(Index size);

    template<typename InputType>
    explicit PartialPivLU(const EigenBase<InputType>& matrix);

    template<typename InputType>
    explicit PartialPivLU(EigenBase<InputType>& matrix);

    template<typename InputType>
    PartialPivLU& compute(const EigenBase<InputType>& matrix) {
      m_lu = matrix.derived();
      compute();
      return *this;
    }

    inline const MatrixType& matrixLU() const
    {
      eigen_assert(m_isInitialized && "PartialPivLU is not initialized.");
      return m_lu;
    }

    inline const PermutationType& permutationP() const
    {
      eigen_assert(m_isInitialized && "PartialPivLU is not initialized.");
      return m_p;
    }

    #ifdef EIGEN_PARSED_BY_DOXYGEN

    template<typename Rhs>
    inline const Solve<PartialPivLU, Rhs>
    solve(const MatrixBase<Rhs>& b) const;
    #endif

    inline RealScalar rcond() const
    {
      eigen_assert(m_isInitialized && "PartialPivLU is not initialized.");
      return internal::rcond_estimate_helper(m_l1_norm, *this);
    }

    inline const Inverse<PartialPivLU> inverse() const
    {
      eigen_assert(m_isInitialized && "PartialPivLU is not initialized.");
      return Inverse<PartialPivLU>(*this);
    }

    Scalar determinant() const;

    MatrixType reconstructedMatrix() const;

    EIGEN_CONSTEXPR inline Index rows() const EIGEN_NOEXCEPT { return m_lu.rows(); }
    EIGEN_CONSTEXPR inline Index cols() const EIGEN_NOEXCEPT { return m_lu.cols(); }

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    template<typename RhsType, typename DstType>
    EIGEN_DEVICE_FUNC
    void _solve_impl(const RhsType &rhs, DstType &dst) const {

      dst = permutationP() * rhs;

      m_lu.template triangularView<UnitLower>().solveInPlace(dst);

      m_lu.template triangularView<Upper>().solveInPlace(dst);
    }

    template<bool Conjugate, typename RhsType, typename DstType>
    EIGEN_DEVICE_FUNC
    void _solve_impl_transposed(const RhsType &rhs, DstType &dst) const {

      eigen_assert(rhs.rows() == m_lu.cols());

      dst = m_lu.template triangularView<Upper>().transpose()
                .template conjugateIf<Conjugate>().solve(rhs);
      
      m_lu.template triangularView<UnitLower>().transpose()
          .template conjugateIf<Conjugate>().solveInPlace(dst);
      
      dst = permutationP().transpose() * dst;
    }
    #endif

  protected:

    static void check_template_parameters()
    {
      EIGEN_STATIC_ASSERT_NON_INTEGER(Scalar);
    }

    void compute();

    MatrixType m_lu;
    PermutationType m_p;
    TranspositionType m_rowsTranspositions;
    RealScalar m_l1_norm;
    signed char m_det_p;
    bool m_isInitialized;
};

template<typename MatrixType>
PartialPivLU<MatrixType>::PartialPivLU()
  : m_lu(),
    m_p(),
    m_rowsTranspositions(),
    m_l1_norm(0),
    m_det_p(0),
    m_isInitialized(false)
{
}

template<typename MatrixType>
PartialPivLU<MatrixType>::PartialPivLU(Index size)
  : m_lu(size, size),
    m_p(size),
    m_rowsTranspositions(size),
    m_l1_norm(0),
    m_det_p(0),
    m_isInitialized(false)
{
}

template<typename MatrixType>
template<typename InputType>
PartialPivLU<MatrixType>::PartialPivLU(const EigenBase<InputType>& matrix)
  : m_lu(matrix.rows(),matrix.cols()),
    m_p(matrix.rows()),
    m_rowsTranspositions(matrix.rows()),
    m_l1_norm(0),
    m_det_p(0),
    m_isInitialized(false)
{
  compute(matrix.derived());
}

template<typename MatrixType>
template<typename InputType>
PartialPivLU<MatrixType>::PartialPivLU(EigenBase<InputType>& matrix)
  : m_lu(matrix.derived()),
    m_p(matrix.rows()),
    m_rowsTranspositions(matrix.rows()),
    m_l1_norm(0),
    m_det_p(0),
    m_isInitialized(false)
{
  compute();
}

namespace internal {

template<typename Scalar, int StorageOrder, typename PivIndex, int SizeAtCompileTime=Dynamic>
struct partial_lu_impl
{
  static const int UnBlockedBound = 16;
  static const bool UnBlockedAtCompileTime = SizeAtCompileTime!=Dynamic && SizeAtCompileTime<=UnBlockedBound;
  static const int ActualSizeAtCompileTime = UnBlockedAtCompileTime ? SizeAtCompileTime : Dynamic;
  
  static const int RRows = SizeAtCompileTime==2 ? 1 : Dynamic;
  static const int RCols = SizeAtCompileTime==2 ? 1 : Dynamic;
  typedef Matrix<Scalar, ActualSizeAtCompileTime, ActualSizeAtCompileTime, StorageOrder> MatrixType;
  typedef Ref<MatrixType> MatrixTypeRef;
  typedef Ref<Matrix<Scalar, Dynamic, Dynamic, StorageOrder> > BlockType;
  typedef typename MatrixType::RealScalar RealScalar;

  static Index unblocked_lu(MatrixTypeRef& lu, PivIndex* row_transpositions, PivIndex& nb_transpositions)
  {
    typedef scalar_score_coeff_op<Scalar> Scoring;
    typedef typename Scoring::result_type Score;
    const Index rows = lu.rows();
    const Index cols = lu.cols();
    const Index size = (std::min)(rows,cols);

    const Index endk = UnBlockedAtCompileTime ? size-1 : size;
    nb_transpositions = 0;
    Index first_zero_pivot = -1;
    for(Index k = 0; k < endk; ++k)
    {
      int rrows = internal::convert_index<int>(rows-k-1);
      int rcols = internal::convert_index<int>(cols-k-1);

      Index row_of_biggest_in_col;
      Score biggest_in_corner
        = lu.col(k).tail(rows-k).unaryExpr(Scoring()).maxCoeff(&row_of_biggest_in_col);
      row_of_biggest_in_col += k;

      row_transpositions[k] = PivIndex(row_of_biggest_in_col);

      if(biggest_in_corner != Score(0))
      {
        if(k != row_of_biggest_in_col)
        {
          lu.row(k).swap(lu.row(row_of_biggest_in_col));
          ++nb_transpositions;
        }

        lu.col(k).tail(fix<RRows>(rrows)) /= lu.coeff(k,k);
      }
      else if(first_zero_pivot==-1)
      {

        first_zero_pivot = k;
      }

      if(k<rows-1)
        lu.bottomRightCorner(fix<RRows>(rrows),fix<RCols>(rcols)).noalias() -= lu.col(k).tail(fix<RRows>(rrows)) * lu.row(k).tail(fix<RCols>(rcols));
    }

    if(UnBlockedAtCompileTime)
    {
      Index k = endk;
      row_transpositions[k] = PivIndex(k);
      if (Scoring()(lu(k, k)) == Score(0) && first_zero_pivot == -1)
        first_zero_pivot = k;
    }

    return first_zero_pivot;
  }

  static Index blocked_lu(Index rows, Index cols, Scalar* lu_data, Index luStride, PivIndex* row_transpositions, PivIndex& nb_transpositions, Index maxBlockSize=256)
  {
    MatrixTypeRef lu = MatrixType::Map(lu_data,rows, cols, OuterStride<>(luStride));

    const Index size = (std::min)(rows,cols);

    if(UnBlockedAtCompileTime || size<=UnBlockedBound)
    {
      return unblocked_lu(lu, row_transpositions, nb_transpositions);
    }

    Index blockSize;
    {
      blockSize = size/8;
      blockSize = (blockSize/16)*16;
      blockSize = (std::min)((std::max)(blockSize,Index(8)), maxBlockSize);
    }

    nb_transpositions = 0;
    Index first_zero_pivot = -1;
    for(Index k = 0; k < size; k+=blockSize)
    {
      Index bs = (std::min)(size-k,blockSize); 
      Index trows = rows - k - bs; 
      Index tsize = size - k - bs; 

      BlockType A_0 = lu.block(0,0,rows,k);
      BlockType A_2 = lu.block(0,k+bs,rows,tsize);
      BlockType A11 = lu.block(k,k,bs,bs);
      BlockType A12 = lu.block(k,k+bs,bs,tsize);
      BlockType A21 = lu.block(k+bs,k,trows,bs);
      BlockType A22 = lu.block(k+bs,k+bs,trows,tsize);

      PivIndex nb_transpositions_in_panel;

      Index ret = blocked_lu(trows+bs, bs, &lu.coeffRef(k,k), luStride,
                   row_transpositions+k, nb_transpositions_in_panel, 16);
      if(ret>=0 && first_zero_pivot==-1)
        first_zero_pivot = k+ret;

      nb_transpositions += nb_transpositions_in_panel;
      
      for(Index i=k; i<k+bs; ++i)
      {
        Index piv = (row_transpositions[i] += internal::convert_index<PivIndex>(k));
        A_0.row(i).swap(A_0.row(piv));
      }

      if(trows)
      {
        
        for(Index i=k;i<k+bs; ++i)
          A_2.row(i).swap(A_2.row(row_transpositions[i]));

        A11.template triangularView<UnitLower>().solveInPlace(A12);

        A22.noalias() -= A21 * A12;
      }
    }
    return first_zero_pivot;
  }
};

template<typename MatrixType, typename TranspositionType>
void partial_lu_inplace(MatrixType& lu, TranspositionType& row_transpositions, typename TranspositionType::StorageIndex& nb_transpositions)
{
  
  if (lu.rows() == 0 || lu.cols() == 0) {
    nb_transpositions = 0;
    return;
  }
  eigen_assert(lu.cols() == row_transpositions.size());
  eigen_assert(row_transpositions.size() < 2 || (&row_transpositions.coeffRef(1)-&row_transpositions.coeffRef(0)) == 1);

  partial_lu_impl
    < typename MatrixType::Scalar, MatrixType::Flags&RowMajorBit?RowMajor:ColMajor,
      typename TranspositionType::StorageIndex,
      EIGEN_SIZE_MIN_PREFER_FIXED(MatrixType::RowsAtCompileTime,MatrixType::ColsAtCompileTime)>
    ::blocked_lu(lu.rows(), lu.cols(), &lu.coeffRef(0,0), lu.outerStride(), &row_transpositions.coeffRef(0), nb_transpositions);
}

} 

template<typename MatrixType>
void PartialPivLU<MatrixType>::compute()
{
  check_template_parameters();

  eigen_assert(m_lu.rows()<NumTraits<int>::highest());

  if(m_lu.cols()>0)
    m_l1_norm = m_lu.cwiseAbs().colwise().sum().maxCoeff();
  else
    m_l1_norm = RealScalar(0);

  eigen_assert(m_lu.rows() == m_lu.cols() && "PartialPivLU is only for square (and moreover invertible) matrices");
  const Index size = m_lu.rows();

  m_rowsTranspositions.resize(size);

  typename TranspositionType::StorageIndex nb_transpositions;
  internal::partial_lu_inplace(m_lu, m_rowsTranspositions, nb_transpositions);
  m_det_p = (nb_transpositions%2) ? -1 : 1;

  m_p = m_rowsTranspositions;

  m_isInitialized = true;
}

template<typename MatrixType>
typename PartialPivLU<MatrixType>::Scalar PartialPivLU<MatrixType>::determinant() const
{
  eigen_assert(m_isInitialized && "PartialPivLU is not initialized.");
  return Scalar(m_det_p) * m_lu.diagonal().prod();
}

template<typename MatrixType>
MatrixType PartialPivLU<MatrixType>::reconstructedMatrix() const
{
  eigen_assert(m_isInitialized && "LU is not initialized.");
  
  MatrixType res = m_lu.template triangularView<UnitLower>().toDenseMatrix()
                 * m_lu.template triangularView<Upper>();

  res = m_p.inverse() * res;

  return res;
}

namespace internal {

template<typename DstXprType, typename MatrixType>
struct Assignment<DstXprType, Inverse<PartialPivLU<MatrixType> >, internal::assign_op<typename DstXprType::Scalar,typename PartialPivLU<MatrixType>::Scalar>, Dense2Dense>
{
  typedef PartialPivLU<MatrixType> LuType;
  typedef Inverse<LuType> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<typename DstXprType::Scalar,typename LuType::Scalar> &)
  {
    dst = src.nestedExpression().solve(MatrixType::Identity(src.rows(), src.cols()));
  }
};
} 

template<typename Derived>
inline const PartialPivLU<typename MatrixBase<Derived>::PlainObject>
MatrixBase<Derived>::partialPivLu() const
{
  return PartialPivLU<PlainObject>(eval());
}

template<typename Derived>
inline const PartialPivLU<typename MatrixBase<Derived>::PlainObject>
MatrixBase<Derived>::lu() const
{
  return PartialPivLU<PlainObject>(eval());
}

} 

#endif 
