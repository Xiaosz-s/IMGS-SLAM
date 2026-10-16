#ifndef EIGEN_DIAGONALMATRIX_H
#define EIGEN_DIAGONALMATRIX_H

namespace Eigen { 

#ifndef EIGEN_PARSED_BY_DOXYGEN
template<typename Derived>
class DiagonalBase : public EigenBase<Derived>
{
  public:
    typedef typename internal::traits<Derived>::DiagonalVectorType DiagonalVectorType;
    typedef typename DiagonalVectorType::Scalar Scalar;
    typedef typename DiagonalVectorType::RealScalar RealScalar;
    typedef typename internal::traits<Derived>::StorageKind StorageKind;
    typedef typename internal::traits<Derived>::StorageIndex StorageIndex;

    enum {
      RowsAtCompileTime = DiagonalVectorType::SizeAtCompileTime,
      ColsAtCompileTime = DiagonalVectorType::SizeAtCompileTime,
      MaxRowsAtCompileTime = DiagonalVectorType::MaxSizeAtCompileTime,
      MaxColsAtCompileTime = DiagonalVectorType::MaxSizeAtCompileTime,
      IsVectorAtCompileTime = 0,
      Flags = NoPreferredStorageOrderBit
    };

    typedef Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime, 0, MaxRowsAtCompileTime, MaxColsAtCompileTime> DenseMatrixType;
    typedef DenseMatrixType DenseType;
    typedef DiagonalMatrix<Scalar,DiagonalVectorType::SizeAtCompileTime,DiagonalVectorType::MaxSizeAtCompileTime> PlainObject;

    EIGEN_DEVICE_FUNC
    inline const Derived& derived() const { return *static_cast<const Derived*>(this); }
    EIGEN_DEVICE_FUNC
    inline Derived& derived() { return *static_cast<Derived*>(this); }

    EIGEN_DEVICE_FUNC
    DenseMatrixType toDenseMatrix() const { return derived(); }

    EIGEN_DEVICE_FUNC
    inline const DiagonalVectorType& diagonal() const { return derived().diagonal(); }
    EIGEN_DEVICE_FUNC
    inline DiagonalVectorType& diagonal() { return derived().diagonal(); }

    EIGEN_DEVICE_FUNC
    inline Index rows() const { return diagonal().size(); }
    EIGEN_DEVICE_FUNC
    inline Index cols() const { return diagonal().size(); }

    template<typename MatrixDerived>
    EIGEN_DEVICE_FUNC
    const Product<Derived,MatrixDerived,LazyProduct>
    operator*(const MatrixBase<MatrixDerived> &matrix) const
    {
      return Product<Derived, MatrixDerived, LazyProduct>(derived(),matrix.derived());
    }

    typedef DiagonalWrapper<const CwiseUnaryOp<internal::scalar_inverse_op<Scalar>, const DiagonalVectorType> > InverseReturnType;
    EIGEN_DEVICE_FUNC
    inline const InverseReturnType
    inverse() const
    {
      return InverseReturnType(diagonal().cwiseInverse());
    }
    
    EIGEN_DEVICE_FUNC
    inline const DiagonalWrapper<const EIGEN_EXPR_BINARYOP_SCALAR_RETURN_TYPE(DiagonalVectorType,Scalar,product) >
    operator*(const Scalar& scalar) const
    {
      return DiagonalWrapper<const EIGEN_EXPR_BINARYOP_SCALAR_RETURN_TYPE(DiagonalVectorType,Scalar,product) >(diagonal() * scalar);
    }
    EIGEN_DEVICE_FUNC
    friend inline const DiagonalWrapper<const EIGEN_SCALAR_BINARYOP_EXPR_RETURN_TYPE(Scalar,DiagonalVectorType,product) >
    operator*(const Scalar& scalar, const DiagonalBase& other)
    {
      return DiagonalWrapper<const EIGEN_SCALAR_BINARYOP_EXPR_RETURN_TYPE(Scalar,DiagonalVectorType,product) >(scalar * other.diagonal());
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    #ifdef EIGEN_PARSED_BY_DOXYGEN
    inline unspecified_expression_type
    #else
    inline const DiagonalWrapper<const EIGEN_CWISE_BINARY_RETURN_TYPE(DiagonalVectorType,typename OtherDerived::DiagonalVectorType,sum) >
    #endif
    operator+(const DiagonalBase<OtherDerived>& other) const
    {
      return (diagonal() + other.diagonal()).asDiagonal();
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    #ifdef EIGEN_PARSED_BY_DOXYGEN
    inline unspecified_expression_type
    #else
    inline const DiagonalWrapper<const EIGEN_CWISE_BINARY_RETURN_TYPE(DiagonalVectorType,typename OtherDerived::DiagonalVectorType,difference) >
    #endif
    operator-(const DiagonalBase<OtherDerived>& other) const
    {
      return (diagonal() - other.diagonal()).asDiagonal();
    }
};

#endif

namespace internal {
template<typename Scalar_, int SizeAtCompileTime, int MaxSizeAtCompileTime>
struct traits<DiagonalMatrix<Scalar_,SizeAtCompileTime,MaxSizeAtCompileTime> >
 : traits<Matrix<Scalar_,SizeAtCompileTime,SizeAtCompileTime,0,MaxSizeAtCompileTime,MaxSizeAtCompileTime> >
{
  typedef Matrix<Scalar_,SizeAtCompileTime,1,0,MaxSizeAtCompileTime,1> DiagonalVectorType;
  typedef DiagonalShape StorageKind;
  enum {
    Flags = LvalueBit | NoPreferredStorageOrderBit
  };
};
}
template<typename Scalar_, int SizeAtCompileTime, int MaxSizeAtCompileTime>
class DiagonalMatrix
  : public DiagonalBase<DiagonalMatrix<Scalar_,SizeAtCompileTime,MaxSizeAtCompileTime> >
{
  public:
    #ifndef EIGEN_PARSED_BY_DOXYGEN
    typedef typename internal::traits<DiagonalMatrix>::DiagonalVectorType DiagonalVectorType;
    typedef const DiagonalMatrix& Nested;
    typedef Scalar_ Scalar;
    typedef typename internal::traits<DiagonalMatrix>::StorageKind StorageKind;
    typedef typename internal::traits<DiagonalMatrix>::StorageIndex StorageIndex;
    #endif

  protected:

    DiagonalVectorType m_diagonal;

  public:

    EIGEN_DEVICE_FUNC
    inline const DiagonalVectorType& diagonal() const { return m_diagonal; }
    
    EIGEN_DEVICE_FUNC
    inline DiagonalVectorType& diagonal() { return m_diagonal; }

    EIGEN_DEVICE_FUNC
    inline DiagonalMatrix() {}

    EIGEN_DEVICE_FUNC
    explicit inline DiagonalMatrix(Index dim) : m_diagonal(dim) {}

    EIGEN_DEVICE_FUNC
    inline DiagonalMatrix(const Scalar& x, const Scalar& y) : m_diagonal(x,y) {}

    EIGEN_DEVICE_FUNC
    inline DiagonalMatrix(const Scalar& x, const Scalar& y, const Scalar& z) : m_diagonal(x,y,z) {}

    #if EIGEN_HAS_CXX11

    template <typename... ArgTypes>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    DiagonalMatrix(const Scalar& a0, const Scalar& a1, const Scalar& a2, const ArgTypes&... args)
      : m_diagonal(a0, a1, a2, args...) {}

    EIGEN_DEVICE_FUNC
    explicit EIGEN_STRONG_INLINE DiagonalMatrix(const std::initializer_list<std::initializer_list<Scalar>>& list)
      : m_diagonal(list) {}
    #endif  

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    inline DiagonalMatrix(const DiagonalBase<OtherDerived>& other) : m_diagonal(other.diagonal()) {}

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    
    inline DiagonalMatrix(const DiagonalMatrix& other) : m_diagonal(other.diagonal()) {}
    #endif

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    explicit inline DiagonalMatrix(const MatrixBase<OtherDerived>& other) : m_diagonal(other)
    {}

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    DiagonalMatrix& operator=(const DiagonalBase<OtherDerived>& other)
    {
      m_diagonal = other.diagonal();
      return *this;
    }

    #ifndef EIGEN_PARSED_BY_DOXYGEN

    EIGEN_DEVICE_FUNC
    DiagonalMatrix& operator=(const DiagonalMatrix& other)
    {
      m_diagonal = other.diagonal();
      return *this;
    }
    #endif

    EIGEN_DEVICE_FUNC
    inline void resize(Index size) { m_diagonal.resize(size); }
    
    EIGEN_DEVICE_FUNC
    inline void setZero() { m_diagonal.setZero(); }
    
    EIGEN_DEVICE_FUNC
    inline void setZero(Index size) { m_diagonal.setZero(size); }
    
    EIGEN_DEVICE_FUNC
    inline void setIdentity() { m_diagonal.setOnes(); }
    
    EIGEN_DEVICE_FUNC
    inline void setIdentity(Index size) { m_diagonal.setOnes(size); }
};

namespace internal {
template<typename DiagonalVectorType_>
struct traits<DiagonalWrapper<DiagonalVectorType_> >
{
  typedef DiagonalVectorType_ DiagonalVectorType;
  typedef typename DiagonalVectorType::Scalar Scalar;
  typedef typename DiagonalVectorType::StorageIndex StorageIndex;
  typedef DiagonalShape StorageKind;
  typedef typename traits<DiagonalVectorType>::XprKind XprKind;
  enum {
    RowsAtCompileTime = DiagonalVectorType::SizeAtCompileTime,
    ColsAtCompileTime = DiagonalVectorType::SizeAtCompileTime,
    MaxRowsAtCompileTime = DiagonalVectorType::MaxSizeAtCompileTime,
    MaxColsAtCompileTime = DiagonalVectorType::MaxSizeAtCompileTime,
    Flags =  (traits<DiagonalVectorType>::Flags & LvalueBit) | NoPreferredStorageOrderBit
  };
};
}

template<typename DiagonalVectorType_>
class DiagonalWrapper
  : public DiagonalBase<DiagonalWrapper<DiagonalVectorType_> >, internal::no_assignment_operator
{
  public:
    #ifndef EIGEN_PARSED_BY_DOXYGEN
    typedef DiagonalVectorType_ DiagonalVectorType;
    typedef DiagonalWrapper Nested;
    #endif

    EIGEN_DEVICE_FUNC
    explicit inline DiagonalWrapper(DiagonalVectorType& a_diagonal) : m_diagonal(a_diagonal) {}

    EIGEN_DEVICE_FUNC
    const DiagonalVectorType& diagonal() const { return m_diagonal; }

  protected:
    typename DiagonalVectorType::Nested m_diagonal;
};

template<typename Derived>
EIGEN_DEVICE_FUNC inline const DiagonalWrapper<const Derived>
MatrixBase<Derived>::asDiagonal() const
{
  return DiagonalWrapper<const Derived>(derived());
}

template<typename Derived>
bool MatrixBase<Derived>::isDiagonal(const RealScalar& prec) const
{
  if(cols() != rows()) return false;
  RealScalar maxAbsOnDiagonal = static_cast<RealScalar>(-1);
  for(Index j = 0; j < cols(); ++j)
  {
    RealScalar absOnDiagonal = numext::abs(coeff(j,j));
    if(absOnDiagonal > maxAbsOnDiagonal) maxAbsOnDiagonal = absOnDiagonal;
  }
  for(Index j = 0; j < cols(); ++j)
    for(Index i = 0; i < j; ++i)
    {
      if(!internal::isMuchSmallerThan(coeff(i, j), maxAbsOnDiagonal, prec)) return false;
      if(!internal::isMuchSmallerThan(coeff(j, i), maxAbsOnDiagonal, prec)) return false;
    }
  return true;
}

namespace internal {

template<> struct storage_kind_to_shape<DiagonalShape> { typedef DiagonalShape Shape; };

struct Diagonal2Dense {};

template<> struct AssignmentKind<DenseShape,DiagonalShape> { typedef Diagonal2Dense Kind; };

template< typename DstXprType, typename SrcXprType, typename Functor>
struct Assignment<DstXprType, SrcXprType, Functor, Diagonal2Dense>
{
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<typename DstXprType::Scalar,typename SrcXprType::Scalar> &)
  {
    Index dstRows = src.rows();
    Index dstCols = src.cols();
    if((dst.rows()!=dstRows) || (dst.cols()!=dstCols))
      dst.resize(dstRows, dstCols);
    
    dst.setZero();
    dst.diagonal() = src.diagonal();
  }
  
  static void run(DstXprType &dst, const SrcXprType &src, const internal::add_assign_op<typename DstXprType::Scalar,typename SrcXprType::Scalar> &)
  { dst.diagonal() += src.diagonal(); }
  
  static void run(DstXprType &dst, const SrcXprType &src, const internal::sub_assign_op<typename DstXprType::Scalar,typename SrcXprType::Scalar> &)
  { dst.diagonal() -= src.diagonal(); }
};

} 

} 

#endif 
