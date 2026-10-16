#ifndef EIGEN_DIAGONAL_H
#define EIGEN_DIAGONAL_H

namespace Eigen {

namespace internal {
template<typename MatrixType, int DiagIndex>
struct traits<Diagonal<MatrixType,DiagIndex> >
 : traits<MatrixType>
{
  typedef typename ref_selector<MatrixType>::type MatrixTypeNested;
  typedef typename remove_reference<MatrixTypeNested>::type _MatrixTypeNested;
  typedef typename MatrixType::StorageKind StorageKind;
  enum {
    RowsAtCompileTime = (int(DiagIndex) == DynamicIndex || int(MatrixType::SizeAtCompileTime) == Dynamic) ? Dynamic
                      : (EIGEN_PLAIN_ENUM_MIN(MatrixType::RowsAtCompileTime - EIGEN_PLAIN_ENUM_MAX(-DiagIndex, 0),
                                              MatrixType::ColsAtCompileTime - EIGEN_PLAIN_ENUM_MAX( DiagIndex, 0))),
    ColsAtCompileTime = 1,
    MaxRowsAtCompileTime = int(MatrixType::MaxSizeAtCompileTime) == Dynamic ? Dynamic
                         : DiagIndex == DynamicIndex ? EIGEN_SIZE_MIN_PREFER_FIXED(MatrixType::MaxRowsAtCompileTime,
                                                                              MatrixType::MaxColsAtCompileTime)
                         : (EIGEN_PLAIN_ENUM_MIN(MatrixType::MaxRowsAtCompileTime - EIGEN_PLAIN_ENUM_MAX(-DiagIndex, 0),
                                                 MatrixType::MaxColsAtCompileTime - EIGEN_PLAIN_ENUM_MAX( DiagIndex, 0))),
    MaxColsAtCompileTime = 1,
    MaskLvalueBit = is_lvalue<MatrixType>::value ? LvalueBit : 0,
    Flags = (unsigned int)_MatrixTypeNested::Flags & (RowMajorBit | MaskLvalueBit | DirectAccessBit) & ~RowMajorBit, 
    MatrixTypeOuterStride = outer_stride_at_compile_time<MatrixType>::ret,
    InnerStrideAtCompileTime = MatrixTypeOuterStride == Dynamic ? Dynamic : MatrixTypeOuterStride+1,
    OuterStrideAtCompileTime = 0
  };
};
}

template<typename MatrixType, int DiagIndex_> class Diagonal
   : public internal::dense_xpr_base< Diagonal<MatrixType,DiagIndex_> >::type
{
  public:

    enum { DiagIndex = DiagIndex_ };
    typedef typename internal::dense_xpr_base<Diagonal>::type Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(Diagonal)

    EIGEN_DEVICE_FUNC
    explicit inline Diagonal(MatrixType& matrix, Index a_index = DiagIndex) : m_matrix(matrix), m_index(a_index)
    {
      eigen_assert( a_index <= m_matrix.cols() && -a_index <= m_matrix.rows() );
    }

    EIGEN_INHERIT_ASSIGNMENT_OPERATORS(Diagonal)

    EIGEN_DEVICE_FUNC
    inline Index rows() const
    {
      return m_index.value()<0 ? numext::mini<Index>(m_matrix.cols(),m_matrix.rows()+m_index.value())
                               : numext::mini<Index>(m_matrix.rows(),m_matrix.cols()-m_index.value());
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index cols() const EIGEN_NOEXCEPT { return 1; }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index innerStride() const EIGEN_NOEXCEPT {
      return m_matrix.outerStride() + 1;
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index outerStride() const EIGEN_NOEXCEPT { return 0; }

    typedef typename internal::conditional<
                       internal::is_lvalue<MatrixType>::value,
                       Scalar,
                       const Scalar
                     >::type ScalarWithConstIfNotLvalue;

    EIGEN_DEVICE_FUNC
    inline ScalarWithConstIfNotLvalue* data() { return &(m_matrix.coeffRef(rowOffset(), colOffset())); }
    EIGEN_DEVICE_FUNC
    inline const Scalar* data() const { return &(m_matrix.coeffRef(rowOffset(), colOffset())); }

    EIGEN_DEVICE_FUNC
    inline Scalar& coeffRef(Index row, Index)
    {
      EIGEN_STATIC_ASSERT_LVALUE(MatrixType)
      return m_matrix.coeffRef(row+rowOffset(), row+colOffset());
    }

    EIGEN_DEVICE_FUNC
    inline const Scalar& coeffRef(Index row, Index) const
    {
      return m_matrix.coeffRef(row+rowOffset(), row+colOffset());
    }

    EIGEN_DEVICE_FUNC
    inline CoeffReturnType coeff(Index row, Index) const
    {
      return m_matrix.coeff(row+rowOffset(), row+colOffset());
    }

    EIGEN_DEVICE_FUNC
    inline Scalar& coeffRef(Index idx)
    {
      EIGEN_STATIC_ASSERT_LVALUE(MatrixType)
      return m_matrix.coeffRef(idx+rowOffset(), idx+colOffset());
    }

    EIGEN_DEVICE_FUNC
    inline const Scalar& coeffRef(Index idx) const
    {
      return m_matrix.coeffRef(idx+rowOffset(), idx+colOffset());
    }

    EIGEN_DEVICE_FUNC
    inline CoeffReturnType coeff(Index idx) const
    {
      return m_matrix.coeff(idx+rowOffset(), idx+colOffset());
    }

    EIGEN_DEVICE_FUNC
    inline const typename internal::remove_all<typename MatrixType::Nested>::type&
    nestedExpression() const
    {
      return m_matrix;
    }

    EIGEN_DEVICE_FUNC
    inline Index index() const
    {
      return m_index.value();
    }

  protected:
    typename internal::ref_selector<MatrixType>::non_const_type m_matrix;
    const internal::variable_if_dynamicindex<Index, DiagIndex> m_index;

  private:
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index absDiagIndex() const EIGEN_NOEXCEPT { return m_index.value()>0 ? m_index.value() : -m_index.value(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index rowOffset() const EIGEN_NOEXCEPT { return m_index.value()>0 ? 0 : -m_index.value(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index colOffset() const EIGEN_NOEXCEPT { return m_index.value()>0 ? m_index.value() : 0; }
    
    template<int LoadMode> typename MatrixType::PacketReturnType packet(Index) const;
    template<int LoadMode> typename MatrixType::PacketReturnType packet(Index,Index) const;
};

template<typename Derived>
EIGEN_DEVICE_FUNC inline typename MatrixBase<Derived>::DiagonalReturnType
MatrixBase<Derived>::diagonal()
{
  return DiagonalReturnType(derived());
}

template<typename Derived>
EIGEN_DEVICE_FUNC inline typename MatrixBase<Derived>::ConstDiagonalReturnType
MatrixBase<Derived>::diagonal() const
{
  return ConstDiagonalReturnType(derived());
}

template<typename Derived>
EIGEN_DEVICE_FUNC inline typename MatrixBase<Derived>::DiagonalDynamicIndexReturnType
MatrixBase<Derived>::diagonal(Index index)
{
  return DiagonalDynamicIndexReturnType(derived(), index);
}

template<typename Derived>
EIGEN_DEVICE_FUNC inline typename MatrixBase<Derived>::ConstDiagonalDynamicIndexReturnType
MatrixBase<Derived>::diagonal(Index index) const
{
  return ConstDiagonalDynamicIndexReturnType(derived(), index);
}

template<typename Derived>
template<int Index_>
EIGEN_DEVICE_FUNC
inline typename MatrixBase<Derived>::template DiagonalIndexReturnType<Index_>::Type
MatrixBase<Derived>::diagonal()
{
  return typename DiagonalIndexReturnType<Index_>::Type(derived());
}

template<typename Derived>
template<int Index_>
EIGEN_DEVICE_FUNC
inline typename MatrixBase<Derived>::template ConstDiagonalIndexReturnType<Index_>::Type
MatrixBase<Derived>::diagonal() const
{
  return typename ConstDiagonalIndexReturnType<Index_>::Type(derived());
}

} 

#endif 
