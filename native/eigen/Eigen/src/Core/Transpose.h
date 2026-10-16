#ifndef EIGEN_TRANSPOSE_H
#define EIGEN_TRANSPOSE_H

namespace Eigen {

namespace internal {
template<typename MatrixType>
struct traits<Transpose<MatrixType> > : public traits<MatrixType>
{
  typedef typename ref_selector<MatrixType>::type MatrixTypeNested;
  typedef typename remove_reference<MatrixTypeNested>::type MatrixTypeNestedPlain;
  enum {
    RowsAtCompileTime = MatrixType::ColsAtCompileTime,
    ColsAtCompileTime = MatrixType::RowsAtCompileTime,
    MaxRowsAtCompileTime = MatrixType::MaxColsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    FlagsLvalueBit = is_lvalue<MatrixType>::value ? LvalueBit : 0,
    Flags0 = traits<MatrixTypeNestedPlain>::Flags & ~(LvalueBit | NestByRefBit),
    Flags1 = Flags0 | FlagsLvalueBit,
    Flags = Flags1 ^ RowMajorBit,
    InnerStrideAtCompileTime = inner_stride_at_compile_time<MatrixType>::ret,
    OuterStrideAtCompileTime = outer_stride_at_compile_time<MatrixType>::ret
  };
};
}

template<typename MatrixType, typename StorageKind> class TransposeImpl;

template<typename MatrixType> class Transpose
  : public TransposeImpl<MatrixType,typename internal::traits<MatrixType>::StorageKind>
{
  public:

    typedef typename internal::ref_selector<MatrixType>::non_const_type MatrixTypeNested;

    typedef typename TransposeImpl<MatrixType,typename internal::traits<MatrixType>::StorageKind>::Base Base;
    EIGEN_GENERIC_PUBLIC_INTERFACE(Transpose)
    typedef typename internal::remove_all<MatrixType>::type NestedExpression;

    EIGEN_DEVICE_FUNC
    explicit EIGEN_STRONG_INLINE Transpose(MatrixType& matrix) : m_matrix(matrix) {}

    EIGEN_INHERIT_ASSIGNMENT_OPERATORS(Transpose)

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index rows() const EIGEN_NOEXCEPT { return m_matrix.cols(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index cols() const EIGEN_NOEXCEPT { return m_matrix.rows(); }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    const typename internal::remove_all<MatrixTypeNested>::type&
    nestedExpression() const { return m_matrix; }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    typename internal::remove_reference<MatrixTypeNested>::type&
    nestedExpression() { return m_matrix; }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    void resize(Index nrows, Index ncols) {
      m_matrix.resize(ncols,nrows);
    }

  protected:
    typename internal::ref_selector<MatrixType>::non_const_type m_matrix;
};

namespace internal {

template<typename MatrixType, bool HasDirectAccess = has_direct_access<MatrixType>::ret>
struct TransposeImpl_base
{
  typedef typename dense_xpr_base<Transpose<MatrixType> >::type type;
};

template<typename MatrixType>
struct TransposeImpl_base<MatrixType, false>
{
  typedef typename dense_xpr_base<Transpose<MatrixType> >::type type;
};

} 

template<typename XprType, typename StorageKind>
class TransposeImpl
  : public internal::generic_xpr_base<Transpose<XprType> >::type
{
public:
  typedef typename internal::generic_xpr_base<Transpose<XprType> >::type Base;
};

template<typename MatrixType> class TransposeImpl<MatrixType,Dense>
  : public internal::TransposeImpl_base<MatrixType>::type
{
  public:

    typedef typename internal::TransposeImpl_base<MatrixType>::type Base;
    using Base::coeffRef;
    EIGEN_DENSE_PUBLIC_INTERFACE(Transpose<MatrixType>)
    EIGEN_INHERIT_ASSIGNMENT_OPERATORS(TransposeImpl)

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Index innerStride() const { return derived().nestedExpression().innerStride(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Index outerStride() const { return derived().nestedExpression().outerStride(); }

    typedef typename internal::conditional<
                       internal::is_lvalue<MatrixType>::value,
                       Scalar,
                       const Scalar
                     >::type ScalarWithConstIfNotLvalue;

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    ScalarWithConstIfNotLvalue* data() { return derived().nestedExpression().data(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    const Scalar* data() const { return derived().nestedExpression().data(); }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    const Scalar& coeffRef(Index rowId, Index colId) const
    {
      return derived().nestedExpression().coeffRef(colId, rowId);
    }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    const Scalar& coeffRef(Index index) const
    {
      return derived().nestedExpression().coeffRef(index);
    }
  protected:
    EIGEN_DEFAULT_EMPTY_CONSTRUCTOR_AND_DESTRUCTOR(TransposeImpl)
};

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
Transpose<Derived>
DenseBase<Derived>::transpose()
{
  return TransposeReturnType(derived());
}

template<typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename DenseBase<Derived>::ConstTransposeReturnType
DenseBase<Derived>::transpose() const
{
  return ConstTransposeReturnType(derived());
}

template<typename Derived>
EIGEN_DEVICE_FUNC inline const typename MatrixBase<Derived>::AdjointReturnType
MatrixBase<Derived>::adjoint() const
{
  return AdjointReturnType(this->transpose());
}

namespace internal {

template<typename MatrixType,
  bool IsSquare = (MatrixType::RowsAtCompileTime == MatrixType::ColsAtCompileTime) && MatrixType::RowsAtCompileTime!=Dynamic,
  bool MatchPacketSize =
        (int(MatrixType::RowsAtCompileTime) == int(internal::packet_traits<typename MatrixType::Scalar>::size))
    &&  (internal::evaluator<MatrixType>::Flags&PacketAccessBit) >
struct inplace_transpose_selector;

template<typename MatrixType>
struct inplace_transpose_selector<MatrixType,true,false> { 
  static void run(MatrixType& m) {
    m.matrix().template triangularView<StrictlyUpper>().swap(m.matrix().transpose().template triangularView<StrictlyUpper>());
  }
};

template<typename MatrixType>
struct inplace_transpose_selector<MatrixType,true,true> { 
  static void run(MatrixType& m) {
    typedef typename MatrixType::Scalar Scalar;
    typedef typename internal::packet_traits<typename MatrixType::Scalar>::type Packet;
    const Index PacketSize = internal::packet_traits<Scalar>::size;
    const Index Alignment = internal::evaluator<MatrixType>::Alignment;
    PacketBlock<Packet> A;
    for (Index i=0; i<PacketSize; ++i)
      A.packet[i] = m.template packetByOuterInner<Alignment>(i,0);
    internal::ptranspose(A);
    for (Index i=0; i<PacketSize; ++i)
      m.template writePacket<Alignment>(m.rowIndexByOuterInner(i,0), m.colIndexByOuterInner(i,0), A.packet[i]);
  }
};

template <typename MatrixType, Index Alignment>
void BlockedInPlaceTranspose(MatrixType& m) {
  typedef typename MatrixType::Scalar Scalar;
  typedef typename internal::packet_traits<typename MatrixType::Scalar>::type Packet;
  const Index PacketSize = internal::packet_traits<Scalar>::size;
  eigen_assert(m.rows() == m.cols());
  int row_start = 0;
  for (; row_start + PacketSize <= m.rows(); row_start += PacketSize) {
    for (int col_start = row_start; col_start + PacketSize <= m.cols(); col_start += PacketSize) {
      PacketBlock<Packet> A;
      if (row_start == col_start) {
        for (Index i=0; i<PacketSize; ++i)
          A.packet[i] = m.template packetByOuterInner<Alignment>(row_start + i,col_start);
        internal::ptranspose(A);
        for (Index i=0; i<PacketSize; ++i)
          m.template writePacket<Alignment>(m.rowIndexByOuterInner(row_start + i, col_start), m.colIndexByOuterInner(row_start + i,col_start), A.packet[i]);
      } else {
        PacketBlock<Packet> B;
        for (Index i=0; i<PacketSize; ++i) {
          A.packet[i] = m.template packetByOuterInner<Alignment>(row_start + i,col_start);
          B.packet[i] = m.template packetByOuterInner<Alignment>(col_start + i, row_start);
        }
        internal::ptranspose(A);
        internal::ptranspose(B);
        for (Index i=0; i<PacketSize; ++i) {
          m.template writePacket<Alignment>(m.rowIndexByOuterInner(row_start + i, col_start), m.colIndexByOuterInner(row_start + i,col_start), B.packet[i]);
          m.template writePacket<Alignment>(m.rowIndexByOuterInner(col_start + i, row_start), m.colIndexByOuterInner(col_start + i,row_start), A.packet[i]);
        }
      }
    }
  }
  for (Index row = row_start; row < m.rows(); ++row) {
    m.matrix().row(row).head(row).swap(
        m.matrix().col(row).head(row).transpose());
  }
}

template<typename MatrixType,bool MatchPacketSize>
struct inplace_transpose_selector<MatrixType,false,MatchPacketSize> { 
  static void run(MatrixType& m) {
    typedef typename MatrixType::Scalar Scalar;
    if (m.rows() == m.cols()) {
      const Index PacketSize = internal::packet_traits<Scalar>::size;
      if (!NumTraits<Scalar>::IsComplex && m.rows() >= PacketSize) {
        if ((m.rows() % PacketSize) == 0)
          BlockedInPlaceTranspose<MatrixType,internal::evaluator<MatrixType>::Alignment>(m);
        else
          BlockedInPlaceTranspose<MatrixType,Unaligned>(m);
      }
      else {
        m.matrix().template triangularView<StrictlyUpper>().swap(m.matrix().transpose().template triangularView<StrictlyUpper>());
      }
    } else {
      m = m.transpose().eval();
    }
  }
};

} 

template<typename Derived>
EIGEN_DEVICE_FUNC inline void DenseBase<Derived>::transposeInPlace()
{
  eigen_assert((rows() == cols() || (RowsAtCompileTime == Dynamic && ColsAtCompileTime == Dynamic))
               && "transposeInPlace() called on a non-square non-resizable matrix");
  internal::inplace_transpose_selector<Derived>::run(derived());
}

template<typename Derived>
EIGEN_DEVICE_FUNC inline void MatrixBase<Derived>::adjointInPlace()
{
  derived() = adjoint().eval();
}

#ifndef EIGEN_NO_DEBUG

namespace internal {

template<bool DestIsTransposed, typename OtherDerived>
struct check_transpose_aliasing_compile_time_selector
{
  enum { ret = bool(blas_traits<OtherDerived>::IsTransposed) != DestIsTransposed };
};

template<bool DestIsTransposed, typename BinOp, typename DerivedA, typename DerivedB>
struct check_transpose_aliasing_compile_time_selector<DestIsTransposed,CwiseBinaryOp<BinOp,DerivedA,DerivedB> >
{
  enum { ret =    bool(blas_traits<DerivedA>::IsTransposed) != DestIsTransposed
               || bool(blas_traits<DerivedB>::IsTransposed) != DestIsTransposed
  };
};

template<typename Scalar, bool DestIsTransposed, typename OtherDerived>
struct check_transpose_aliasing_run_time_selector
{
  static bool run(const Scalar* dest, const OtherDerived& src)
  {
    return (bool(blas_traits<OtherDerived>::IsTransposed) != DestIsTransposed) && (dest!=0 && dest==(const Scalar*)extract_data(src));
  }
};

template<typename Scalar, bool DestIsTransposed, typename BinOp, typename DerivedA, typename DerivedB>
struct check_transpose_aliasing_run_time_selector<Scalar,DestIsTransposed,CwiseBinaryOp<BinOp,DerivedA,DerivedB> >
{
  static bool run(const Scalar* dest, const CwiseBinaryOp<BinOp,DerivedA,DerivedB>& src)
  {
    return ((blas_traits<DerivedA>::IsTransposed != DestIsTransposed) && (dest!=0 && dest==(const Scalar*)extract_data(src.lhs())))
        || ((blas_traits<DerivedB>::IsTransposed != DestIsTransposed) && (dest!=0 && dest==(const Scalar*)extract_data(src.rhs())));
  }
};

template<typename Derived, typename OtherDerived,
         bool MightHaveTransposeAliasing
                 = check_transpose_aliasing_compile_time_selector
                     <blas_traits<Derived>::IsTransposed,OtherDerived>::ret
        >
struct checkTransposeAliasing_impl
{
    static void run(const Derived& dst, const OtherDerived& other)
    {
        eigen_assert((!check_transpose_aliasing_run_time_selector
                      <typename Derived::Scalar,blas_traits<Derived>::IsTransposed,OtherDerived>
                      ::run(extract_data(dst), other))
          && "aliasing detected during transposition, use transposeInPlace() "
             "or evaluate the rhs into a temporary using .eval()");

    }
};

template<typename Derived, typename OtherDerived>
struct checkTransposeAliasing_impl<Derived, OtherDerived, false>
{
    static void run(const Derived&, const OtherDerived&)
    {
    }
};

template<typename Dst, typename Src>
void check_for_aliasing(const Dst &dst, const Src &src)
{
  if((!Dst::IsVectorAtCompileTime) && dst.rows()>1 && dst.cols()>1)
    internal::checkTransposeAliasing_impl<Dst, Src>::run(dst, src);
}

} 

#endif 

} 

#endif 
