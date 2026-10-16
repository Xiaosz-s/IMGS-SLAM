#ifndef EIGEN_REVERSE_H
#define EIGEN_REVERSE_H

namespace Eigen {

namespace internal {

template<typename MatrixType, int Direction>
struct traits<Reverse<MatrixType, Direction> >
 : traits<MatrixType>
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename traits<MatrixType>::StorageKind StorageKind;
  typedef typename traits<MatrixType>::XprKind XprKind;
  typedef typename ref_selector<MatrixType>::type MatrixTypeNested;
  typedef typename remove_reference<MatrixTypeNested>::type _MatrixTypeNested;
  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime,
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime,
    Flags = _MatrixTypeNested::Flags & (RowMajorBit | LvalueBit)
  };
};

template<typename PacketType, bool ReversePacket> struct reverse_packet_cond
{
  static inline PacketType run(const PacketType& x) { return preverse(x); }
};

template<typename PacketType> struct reverse_packet_cond<PacketType,false>
{
  static inline PacketType run(const PacketType& x) { return x; }
};

} 

template<typename MatrixType, int Direction> class Reverse
  : public internal::dense_xpr_base< Reverse<MatrixType, Direction> >::type
{
  public:

    typedef typename internal::dense_xpr_base<Reverse>::type Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(Reverse)
    typedef typename internal::remove_all<MatrixType>::type NestedExpression;
    using Base::IsRowMajor;

  protected:
    enum {
      PacketSize = internal::packet_traits<Scalar>::size,
      IsColMajor = !IsRowMajor,
      ReverseRow = (Direction == Vertical)   || (Direction == BothDirections),
      ReverseCol = (Direction == Horizontal) || (Direction == BothDirections),
      OffsetRow  = ReverseRow && IsColMajor ? PacketSize : 1,
      OffsetCol  = ReverseCol && IsRowMajor ? PacketSize : 1,
      ReversePacket = (Direction == BothDirections)
                    || ((Direction == Vertical)   && IsColMajor)
                    || ((Direction == Horizontal) && IsRowMajor)
    };
    typedef internal::reverse_packet_cond<PacketScalar,ReversePacket> reverse_packet;
  public:

    EIGEN_DEVICE_FUNC explicit inline Reverse(const MatrixType& matrix) : m_matrix(matrix) { }

    EIGEN_INHERIT_ASSIGNMENT_OPERATORS(Reverse)

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index rows() const EIGEN_NOEXCEPT { return m_matrix.rows(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index cols() const EIGEN_NOEXCEPT { return m_matrix.cols(); }

    EIGEN_DEVICE_FUNC inline Index innerStride() const
    {
      return -m_matrix.innerStride();
    }

    EIGEN_DEVICE_FUNC const typename internal::remove_all<typename MatrixType::Nested>::type&
    nestedExpression() const
    {
      return m_matrix;
    }

  protected:
    typename MatrixType::Nested m_matrix;
};

template<typename Derived>
EIGEN_DEVICE_FUNC inline typename DenseBase<Derived>::ReverseReturnType
DenseBase<Derived>::reverse()
{
  return ReverseReturnType(derived());
}

template<typename Derived>
EIGEN_DEVICE_FUNC inline void DenseBase<Derived>::reverseInPlace()
{
  if(cols()>rows())
  {
    Index half = cols()/2;
    leftCols(half).swap(rightCols(half).reverse());
    if((cols()%2)==1)
    {
      Index half2 = rows()/2;
      col(half).head(half2).swap(col(half).tail(half2).reverse());
    }
  }
  else
  {
    Index half = rows()/2;
    topRows(half).swap(bottomRows(half).reverse());
    if((rows()%2)==1)
    {
      Index half2 = cols()/2;
      row(half).head(half2).swap(row(half).tail(half2).reverse());
    }
  }
}

namespace internal {

template<int Direction>
struct vectorwise_reverse_inplace_impl;

template<>
struct vectorwise_reverse_inplace_impl<Vertical>
{
  template<typename ExpressionType>
  static void run(ExpressionType &xpr)
  {
    const int HalfAtCompileTime = ExpressionType::RowsAtCompileTime==Dynamic?Dynamic:ExpressionType::RowsAtCompileTime/2;
    Index half = xpr.rows()/2;
    xpr.topRows(fix<HalfAtCompileTime>(half))
       .swap(xpr.bottomRows(fix<HalfAtCompileTime>(half)).colwise().reverse());
  }
};

template<>
struct vectorwise_reverse_inplace_impl<Horizontal>
{
  template<typename ExpressionType>
  static void run(ExpressionType &xpr)
  {
    const int HalfAtCompileTime = ExpressionType::ColsAtCompileTime==Dynamic?Dynamic:ExpressionType::ColsAtCompileTime/2;
    Index half = xpr.cols()/2;
    xpr.leftCols(fix<HalfAtCompileTime>(half))
       .swap(xpr.rightCols(fix<HalfAtCompileTime>(half)).rowwise().reverse());
  }
};

} 

template<typename ExpressionType, int Direction>
EIGEN_DEVICE_FUNC void VectorwiseOp<ExpressionType,Direction>::reverseInPlace()
{
  internal::vectorwise_reverse_inplace_impl<Direction>::run(m_matrix);
}

} 

#endif 
