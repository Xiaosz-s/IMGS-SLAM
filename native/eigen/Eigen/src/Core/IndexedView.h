#ifndef EIGEN_INDEXED_VIEW_H
#define EIGEN_INDEXED_VIEW_H

namespace Eigen {

namespace internal {

template<typename XprType, typename RowIndices, typename ColIndices>
struct traits<IndexedView<XprType, RowIndices, ColIndices> >
 : traits<XprType>
{
  enum {
    RowsAtCompileTime = int(array_size<RowIndices>::value),
    ColsAtCompileTime = int(array_size<ColIndices>::value),
    MaxRowsAtCompileTime = RowsAtCompileTime != Dynamic ? int(RowsAtCompileTime) : Dynamic,
    MaxColsAtCompileTime = ColsAtCompileTime != Dynamic ? int(ColsAtCompileTime) : Dynamic,

    XprTypeIsRowMajor = (int(traits<XprType>::Flags)&RowMajorBit) != 0,
    IsRowMajor = (MaxRowsAtCompileTime==1&&MaxColsAtCompileTime!=1) ? 1
               : (MaxColsAtCompileTime==1&&MaxRowsAtCompileTime!=1) ? 0
               : XprTypeIsRowMajor,

    RowIncr = int(get_compile_time_incr<RowIndices>::value),
    ColIncr = int(get_compile_time_incr<ColIndices>::value),
    InnerIncr = IsRowMajor ? ColIncr : RowIncr,
    OuterIncr = IsRowMajor ? RowIncr : ColIncr,

    HasSameStorageOrderAsXprType = (IsRowMajor == XprTypeIsRowMajor),
    XprInnerStride = HasSameStorageOrderAsXprType ? int(inner_stride_at_compile_time<XprType>::ret) : int(outer_stride_at_compile_time<XprType>::ret),
    XprOuterstride = HasSameStorageOrderAsXprType ? int(outer_stride_at_compile_time<XprType>::ret) : int(inner_stride_at_compile_time<XprType>::ret),

    InnerSize = XprTypeIsRowMajor ? ColsAtCompileTime : RowsAtCompileTime,
    IsBlockAlike = InnerIncr==1 && OuterIncr==1,
    IsInnerPannel = HasSameStorageOrderAsXprType && is_same<AllRange<InnerSize>,typename conditional<XprTypeIsRowMajor,ColIndices,RowIndices>::type>::value,

    InnerStrideAtCompileTime = InnerIncr<0 || InnerIncr==DynamicIndex || XprInnerStride==Dynamic ? Dynamic : XprInnerStride * InnerIncr,
    OuterStrideAtCompileTime = OuterIncr<0 || OuterIncr==DynamicIndex || XprOuterstride==Dynamic ? Dynamic : XprOuterstride * OuterIncr,

    ReturnAsScalar = is_same<RowIndices,SingleRange>::value && is_same<ColIndices,SingleRange>::value,
    ReturnAsBlock = (!ReturnAsScalar) && IsBlockAlike,
    ReturnAsIndexedView = (!ReturnAsScalar) && (!ReturnAsBlock),

    DirectAccessMask = (int(InnerIncr)!=UndefinedIncr && int(OuterIncr)!=UndefinedIncr && InnerIncr>=0 && OuterIncr>=0) ? DirectAccessBit : 0,
    FlagsRowMajorBit = IsRowMajor ? RowMajorBit : 0,
    FlagsLvalueBit = is_lvalue<XprType>::value ? LvalueBit : 0,
    FlagsLinearAccessBit = (RowsAtCompileTime == 1 || ColsAtCompileTime == 1) ? LinearAccessBit : 0,
    Flags = (traits<XprType>::Flags & (HereditaryBits | DirectAccessMask )) | FlagsLvalueBit | FlagsRowMajorBit | FlagsLinearAccessBit
  };

  typedef Block<XprType,RowsAtCompileTime,ColsAtCompileTime,IsInnerPannel> BlockType;
};

}

template<typename XprType, typename RowIndices, typename ColIndices, typename StorageKind>
class IndexedViewImpl;

template<typename XprType, typename RowIndices, typename ColIndices>
class IndexedView : public IndexedViewImpl<XprType, RowIndices, ColIndices, typename internal::traits<XprType>::StorageKind>
{
public:
  typedef typename IndexedViewImpl<XprType, RowIndices, ColIndices, typename internal::traits<XprType>::StorageKind>::Base Base;
  EIGEN_GENERIC_PUBLIC_INTERFACE(IndexedView)
  EIGEN_INHERIT_ASSIGNMENT_OPERATORS(IndexedView)

  typedef typename internal::ref_selector<XprType>::non_const_type MatrixTypeNested;
  typedef typename internal::remove_all<XprType>::type NestedExpression;

  template<typename T0, typename T1>
  IndexedView(XprType& xpr, const T0& rowIndices, const T1& colIndices)
    : m_xpr(xpr), m_rowIndices(rowIndices), m_colIndices(colIndices)
  {}

  Index rows() const { return internal::size(m_rowIndices); }

  Index cols() const { return internal::size(m_colIndices); }

  const typename internal::remove_all<XprType>::type&
  nestedExpression() const { return m_xpr; }

  typename internal::remove_reference<XprType>::type&
  nestedExpression() { return m_xpr; }

  const RowIndices& rowIndices() const { return m_rowIndices; }

  const ColIndices& colIndices() const { return m_colIndices; }

protected:
  MatrixTypeNested m_xpr;
  RowIndices m_rowIndices;
  ColIndices m_colIndices;
};

template<typename XprType, typename RowIndices, typename ColIndices, typename StorageKind>
class IndexedViewImpl
  : public internal::generic_xpr_base<IndexedView<XprType, RowIndices, ColIndices> >::type
{
public:
  typedef typename internal::generic_xpr_base<IndexedView<XprType, RowIndices, ColIndices> >::type Base;
};

namespace internal {

template<typename ArgType, typename RowIndices, typename ColIndices>
struct unary_evaluator<IndexedView<ArgType, RowIndices, ColIndices>, IndexBased>
  : evaluator_base<IndexedView<ArgType, RowIndices, ColIndices> >
{
  typedef IndexedView<ArgType, RowIndices, ColIndices> XprType;

  enum {
    CoeffReadCost = evaluator<ArgType>::CoeffReadCost ,

    FlagsLinearAccessBit = (traits<XprType>::RowsAtCompileTime == 1 || traits<XprType>::ColsAtCompileTime == 1) ? LinearAccessBit : 0,

    FlagsRowMajorBit = traits<XprType>::FlagsRowMajorBit, 

    Flags = (evaluator<ArgType>::Flags & (HereditaryBits & ~RowMajorBit )) | FlagsLinearAccessBit | FlagsRowMajorBit,

    Alignment = 0
  };

  EIGEN_DEVICE_FUNC explicit unary_evaluator(const XprType& xpr) : m_argImpl(xpr.nestedExpression()), m_xpr(xpr)
  {
    EIGEN_INTERNAL_CHECK_COST_VALUE(CoeffReadCost);
  }

  typedef typename XprType::Scalar Scalar;
  typedef typename XprType::CoeffReturnType CoeffReturnType;

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  CoeffReturnType coeff(Index row, Index col) const
  {
    return m_argImpl.coeff(m_xpr.rowIndices()[row], m_xpr.colIndices()[col]);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  Scalar& coeffRef(Index row, Index col)
  {
    return m_argImpl.coeffRef(m_xpr.rowIndices()[row], m_xpr.colIndices()[col]);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  Scalar& coeffRef(Index index)
  {
    EIGEN_STATIC_ASSERT_LVALUE(XprType)
    Index row = XprType::RowsAtCompileTime == 1 ? 0 : index;
    Index col = XprType::RowsAtCompileTime == 1 ? index : 0;
    return m_argImpl.coeffRef( m_xpr.rowIndices()[row], m_xpr.colIndices()[col]);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  const Scalar& coeffRef(Index index) const
  {
    Index row = XprType::RowsAtCompileTime == 1 ? 0 : index;
    Index col = XprType::RowsAtCompileTime == 1 ? index : 0;
    return m_argImpl.coeffRef( m_xpr.rowIndices()[row], m_xpr.colIndices()[col]);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  const CoeffReturnType coeff(Index index) const
  {
    Index row = XprType::RowsAtCompileTime == 1 ? 0 : index;
    Index col = XprType::RowsAtCompileTime == 1 ? index : 0;
    return m_argImpl.coeff( m_xpr.rowIndices()[row], m_xpr.colIndices()[col]);
  }

protected:

  evaluator<ArgType> m_argImpl;
  const XprType& m_xpr;

};

} 

} 

#endif 
