#ifndef EIGEN_CWISE_UNARY_OP_H
#define EIGEN_CWISE_UNARY_OP_H

namespace Eigen {

namespace internal {
template<typename UnaryOp, typename XprType>
struct traits<CwiseUnaryOp<UnaryOp, XprType> >
 : traits<XprType>
{
  typedef typename result_of<
                     UnaryOp(const typename XprType::Scalar&)
                   >::type Scalar;
  typedef typename XprType::Nested XprTypeNested;
  typedef typename remove_reference<XprTypeNested>::type _XprTypeNested;
  enum {
    Flags = _XprTypeNested::Flags & RowMajorBit
  };
};
}

template<typename UnaryOp, typename XprType, typename StorageKind>
class CwiseUnaryOpImpl;

template<typename UnaryOp, typename XprType>
class CwiseUnaryOp : public CwiseUnaryOpImpl<UnaryOp, XprType, typename internal::traits<XprType>::StorageKind>, internal::no_assignment_operator
{
  public:

    typedef typename CwiseUnaryOpImpl<UnaryOp, XprType,typename internal::traits<XprType>::StorageKind>::Base Base;
    EIGEN_GENERIC_PUBLIC_INTERFACE(CwiseUnaryOp)
    typedef typename internal::ref_selector<XprType>::type XprTypeNested;
    typedef typename internal::remove_all<XprType>::type NestedExpression;

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    explicit CwiseUnaryOp(const XprType& xpr, const UnaryOp& func = UnaryOp())
      : m_xpr(xpr), m_functor(func) {}

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index rows() const EIGEN_NOEXCEPT { return m_xpr.rows(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index cols() const EIGEN_NOEXCEPT { return m_xpr.cols(); }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    const UnaryOp& functor() const { return m_functor; }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    const typename internal::remove_all<XprTypeNested>::type&
    nestedExpression() const { return m_xpr; }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    typename internal::remove_all<XprTypeNested>::type&
    nestedExpression() { return m_xpr; }

  protected:
    XprTypeNested m_xpr;
    const UnaryOp m_functor;
};

template<typename UnaryOp, typename XprType, typename StorageKind>
class CwiseUnaryOpImpl
  : public internal::generic_xpr_base<CwiseUnaryOp<UnaryOp, XprType> >::type
{
public:
  typedef typename internal::generic_xpr_base<CwiseUnaryOp<UnaryOp, XprType> >::type Base;
};

} 

#endif 
