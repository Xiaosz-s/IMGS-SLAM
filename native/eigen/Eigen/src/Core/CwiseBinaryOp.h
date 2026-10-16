#ifndef EIGEN_CWISE_BINARY_OP_H
#define EIGEN_CWISE_BINARY_OP_H

namespace Eigen {

namespace internal {
template<typename BinaryOp, typename Lhs, typename Rhs>
struct traits<CwiseBinaryOp<BinaryOp, Lhs, Rhs> >
{

  typedef typename remove_all<Lhs>::type Ancestor;
  typedef typename traits<Ancestor>::XprKind XprKind;
  enum {
    RowsAtCompileTime = traits<Ancestor>::RowsAtCompileTime,
    ColsAtCompileTime = traits<Ancestor>::ColsAtCompileTime,
    MaxRowsAtCompileTime = traits<Ancestor>::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = traits<Ancestor>::MaxColsAtCompileTime
  };

  typedef typename result_of<
                     BinaryOp(
                       const typename Lhs::Scalar&,
                       const typename Rhs::Scalar&
                     )
                   >::type Scalar;
  typedef typename cwise_promote_storage_type<typename traits<Lhs>::StorageKind,
                                              typename traits<Rhs>::StorageKind,
                                              BinaryOp>::ret StorageKind;
  typedef typename promote_index_type<typename traits<Lhs>::StorageIndex,
                                      typename traits<Rhs>::StorageIndex>::type StorageIndex;
  typedef typename Lhs::Nested LhsNested;
  typedef typename Rhs::Nested RhsNested;
  typedef typename remove_reference<LhsNested>::type _LhsNested;
  typedef typename remove_reference<RhsNested>::type _RhsNested;
  enum {
    Flags = cwise_promote_storage_order<typename traits<Lhs>::StorageKind,typename traits<Rhs>::StorageKind,_LhsNested::Flags & RowMajorBit,_RhsNested::Flags & RowMajorBit>::value
  };
};
} 

template<typename BinaryOp, typename Lhs, typename Rhs, typename StorageKind>
class CwiseBinaryOpImpl;

template<typename BinaryOp, typename LhsType, typename RhsType>
class CwiseBinaryOp :
  public CwiseBinaryOpImpl<
          BinaryOp, LhsType, RhsType,
          typename internal::cwise_promote_storage_type<typename internal::traits<LhsType>::StorageKind,
                                                        typename internal::traits<RhsType>::StorageKind,
                                                        BinaryOp>::ret>,
  internal::no_assignment_operator
{
  public:

    typedef typename internal::remove_all<BinaryOp>::type Functor;
    typedef typename internal::remove_all<LhsType>::type Lhs;
    typedef typename internal::remove_all<RhsType>::type Rhs;

    typedef typename CwiseBinaryOpImpl<
        BinaryOp, LhsType, RhsType,
        typename internal::cwise_promote_storage_type<typename internal::traits<LhsType>::StorageKind,
                                                      typename internal::traits<Rhs>::StorageKind,
                                                      BinaryOp>::ret>::Base Base;
    EIGEN_GENERIC_PUBLIC_INTERFACE(CwiseBinaryOp)

    typedef typename internal::ref_selector<LhsType>::type LhsNested;
    typedef typename internal::ref_selector<RhsType>::type RhsNested;
    typedef typename internal::remove_reference<LhsNested>::type _LhsNested;
    typedef typename internal::remove_reference<RhsNested>::type _RhsNested;

#if EIGEN_COMP_MSVC && EIGEN_HAS_CXX11
    
    EIGEN_STRONG_INLINE
    CwiseBinaryOp(const CwiseBinaryOp<BinaryOp,LhsType,RhsType>&) = default;
#endif

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    CwiseBinaryOp(const Lhs& aLhs, const Rhs& aRhs, const BinaryOp& func = BinaryOp())
      : m_lhs(aLhs), m_rhs(aRhs), m_functor(func)
    {
      EIGEN_CHECK_BINARY_COMPATIBILIY(BinaryOp,typename Lhs::Scalar,typename Rhs::Scalar);
      
      EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Lhs, Rhs)
      eigen_assert(aLhs.rows() == aRhs.rows() && aLhs.cols() == aRhs.cols());
    }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index rows() const EIGEN_NOEXCEPT {
      
      return internal::traits<typename internal::remove_all<LhsNested>::type>::RowsAtCompileTime==Dynamic ? m_rhs.rows() : m_lhs.rows();
    }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index cols() const EIGEN_NOEXCEPT {
      
      return internal::traits<typename internal::remove_all<LhsNested>::type>::ColsAtCompileTime==Dynamic ? m_rhs.cols() : m_lhs.cols();
    }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    const _LhsNested& lhs() const { return m_lhs; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    const _RhsNested& rhs() const { return m_rhs; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    const BinaryOp& functor() const { return m_functor; }

  protected:
    LhsNested m_lhs;
    RhsNested m_rhs;
    const BinaryOp m_functor;
};

template<typename BinaryOp, typename Lhs, typename Rhs, typename StorageKind>
class CwiseBinaryOpImpl
  : public internal::generic_xpr_base<CwiseBinaryOp<BinaryOp, Lhs, Rhs> >::type
{
public:
  typedef typename internal::generic_xpr_base<CwiseBinaryOp<BinaryOp, Lhs, Rhs> >::type Base;
};

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived &
MatrixBase<Derived>::operator-=(const MatrixBase<OtherDerived> &other)
{
  call_assignment(derived(), other.derived(), internal::sub_assign_op<Scalar,typename OtherDerived::Scalar>());
  return derived();
}

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Derived &
MatrixBase<Derived>::operator+=(const MatrixBase<OtherDerived>& other)
{
  call_assignment(derived(), other.derived(), internal::add_assign_op<Scalar,typename OtherDerived::Scalar>());
  return derived();
}

} 

#endif 
