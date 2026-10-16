#ifndef EIGEN_DENSECOEFFSBASE_H
#define EIGEN_DENSECOEFFSBASE_H

namespace Eigen {

namespace internal {
template<typename T> struct add_const_on_value_type_if_arithmetic
{
  typedef typename conditional<is_arithmetic<T>::value, T, typename add_const_on_value_type<T>::type>::type type;
};
}

template<typename Derived>
class DenseCoeffsBase<Derived,ReadOnlyAccessors> : public EigenBase<Derived>
{
  public:
    typedef typename internal::traits<Derived>::StorageKind StorageKind;
    typedef typename internal::traits<Derived>::Scalar Scalar;
    typedef typename internal::packet_traits<Scalar>::type PacketScalar;

    typedef typename internal::conditional<bool(internal::traits<Derived>::Flags&LvalueBit),
                         const Scalar&,
                         typename internal::conditional<internal::is_arithmetic<Scalar>::value, Scalar, const Scalar>::type
                     >::type CoeffReturnType;

    typedef typename internal::add_const_on_value_type_if_arithmetic<
                         typename internal::packet_traits<Scalar>::type
                     >::type PacketReturnType;

    typedef EigenBase<Derived> Base;
    using Base::rows;
    using Base::cols;
    using Base::size;
    using Base::derived;

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Index rowIndexByOuterInner(Index outer, Index inner) const
    {
      return int(Derived::RowsAtCompileTime) == 1 ? 0
          : int(Derived::ColsAtCompileTime) == 1 ? inner
          : int(Derived::Flags)&RowMajorBit ? outer
          : inner;
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Index colIndexByOuterInner(Index outer, Index inner) const
    {
      return int(Derived::ColsAtCompileTime) == 1 ? 0
          : int(Derived::RowsAtCompileTime) == 1 ? inner
          : int(Derived::Flags)&RowMajorBit ? inner
          : outer;
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE CoeffReturnType coeff(Index row, Index col) const
    {
      eigen_internal_assert(row >= 0 && row < rows()
                         && col >= 0 && col < cols());
      return internal::evaluator<Derived>(derived()).coeff(row,col);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE CoeffReturnType coeffByOuterInner(Index outer, Index inner) const
    {
      return coeff(rowIndexByOuterInner(outer, inner),
                   colIndexByOuterInner(outer, inner));
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE CoeffReturnType operator()(Index row, Index col) const
    {
      eigen_assert(row >= 0 && row < rows()
          && col >= 0 && col < cols());
      return coeff(row, col);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE CoeffReturnType
    coeff(Index index) const
    {
      EIGEN_STATIC_ASSERT(internal::evaluator<Derived>::Flags & LinearAccessBit,
                          THIS_COEFFICIENT_ACCESSOR_TAKING_ONE_ACCESS_IS_ONLY_FOR_EXPRESSIONS_ALLOWING_LINEAR_ACCESS)
      eigen_internal_assert(index >= 0 && index < size());
      return internal::evaluator<Derived>(derived()).coeff(index);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE CoeffReturnType
    operator[](Index index) const
    {
      EIGEN_STATIC_ASSERT(Derived::IsVectorAtCompileTime,
                          THE_BRACKET_OPERATOR_IS_ONLY_FOR_VECTORS__USE_THE_PARENTHESIS_OPERATOR_INSTEAD)
      eigen_assert(index >= 0 && index < size());
      return coeff(index);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE CoeffReturnType
    operator()(Index index) const
    {
      eigen_assert(index >= 0 && index < size());
      return coeff(index);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE CoeffReturnType
    x() const { return (*this)[0]; }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE CoeffReturnType
    y() const
    {
      EIGEN_STATIC_ASSERT(Derived::SizeAtCompileTime==-1 || Derived::SizeAtCompileTime>=2, OUT_OF_RANGE_ACCESS);
      return (*this)[1];
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE CoeffReturnType
    z() const
    {
      EIGEN_STATIC_ASSERT(Derived::SizeAtCompileTime==-1 || Derived::SizeAtCompileTime>=3, OUT_OF_RANGE_ACCESS);
      return (*this)[2];
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE CoeffReturnType
    w() const
    {
      EIGEN_STATIC_ASSERT(Derived::SizeAtCompileTime==-1 || Derived::SizeAtCompileTime>=4, OUT_OF_RANGE_ACCESS);
      return (*this)[3];
    }

    template<int LoadMode>
    EIGEN_STRONG_INLINE PacketReturnType packet(Index row, Index col) const
    {
      typedef typename internal::packet_traits<Scalar>::type DefaultPacketType;
      eigen_internal_assert(row >= 0 && row < rows() && col >= 0 && col < cols());
      return internal::evaluator<Derived>(derived()).template packet<LoadMode,DefaultPacketType>(row,col);
    }

    template<int LoadMode>
    EIGEN_STRONG_INLINE PacketReturnType packetByOuterInner(Index outer, Index inner) const
    {
      return packet<LoadMode>(rowIndexByOuterInner(outer, inner),
                              colIndexByOuterInner(outer, inner));
    }

    template<int LoadMode>
    EIGEN_STRONG_INLINE PacketReturnType packet(Index index) const
    {
      EIGEN_STATIC_ASSERT(internal::evaluator<Derived>::Flags & LinearAccessBit,
                          THIS_COEFFICIENT_ACCESSOR_TAKING_ONE_ACCESS_IS_ONLY_FOR_EXPRESSIONS_ALLOWING_LINEAR_ACCESS)
      typedef typename internal::packet_traits<Scalar>::type DefaultPacketType;
      eigen_internal_assert(index >= 0 && index < size());
      return internal::evaluator<Derived>(derived()).template packet<LoadMode,DefaultPacketType>(index);
    }

  protected:

    void coeffRef();
    void coeffRefByOuterInner();
    void writePacket();
    void writePacketByOuterInner();
    void copyCoeff();
    void copyCoeffByOuterInner();
    void copyPacket();
    void copyPacketByOuterInner();
    void stride();
    void innerStride();
    void outerStride();
    void rowStride();
    void colStride();
};

template<typename Derived>
class DenseCoeffsBase<Derived, WriteAccessors> : public DenseCoeffsBase<Derived, ReadOnlyAccessors>
{
  public:

    typedef DenseCoeffsBase<Derived, ReadOnlyAccessors> Base;

    typedef typename internal::traits<Derived>::StorageKind StorageKind;
    typedef typename internal::traits<Derived>::Scalar Scalar;
    typedef typename internal::packet_traits<Scalar>::type PacketScalar;
    typedef typename NumTraits<Scalar>::Real RealScalar;

    using Base::coeff;
    using Base::rows;
    using Base::cols;
    using Base::size;
    using Base::derived;
    using Base::rowIndexByOuterInner;
    using Base::colIndexByOuterInner;
    using Base::operator[];
    using Base::operator();
    using Base::x;
    using Base::y;
    using Base::z;
    using Base::w;

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar& coeffRef(Index row, Index col)
    {
      eigen_internal_assert(row >= 0 && row < rows()
                         && col >= 0 && col < cols());
      return internal::evaluator<Derived>(derived()).coeffRef(row,col);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar&
    coeffRefByOuterInner(Index outer, Index inner)
    {
      return coeffRef(rowIndexByOuterInner(outer, inner),
                      colIndexByOuterInner(outer, inner));
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar&
    operator()(Index row, Index col)
    {
      eigen_assert(row >= 0 && row < rows()
          && col >= 0 && col < cols());
      return coeffRef(row, col);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar&
    coeffRef(Index index)
    {
      EIGEN_STATIC_ASSERT(internal::evaluator<Derived>::Flags & LinearAccessBit,
                          THIS_COEFFICIENT_ACCESSOR_TAKING_ONE_ACCESS_IS_ONLY_FOR_EXPRESSIONS_ALLOWING_LINEAR_ACCESS)
      eigen_internal_assert(index >= 0 && index < size());
      return internal::evaluator<Derived>(derived()).coeffRef(index);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar&
    operator[](Index index)
    {
      EIGEN_STATIC_ASSERT(Derived::IsVectorAtCompileTime,
                          THE_BRACKET_OPERATOR_IS_ONLY_FOR_VECTORS__USE_THE_PARENTHESIS_OPERATOR_INSTEAD)
      eigen_assert(index >= 0 && index < size());
      return coeffRef(index);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar&
    operator()(Index index)
    {
      eigen_assert(index >= 0 && index < size());
      return coeffRef(index);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar&
    x() { return (*this)[0]; }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar&
    y()
    {
      EIGEN_STATIC_ASSERT(Derived::SizeAtCompileTime==-1 || Derived::SizeAtCompileTime>=2, OUT_OF_RANGE_ACCESS);
      return (*this)[1];
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar&
    z()
    {
      EIGEN_STATIC_ASSERT(Derived::SizeAtCompileTime==-1 || Derived::SizeAtCompileTime>=3, OUT_OF_RANGE_ACCESS);
      return (*this)[2];
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar&
    w()
    {
      EIGEN_STATIC_ASSERT(Derived::SizeAtCompileTime==-1 || Derived::SizeAtCompileTime>=4, OUT_OF_RANGE_ACCESS);
      return (*this)[3];
    }
};

template<typename Derived>
class DenseCoeffsBase<Derived, DirectAccessors> : public DenseCoeffsBase<Derived, ReadOnlyAccessors>
{
  public:

    typedef DenseCoeffsBase<Derived, ReadOnlyAccessors> Base;
    typedef typename internal::traits<Derived>::Scalar Scalar;
    typedef typename NumTraits<Scalar>::Real RealScalar;

    using Base::rows;
    using Base::cols;
    using Base::size;
    using Base::derived;

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index innerStride() const
    {
      return derived().innerStride();
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index outerStride() const
    {
      return derived().outerStride();
    }

    EIGEN_CONSTEXPR inline Index stride() const
    {
      return Derived::IsVectorAtCompileTime ? innerStride() : outerStride();
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index rowStride() const
    {
      return Derived::IsRowMajor ? outerStride() : innerStride();
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index colStride() const
    {
      return Derived::IsRowMajor ? innerStride() : outerStride();
    }
};

template<typename Derived>
class DenseCoeffsBase<Derived, DirectWriteAccessors>
  : public DenseCoeffsBase<Derived, WriteAccessors>
{
  public:

    typedef DenseCoeffsBase<Derived, WriteAccessors> Base;
    typedef typename internal::traits<Derived>::Scalar Scalar;
    typedef typename NumTraits<Scalar>::Real RealScalar;

    using Base::rows;
    using Base::cols;
    using Base::size;
    using Base::derived;

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index innerStride() const EIGEN_NOEXCEPT
    {
      return derived().innerStride();
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index outerStride() const EIGEN_NOEXCEPT
    {
      return derived().outerStride();
    }

    EIGEN_CONSTEXPR inline Index stride() const EIGEN_NOEXCEPT
    {
      return Derived::IsVectorAtCompileTime ? innerStride() : outerStride();
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index rowStride() const EIGEN_NOEXCEPT
    {
      return Derived::IsRowMajor ? outerStride() : innerStride();
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index colStride() const EIGEN_NOEXCEPT
    {
      return Derived::IsRowMajor ? innerStride() : outerStride();
    }
};

namespace internal {

template<int Alignment, typename Derived, bool JustReturnZero>
struct first_aligned_impl
{
  static EIGEN_CONSTEXPR inline Index run(const Derived&) EIGEN_NOEXCEPT
  { return 0; }
};

template<int Alignment, typename Derived>
struct first_aligned_impl<Alignment, Derived, false>
{
  static inline Index run(const Derived& m)
  {
    return internal::first_aligned<Alignment>(m.data(), m.size());
  }
};

template<int Alignment, typename Derived>
static inline Index first_aligned(const DenseBase<Derived>& m)
{
  enum { ReturnZero = (int(evaluator<Derived>::Alignment) >= Alignment) || !(Derived::Flags & DirectAccessBit) };
  return first_aligned_impl<Alignment, Derived, ReturnZero>::run(m.derived());
}

template<typename Derived>
static inline Index first_default_aligned(const DenseBase<Derived>& m)
{
  typedef typename Derived::Scalar Scalar;
  typedef typename packet_traits<Scalar>::type DefaultPacketType;
  return internal::first_aligned<int(unpacket_traits<DefaultPacketType>::alignment),Derived>(m);
}

template<typename Derived, bool HasDirectAccess = has_direct_access<Derived>::ret>
struct inner_stride_at_compile_time
{
  enum { ret = traits<Derived>::InnerStrideAtCompileTime };
};

template<typename Derived>
struct inner_stride_at_compile_time<Derived, false>
{
  enum { ret = 0 };
};

template<typename Derived, bool HasDirectAccess = has_direct_access<Derived>::ret>
struct outer_stride_at_compile_time
{
  enum { ret = traits<Derived>::OuterStrideAtCompileTime };
};

template<typename Derived>
struct outer_stride_at_compile_time<Derived, false>
{
  enum { ret = 0 };
};

} 

} 

#endif 
