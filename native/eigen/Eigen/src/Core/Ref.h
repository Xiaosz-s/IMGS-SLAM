#ifndef EIGEN_REF_H
#define EIGEN_REF_H

namespace Eigen {

namespace internal {

template<typename _PlainObjectType, int Options_, typename _StrideType>
struct traits<Ref<_PlainObjectType, Options_, _StrideType> >
  : public traits<Map<_PlainObjectType, Options_, _StrideType> >
{
  typedef _PlainObjectType PlainObjectType;
  typedef _StrideType StrideType;
  enum {
    Options = Options_,
    Flags = traits<Map<_PlainObjectType, Options_, _StrideType> >::Flags | NestByRefBit,
    Alignment = traits<Map<_PlainObjectType, Options_, _StrideType> >::Alignment
  };

  template<typename Derived> struct match {
    enum {
      IsVectorAtCompileTime = PlainObjectType::IsVectorAtCompileTime || Derived::IsVectorAtCompileTime,
      HasDirectAccess = internal::has_direct_access<Derived>::ret,
      StorageOrderMatch = IsVectorAtCompileTime || ((PlainObjectType::Flags&RowMajorBit)==(Derived::Flags&RowMajorBit)),
      InnerStrideMatch = int(StrideType::InnerStrideAtCompileTime)==int(Dynamic)
                      || int(StrideType::InnerStrideAtCompileTime)==int(Derived::InnerStrideAtCompileTime)
                      || (int(StrideType::InnerStrideAtCompileTime)==0 && int(Derived::InnerStrideAtCompileTime)==1),
      OuterStrideMatch = IsVectorAtCompileTime
                      || int(StrideType::OuterStrideAtCompileTime)==int(Dynamic) || int(StrideType::OuterStrideAtCompileTime)==int(Derived::OuterStrideAtCompileTime),

      DerivedAlignment = int(evaluator<Derived>::Alignment),
      AlignmentMatch = (int(traits<PlainObjectType>::Alignment)==int(Unaligned)) || (DerivedAlignment >= int(Alignment)), 
      ScalarTypeMatch = internal::is_same<typename PlainObjectType::Scalar, typename Derived::Scalar>::value,
      MatchAtCompileTime = HasDirectAccess && StorageOrderMatch && InnerStrideMatch && OuterStrideMatch && AlignmentMatch && ScalarTypeMatch
    };
    typedef typename internal::conditional<MatchAtCompileTime,internal::true_type,internal::false_type>::type type;
  };

};

template<typename Derived>
struct traits<RefBase<Derived> > : public traits<Derived> {};

}

template<typename Derived> class RefBase
 : public MapBase<Derived>
{
  typedef typename internal::traits<Derived>::PlainObjectType PlainObjectType;
  typedef typename internal::traits<Derived>::StrideType StrideType;

public:

  typedef MapBase<Derived> Base;
  EIGEN_DENSE_PUBLIC_INTERFACE(RefBase)

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR inline Index innerStride() const
  {
    return StrideType::InnerStrideAtCompileTime != 0 ? m_stride.inner() : 1;
  }

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR inline Index outerStride() const
  {
    return StrideType::OuterStrideAtCompileTime != 0 ? m_stride.outer()
         : IsVectorAtCompileTime ? this->size()
         : int(Flags)&RowMajorBit ? this->cols()
         : this->rows();
  }

  EIGEN_DEVICE_FUNC RefBase()
    : Base(0,RowsAtCompileTime==Dynamic?0:RowsAtCompileTime,ColsAtCompileTime==Dynamic?0:ColsAtCompileTime),
      
      m_stride(StrideType::OuterStrideAtCompileTime==Dynamic?0:StrideType::OuterStrideAtCompileTime,
               StrideType::InnerStrideAtCompileTime==Dynamic?0:StrideType::InnerStrideAtCompileTime)
  {}

  EIGEN_INHERIT_ASSIGNMENT_OPERATORS(RefBase)

protected:

  typedef Stride<StrideType::OuterStrideAtCompileTime,StrideType::InnerStrideAtCompileTime> StrideBase;

  static EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR Index resolveInnerStride(Index inner) {
    return inner == 0 ? 1 : inner;
  }

  static EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR Index resolveOuterStride(Index inner, Index outer, Index rows, Index cols, bool isVectorAtCompileTime, bool isRowMajor) {
    return outer == 0 ? isVectorAtCompileTime ? inner * rows * cols : isRowMajor ? inner * cols : inner * rows : outer;
  }

  template<typename Expression>
  EIGEN_DEVICE_FUNC bool construct(Expression& expr)
  {

    EIGEN_STATIC_ASSERT(
      EIGEN_PREDICATE_SAME_MATRIX_SIZE(PlainObjectType, Expression)
      
      || ( PlainObjectType::IsVectorAtCompileTime
            && ((int(PlainObjectType::RowsAtCompileTime)==Eigen::Dynamic
              || int(Expression::ColsAtCompileTime)==Eigen::Dynamic
              || int(PlainObjectType::RowsAtCompileTime)==int(Expression::ColsAtCompileTime))
            &&  (int(PlainObjectType::ColsAtCompileTime)==Eigen::Dynamic
              || int(Expression::RowsAtCompileTime)==Eigen::Dynamic
              || int(PlainObjectType::ColsAtCompileTime)==int(Expression::RowsAtCompileTime)))),
      YOU_MIXED_MATRICES_OF_DIFFERENT_SIZES
    )

    Index rows = expr.rows();
    Index cols = expr.cols();
    if(PlainObjectType::RowsAtCompileTime==1)
    {
      eigen_assert(expr.rows()==1 || expr.cols()==1);
      rows = 1;
      cols = expr.size();
    }
    else if(PlainObjectType::ColsAtCompileTime==1)
    {
      eigen_assert(expr.rows()==1 || expr.cols()==1);
      rows = expr.size();
      cols = 1;
    }
    
    eigen_assert(
      (PlainObjectType::RowsAtCompileTime == Dynamic) || (PlainObjectType::RowsAtCompileTime == rows));
    eigen_assert(
      (PlainObjectType::ColsAtCompileTime == Dynamic) || (PlainObjectType::ColsAtCompileTime == cols));

    const bool transpose = PlainObjectType::IsVectorAtCompileTime && (rows != expr.rows());
    
    const bool row_major = ((PlainObjectType::Flags)&RowMajorBit) != 0;
    const bool expr_row_major = (Expression::Flags&RowMajorBit) != 0;
    const bool storage_differs =  (row_major != expr_row_major);

    const bool swap_stride = (transpose != storage_differs);

    const Index expr_inner_actual = resolveInnerStride(expr.innerStride());
    const Index expr_outer_actual = resolveOuterStride(expr_inner_actual,
                                                       expr.outerStride(),
                                                       expr.rows(),
                                                       expr.cols(),
                                                       Expression::IsVectorAtCompileTime != 0,
                                                       expr_row_major);

    const bool row_vector = (rows == 1);
    const bool col_vector = (cols == 1);
    const Index inner_stride =
        ( (!row_major && row_vector) || (row_major && col_vector) ) ?
            ( StrideType::InnerStrideAtCompileTime > 0 ? Index(StrideType::InnerStrideAtCompileTime) : 1)
            : swap_stride ? expr_outer_actual : expr_inner_actual;

    const Index outer_stride =
      ( (!row_major && col_vector) || (row_major && row_vector) ) ?
          ( StrideType::OuterStrideAtCompileTime > 0 ? Index(StrideType::OuterStrideAtCompileTime) : rows * cols * inner_stride)
          : swap_stride ? expr_inner_actual : expr_outer_actual;

    const bool inner_valid = (StrideType::InnerStrideAtCompileTime == Dynamic)
        || (resolveInnerStride(Index(StrideType::InnerStrideAtCompileTime)) == inner_stride);
    if (!inner_valid) {
      return false;
    }

    const bool outer_valid = (StrideType::OuterStrideAtCompileTime == Dynamic)
        || (resolveOuterStride(
              inner_stride,
              Index(StrideType::OuterStrideAtCompileTime),
              rows, cols, PlainObjectType::IsVectorAtCompileTime != 0,
              row_major)
            == outer_stride);
    if (!outer_valid) {
      return false;
    }

    ::new (static_cast<Base*>(this)) Base(expr.data(), rows, cols);
    ::new (&m_stride) StrideBase(
      (StrideType::OuterStrideAtCompileTime == 0) ? 0 : outer_stride,
      (StrideType::InnerStrideAtCompileTime == 0) ? 0 : inner_stride );
    return true;
  }

  StrideBase m_stride;
};

template<typename PlainObjectType, int Options, typename StrideType> class Ref
  : public RefBase<Ref<PlainObjectType, Options, StrideType> >
{
  private:
    typedef internal::traits<Ref> Traits;
    template<typename Derived>
    EIGEN_DEVICE_FUNC inline Ref(const PlainObjectBase<Derived>& expr,
                                 typename internal::enable_if<bool(Traits::template match<Derived>::MatchAtCompileTime),Derived>::type* = 0);
  public:

    typedef RefBase<Ref> Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(Ref)

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    template<typename Derived>
    EIGEN_DEVICE_FUNC inline Ref(PlainObjectBase<Derived>& expr,
                                 typename internal::enable_if<bool(Traits::template match<Derived>::MatchAtCompileTime),Derived>::type* = 0)
    {
      EIGEN_STATIC_ASSERT(bool(Traits::template match<Derived>::MatchAtCompileTime), STORAGE_LAYOUT_DOES_NOT_MATCH);
      
      const bool success = Base::construct(expr.derived());
      EIGEN_UNUSED_VARIABLE(success)
      eigen_assert(success);
    }
    template<typename Derived>
    EIGEN_DEVICE_FUNC inline Ref(const DenseBase<Derived>& expr,
                                 typename internal::enable_if<bool(Traits::template match<Derived>::MatchAtCompileTime),Derived>::type* = 0)
    #else
    
    template<typename Derived>
    inline Ref(DenseBase<Derived>& expr)
    #endif
    {
      EIGEN_STATIC_ASSERT(bool(internal::is_lvalue<Derived>::value), THIS_EXPRESSION_IS_NOT_A_LVALUE__IT_IS_READ_ONLY);
      EIGEN_STATIC_ASSERT(bool(Traits::template match<Derived>::MatchAtCompileTime), STORAGE_LAYOUT_DOES_NOT_MATCH);
      EIGEN_STATIC_ASSERT(!Derived::IsPlainObjectBase,THIS_EXPRESSION_IS_NOT_A_LVALUE__IT_IS_READ_ONLY);
      
      const bool success = Base::construct(expr.const_cast_derived());
      EIGEN_UNUSED_VARIABLE(success)
      eigen_assert(success);
    }

    EIGEN_INHERIT_ASSIGNMENT_OPERATORS(Ref)

};

template<typename TPlainObjectType, int Options, typename StrideType> class Ref<const TPlainObjectType, Options, StrideType>
  : public RefBase<Ref<const TPlainObjectType, Options, StrideType> >
{
    typedef internal::traits<Ref> Traits;
  public:

    typedef RefBase<Ref> Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(Ref)

    template<typename Derived>
    EIGEN_DEVICE_FUNC inline Ref(const DenseBase<Derived>& expr,
                                 typename internal::enable_if<bool(Traits::template match<Derived>::ScalarTypeMatch),Derived>::type* = 0)
    {

      construct(expr.derived(), typename Traits::template match<Derived>::type());
    }

    EIGEN_DEVICE_FUNC inline Ref(const Ref& other) : Base(other) {
      
    }

    template<typename OtherRef>
    EIGEN_DEVICE_FUNC inline Ref(const RefBase<OtherRef>& other) {
      construct(other.derived(), typename Traits::template match<OtherRef>::type());
    }

  protected:

    template<typename Expression>
    EIGEN_DEVICE_FUNC void construct(const Expression& expr,internal::true_type)
    {
      
      if (!Base::construct(expr)) {
        construct(expr, internal::false_type());
      }
    }

    template<typename Expression>
    EIGEN_DEVICE_FUNC void construct(const Expression& expr, internal::false_type)
    {
      internal::call_assignment_no_alias(m_object,expr,internal::assign_op<Scalar,Scalar>());
      Base::construct(m_object);
    }

  protected:
    TPlainObjectType m_object;
};

} 

#endif 
