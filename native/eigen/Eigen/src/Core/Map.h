#ifndef EIGEN_MAP_H
#define EIGEN_MAP_H

namespace Eigen {

namespace internal {
template<typename PlainObjectType, int MapOptions, typename StrideType>
struct traits<Map<PlainObjectType, MapOptions, StrideType> >
  : public traits<PlainObjectType>
{
  typedef traits<PlainObjectType> TraitsBase;
  enum {
    PlainObjectTypeInnerSize = ((traits<PlainObjectType>::Flags&RowMajorBit)==RowMajorBit)
                             ? PlainObjectType::ColsAtCompileTime
                             : PlainObjectType::RowsAtCompileTime,

    InnerStrideAtCompileTime = StrideType::InnerStrideAtCompileTime == 0
                             ? int(PlainObjectType::InnerStrideAtCompileTime)
                             : int(StrideType::InnerStrideAtCompileTime),
    OuterStrideAtCompileTime = StrideType::OuterStrideAtCompileTime == 0
                             ? (InnerStrideAtCompileTime==Dynamic || PlainObjectTypeInnerSize==Dynamic
                                ? Dynamic
                                : int(InnerStrideAtCompileTime) * int(PlainObjectTypeInnerSize))
                             : int(StrideType::OuterStrideAtCompileTime),
    Alignment = int(MapOptions)&int(AlignedMask),
    Flags0 = TraitsBase::Flags & (~NestByRefBit),
    Flags = is_lvalue<PlainObjectType>::value ? int(Flags0) : (int(Flags0) & ~LvalueBit)
  };
private:
  enum { Options }; 
};
}

template<typename PlainObjectType, int MapOptions, typename StrideType> class Map
  : public MapBase<Map<PlainObjectType, MapOptions, StrideType> >
{
  public:

    typedef MapBase<Map> Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(Map)

    typedef typename Base::PointerType PointerType;
    typedef PointerType PointerArgType;
    EIGEN_DEVICE_FUNC
    inline PointerType cast_to_pointer_type(PointerArgType ptr) { return ptr; }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index innerStride() const
    {
      return StrideType::InnerStrideAtCompileTime != 0 ? m_stride.inner() : 1;
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index outerStride() const
    {
      return StrideType::OuterStrideAtCompileTime != 0 ? m_stride.outer()
           : internal::traits<Map>::OuterStrideAtCompileTime != Dynamic ? Index(internal::traits<Map>::OuterStrideAtCompileTime)
           : IsVectorAtCompileTime ? (this->size() * innerStride())
           : int(Flags)&RowMajorBit ? (this->cols() * innerStride())
           : (this->rows() * innerStride());
    }

    EIGEN_DEVICE_FUNC
    explicit inline Map(PointerArgType dataPtr, const StrideType& stride = StrideType())
      : Base(cast_to_pointer_type(dataPtr)), m_stride(stride)
    {
      PlainObjectType::Base::_check_template_params();
    }

    EIGEN_DEVICE_FUNC
    inline Map(PointerArgType dataPtr, Index size, const StrideType& stride = StrideType())
      : Base(cast_to_pointer_type(dataPtr), size), m_stride(stride)
    {
      PlainObjectType::Base::_check_template_params();
    }

    EIGEN_DEVICE_FUNC
    inline Map(PointerArgType dataPtr, Index rows, Index cols, const StrideType& stride = StrideType())
      : Base(cast_to_pointer_type(dataPtr), rows, cols), m_stride(stride)
    {
      PlainObjectType::Base::_check_template_params();
    }

    EIGEN_INHERIT_ASSIGNMENT_OPERATORS(Map)

  protected:
    StrideType m_stride;
};

} 

#endif 
