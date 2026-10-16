#ifndef EIGEN_VECTORBLOCK_H
#define EIGEN_VECTORBLOCK_H

namespace Eigen { 

namespace internal {
template<typename VectorType, int Size>
struct traits<VectorBlock<VectorType, Size> >
  : public traits<Block<VectorType,
                     traits<VectorType>::Flags & RowMajorBit ? 1 : Size,
                     traits<VectorType>::Flags & RowMajorBit ? Size : 1> >
{
};
}

template<typename VectorType, int Size> class VectorBlock
  : public Block<VectorType,
                     internal::traits<VectorType>::Flags & RowMajorBit ? 1 : Size,
                     internal::traits<VectorType>::Flags & RowMajorBit ? Size : 1>
{
    typedef Block<VectorType,
                     internal::traits<VectorType>::Flags & RowMajorBit ? 1 : Size,
                     internal::traits<VectorType>::Flags & RowMajorBit ? Size : 1> Base;
    enum {
      IsColVector = !(internal::traits<VectorType>::Flags & RowMajorBit)
    };
  public:
    EIGEN_DENSE_PUBLIC_INTERFACE(VectorBlock)

    using Base::operator=;

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    VectorBlock(VectorType& vector, Index start, Index size)
      : Base(vector,
             IsColVector ? start : 0, IsColVector ? 0 : start,
             IsColVector ? size  : 1, IsColVector ? 1 : size)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBlock);
    }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    VectorBlock(VectorType& vector, Index start)
      : Base(vector, IsColVector ? start : 0, IsColVector ? 0 : start)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBlock);
    }
};

} 

#endif 
