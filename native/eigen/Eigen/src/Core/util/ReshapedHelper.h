#ifndef EIGEN_RESHAPED_HELPER_H
#define EIGEN_RESHAPED_HELPER_H

namespace Eigen {

enum AutoSize_t   { AutoSize };
const int AutoOrder = 2;

namespace internal {

template<typename SizeType,typename OtherSize, int TotalSize>
struct get_compiletime_reshape_size {
  enum { value = get_fixed_value<SizeType>::value };
};

template<typename SizeType>
Index get_runtime_reshape_size(SizeType size, Index , Index ) {
  return internal::get_runtime_value(size);
}

template<typename OtherSize, int TotalSize>
struct get_compiletime_reshape_size<AutoSize_t,OtherSize,TotalSize> {
  enum {
    other_size = get_fixed_value<OtherSize>::value,
    value = (TotalSize==Dynamic || other_size==Dynamic) ? Dynamic : TotalSize / other_size };
};

inline Index get_runtime_reshape_size(AutoSize_t , Index other, Index total) {
  return total/other;
}

template<int Flags, int Order>
struct get_compiletime_reshape_order {
  enum { value = Order == AutoOrder ? Flags & RowMajorBit : Order };
};

}

} 

#endif 
