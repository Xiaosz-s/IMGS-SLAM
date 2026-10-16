#ifdef EIGEN_PARSED_BY_DOXYGEN

template<int Order = ColMajor, typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC
inline Reshaped<Derived,...>
reshaped(NRowsType nRows, NColsType nCols);

template<int Order = ColMajor, typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC
inline const Reshaped<const Derived,...>
reshaped(NRowsType nRows, NColsType nCols) const;

template<int Order = ColMajor>
EIGEN_DEVICE_FUNC
inline Reshaped<Derived,...>
reshaped();

template<int Order = ColMajor>
EIGEN_DEVICE_FUNC
inline const Reshaped<const Derived,...>
reshaped() const;

#else

#ifndef EIGEN_RESHAPED_METHOD_2ND_PASS
#define EIGEN_RESHAPED_METHOD_CONST const
#else
#define EIGEN_RESHAPED_METHOD_CONST
#endif

#ifndef EIGEN_RESHAPED_METHOD_2ND_PASS

#endif

template<typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC
inline Reshaped<EIGEN_RESHAPED_METHOD_CONST Derived,
                internal::get_compiletime_reshape_size<NRowsType,NColsType,SizeAtCompileTime>::value,
                internal::get_compiletime_reshape_size<NColsType,NRowsType,SizeAtCompileTime>::value>
reshaped(NRowsType nRows, NColsType nCols) EIGEN_RESHAPED_METHOD_CONST
{
  return Reshaped<EIGEN_RESHAPED_METHOD_CONST Derived,
                  internal::get_compiletime_reshape_size<NRowsType,NColsType,SizeAtCompileTime>::value,
                  internal::get_compiletime_reshape_size<NColsType,NRowsType,SizeAtCompileTime>::value>
                (derived(),
                 internal::get_runtime_reshape_size(nRows,internal::get_runtime_value(nCols),size()),
                 internal::get_runtime_reshape_size(nCols,internal::get_runtime_value(nRows),size()));
}

template<int Order, typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC
inline Reshaped<EIGEN_RESHAPED_METHOD_CONST Derived,
                internal::get_compiletime_reshape_size<NRowsType,NColsType,SizeAtCompileTime>::value,
                internal::get_compiletime_reshape_size<NColsType,NRowsType,SizeAtCompileTime>::value,
                internal::get_compiletime_reshape_order<Flags,Order>::value>
reshaped(NRowsType nRows, NColsType nCols) EIGEN_RESHAPED_METHOD_CONST
{
  return Reshaped<EIGEN_RESHAPED_METHOD_CONST Derived,
                  internal::get_compiletime_reshape_size<NRowsType,NColsType,SizeAtCompileTime>::value,
                  internal::get_compiletime_reshape_size<NColsType,NRowsType,SizeAtCompileTime>::value,
                  internal::get_compiletime_reshape_order<Flags,Order>::value>
                (derived(),
                 internal::get_runtime_reshape_size(nRows,internal::get_runtime_value(nCols),size()),
                 internal::get_runtime_reshape_size(nCols,internal::get_runtime_value(nRows),size()));
}

EIGEN_DEVICE_FUNC
inline Reshaped<EIGEN_RESHAPED_METHOD_CONST Derived,SizeAtCompileTime,1>
reshaped() EIGEN_RESHAPED_METHOD_CONST
{
  return Reshaped<EIGEN_RESHAPED_METHOD_CONST Derived,SizeAtCompileTime,1>(derived(),size(),1);
}

template<int Order>
EIGEN_DEVICE_FUNC
inline Reshaped<EIGEN_RESHAPED_METHOD_CONST Derived, SizeAtCompileTime, 1,
                internal::get_compiletime_reshape_order<Flags,Order>::value>
reshaped() EIGEN_RESHAPED_METHOD_CONST
{
  EIGEN_STATIC_ASSERT(Order==RowMajor || Order==ColMajor || Order==AutoOrder, INVALID_TEMPLATE_PARAMETER);
  return Reshaped<EIGEN_RESHAPED_METHOD_CONST Derived, SizeAtCompileTime, 1,
                  internal::get_compiletime_reshape_order<Flags,Order>::value>
                (derived(), size(), 1);
}

#undef EIGEN_RESHAPED_METHOD_CONST

#ifndef EIGEN_RESHAPED_METHOD_2ND_PASS
#define EIGEN_RESHAPED_METHOD_2ND_PASS
#include "ReshapedMethods.h"
#undef EIGEN_RESHAPED_METHOD_2ND_PASS
#endif

#endif 
