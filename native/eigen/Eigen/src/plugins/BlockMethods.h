#ifndef EIGEN_PARSED_BY_DOXYGEN

typedef Block<Derived, internal::traits<Derived>::RowsAtCompileTime, 1, !IsRowMajor> ColXpr;
typedef const Block<const Derived, internal::traits<Derived>::RowsAtCompileTime, 1, !IsRowMajor> ConstColXpr;

typedef Block<Derived, 1, internal::traits<Derived>::ColsAtCompileTime, IsRowMajor> RowXpr;
typedef const Block<const Derived, 1, internal::traits<Derived>::ColsAtCompileTime, IsRowMajor> ConstRowXpr;

typedef Block<Derived, internal::traits<Derived>::RowsAtCompileTime, Dynamic, !IsRowMajor> ColsBlockXpr;
typedef const Block<const Derived, internal::traits<Derived>::RowsAtCompileTime, Dynamic, !IsRowMajor> ConstColsBlockXpr;

typedef Block<Derived, Dynamic, internal::traits<Derived>::ColsAtCompileTime, IsRowMajor> RowsBlockXpr;
typedef const Block<const Derived, Dynamic, internal::traits<Derived>::ColsAtCompileTime, IsRowMajor> ConstRowsBlockXpr;

template<int N> struct NColsBlockXpr { typedef Block<Derived, internal::traits<Derived>::RowsAtCompileTime, N, !IsRowMajor> Type; };
template<int N> struct ConstNColsBlockXpr { typedef const Block<const Derived, internal::traits<Derived>::RowsAtCompileTime, N, !IsRowMajor> Type; };

template<int N> struct NRowsBlockXpr { typedef Block<Derived, N, internal::traits<Derived>::ColsAtCompileTime, IsRowMajor> Type; };
template<int N> struct ConstNRowsBlockXpr { typedef const Block<const Derived, N, internal::traits<Derived>::ColsAtCompileTime, IsRowMajor> Type; };

typedef Block<Derived> BlockXpr;
typedef const Block<const Derived> ConstBlockXpr;

template<int Rows, int Cols> struct FixedBlockXpr { typedef Block<Derived,Rows,Cols> Type; };
template<int Rows, int Cols> struct ConstFixedBlockXpr { typedef Block<const Derived,Rows,Cols> Type; };

typedef VectorBlock<Derived> SegmentReturnType;
typedef const VectorBlock<const Derived> ConstSegmentReturnType;
template<int Size> struct FixedSegmentReturnType { typedef VectorBlock<Derived, Size> Type; };
template<int Size> struct ConstFixedSegmentReturnType { typedef const VectorBlock<const Derived, Size> Type; };

typedef Block<Derived,IsRowMajor?1:Dynamic,IsRowMajor?Dynamic:1,true>       InnerVectorReturnType;
typedef Block<const Derived,IsRowMajor?1:Dynamic,IsRowMajor?Dynamic:1,true> ConstInnerVectorReturnType;

typedef Block<Derived,Dynamic,Dynamic,true> InnerVectorsReturnType;
typedef Block<const Derived,Dynamic,Dynamic,true> ConstInnerVectorsReturnType;

#endif 

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename FixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
#else
typename FixedBlockXpr<...,...>::Type
#endif
block(Index startRow, Index startCol, NRowsType blockRows, NColsType blockCols)
{
  return typename FixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type(
            derived(), startRow, startCol, internal::get_runtime_value(blockRows), internal::get_runtime_value(blockCols));
}

template<typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstFixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
#else
const typename ConstFixedBlockXpr<...,...>::Type
#endif
block(Index startRow, Index startCol, NRowsType blockRows, NColsType blockCols) const
{
  return typename ConstFixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type(
            derived(), startRow, startCol, internal::get_runtime_value(blockRows), internal::get_runtime_value(blockCols));
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename FixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
#else
typename FixedBlockXpr<...,...>::Type
#endif
topRightCorner(NRowsType cRows, NColsType cCols)
{
  return typename FixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
            (derived(), 0, cols() - internal::get_runtime_value(cCols), internal::get_runtime_value(cRows), internal::get_runtime_value(cCols));
}

template<typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstFixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
#else
const typename ConstFixedBlockXpr<...,...>::Type
#endif
topRightCorner(NRowsType cRows, NColsType cCols) const
{
  return typename ConstFixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
            (derived(), 0, cols() - internal::get_runtime_value(cCols), internal::get_runtime_value(cRows), internal::get_runtime_value(cCols));
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedBlockXpr<CRows,CCols>::Type topRightCorner()
{
  return typename FixedBlockXpr<CRows,CCols>::Type(derived(), 0, cols() - CCols);
}

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const typename ConstFixedBlockXpr<CRows,CCols>::Type topRightCorner() const
{
  return typename ConstFixedBlockXpr<CRows,CCols>::Type(derived(), 0, cols() - CCols);
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedBlockXpr<CRows,CCols>::Type topRightCorner(Index cRows, Index cCols)
{
  return typename FixedBlockXpr<CRows,CCols>::Type(derived(), 0, cols() - cCols, cRows, cCols);
}

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const typename ConstFixedBlockXpr<CRows,CCols>::Type topRightCorner(Index cRows, Index cCols) const
{
  return typename ConstFixedBlockXpr<CRows,CCols>::Type(derived(), 0, cols() - cCols, cRows, cCols);
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename FixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
#else
typename FixedBlockXpr<...,...>::Type
#endif
topLeftCorner(NRowsType cRows, NColsType cCols)
{
  return typename FixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
            (derived(), 0, 0, internal::get_runtime_value(cRows), internal::get_runtime_value(cCols));
}

template<typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstFixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
#else
const typename ConstFixedBlockXpr<...,...>::Type
#endif
topLeftCorner(NRowsType cRows, NColsType cCols) const
{
  return typename ConstFixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
            (derived(), 0, 0, internal::get_runtime_value(cRows), internal::get_runtime_value(cCols));
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedBlockXpr<CRows,CCols>::Type topLeftCorner()
{
  return typename FixedBlockXpr<CRows,CCols>::Type(derived(), 0, 0);
}

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const typename ConstFixedBlockXpr<CRows,CCols>::Type topLeftCorner() const
{
  return typename ConstFixedBlockXpr<CRows,CCols>::Type(derived(), 0, 0);
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedBlockXpr<CRows,CCols>::Type topLeftCorner(Index cRows, Index cCols)
{
  return typename FixedBlockXpr<CRows,CCols>::Type(derived(), 0, 0, cRows, cCols);
}

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const typename ConstFixedBlockXpr<CRows,CCols>::Type topLeftCorner(Index cRows, Index cCols) const
{
  return typename ConstFixedBlockXpr<CRows,CCols>::Type(derived(), 0, 0, cRows, cCols);
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename FixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
#else
typename FixedBlockXpr<...,...>::Type
#endif
bottomRightCorner(NRowsType cRows, NColsType cCols)
{
  return typename FixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
            (derived(), rows() - internal::get_runtime_value(cRows), cols() - internal::get_runtime_value(cCols),
                        internal::get_runtime_value(cRows), internal::get_runtime_value(cCols));
}

template<typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstFixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
#else
const typename ConstFixedBlockXpr<...,...>::Type
#endif
bottomRightCorner(NRowsType cRows, NColsType cCols) const
{
  return typename ConstFixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
            (derived(), rows() - internal::get_runtime_value(cRows), cols() - internal::get_runtime_value(cCols),
                        internal::get_runtime_value(cRows), internal::get_runtime_value(cCols));
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedBlockXpr<CRows,CCols>::Type bottomRightCorner()
{
  return typename FixedBlockXpr<CRows,CCols>::Type(derived(), rows() - CRows, cols() - CCols);
}

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const typename ConstFixedBlockXpr<CRows,CCols>::Type bottomRightCorner() const
{
  return typename ConstFixedBlockXpr<CRows,CCols>::Type(derived(), rows() - CRows, cols() - CCols);
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedBlockXpr<CRows,CCols>::Type bottomRightCorner(Index cRows, Index cCols)
{
  return typename FixedBlockXpr<CRows,CCols>::Type(derived(), rows() - cRows, cols() - cCols, cRows, cCols);
}

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const typename ConstFixedBlockXpr<CRows,CCols>::Type bottomRightCorner(Index cRows, Index cCols) const
{
  return typename ConstFixedBlockXpr<CRows,CCols>::Type(derived(), rows() - cRows, cols() - cCols, cRows, cCols);
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename FixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
#else
typename FixedBlockXpr<...,...>::Type
#endif
bottomLeftCorner(NRowsType cRows, NColsType cCols)
{
  return typename FixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
            (derived(), rows() - internal::get_runtime_value(cRows), 0,
                        internal::get_runtime_value(cRows), internal::get_runtime_value(cCols));
}

template<typename NRowsType, typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename ConstFixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
#else
typename ConstFixedBlockXpr<...,...>::Type
#endif
bottomLeftCorner(NRowsType cRows, NColsType cCols) const
{
  return typename ConstFixedBlockXpr<internal::get_fixed_value<NRowsType>::value,internal::get_fixed_value<NColsType>::value>::Type
            (derived(), rows() - internal::get_runtime_value(cRows), 0,
                        internal::get_runtime_value(cRows), internal::get_runtime_value(cCols));
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedBlockXpr<CRows,CCols>::Type bottomLeftCorner()
{
  return typename FixedBlockXpr<CRows,CCols>::Type(derived(), rows() - CRows, 0);
}

template<int CRows, int CCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const typename ConstFixedBlockXpr<CRows,CCols>::Type bottomLeftCorner() const
{
  return typename ConstFixedBlockXpr<CRows,CCols>::Type(derived(), rows() - CRows, 0);
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<int CRows, int CCols>
EIGEN_STRONG_INLINE
typename FixedBlockXpr<CRows,CCols>::Type bottomLeftCorner(Index cRows, Index cCols)
{
  return typename FixedBlockXpr<CRows,CCols>::Type(derived(), rows() - cRows, 0, cRows, cCols);
}

template<int CRows, int CCols>
EIGEN_STRONG_INLINE
const typename ConstFixedBlockXpr<CRows,CCols>::Type bottomLeftCorner(Index cRows, Index cCols) const
{
  return typename ConstFixedBlockXpr<CRows,CCols>::Type(derived(), rows() - cRows, 0, cRows, cCols);
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(row-major)

template<typename NRowsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename NRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
#else
typename NRowsBlockXpr<...>::Type
#endif
topRows(NRowsType n)
{
  return typename NRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
            (derived(), 0, 0, internal::get_runtime_value(n), cols());
}

template<typename NRowsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstNRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
#else
const typename ConstNRowsBlockXpr<...>::Type
#endif
topRows(NRowsType n) const
{
  return typename ConstNRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
            (derived(), 0, 0, internal::get_runtime_value(n), cols());
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(row-major)

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename NRowsBlockXpr<N>::Type topRows(Index n = N)
{
  return typename NRowsBlockXpr<N>::Type(derived(), 0, 0, n, cols());
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename ConstNRowsBlockXpr<N>::Type topRows(Index n = N) const
{
  return typename ConstNRowsBlockXpr<N>::Type(derived(), 0, 0, n, cols());
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(row-major)

template<typename NRowsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename NRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
#else
typename NRowsBlockXpr<...>::Type
#endif
bottomRows(NRowsType n)
{
  return typename NRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
            (derived(), rows() - internal::get_runtime_value(n), 0, internal::get_runtime_value(n), cols());
}

template<typename NRowsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstNRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
#else
const typename ConstNRowsBlockXpr<...>::Type
#endif
bottomRows(NRowsType n) const
{
  return typename ConstNRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
            (derived(), rows() - internal::get_runtime_value(n), 0, internal::get_runtime_value(n), cols());
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(row-major)

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename NRowsBlockXpr<N>::Type bottomRows(Index n = N)
{
  return typename NRowsBlockXpr<N>::Type(derived(), rows() - n, 0, n, cols());
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename ConstNRowsBlockXpr<N>::Type bottomRows(Index n = N) const
{
  return typename ConstNRowsBlockXpr<N>::Type(derived(), rows() - n, 0, n, cols());
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(row-major)

template<typename NRowsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename NRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
#else
typename NRowsBlockXpr<...>::Type
#endif
middleRows(Index startRow, NRowsType n)
{
  return typename NRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
            (derived(), startRow, 0, internal::get_runtime_value(n), cols());
}

template<typename NRowsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstNRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
#else
const typename ConstNRowsBlockXpr<...>::Type
#endif
middleRows(Index startRow, NRowsType n) const
{
  return typename ConstNRowsBlockXpr<internal::get_fixed_value<NRowsType>::value>::Type
            (derived(), startRow, 0, internal::get_runtime_value(n), cols());
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(row-major)

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename NRowsBlockXpr<N>::Type middleRows(Index startRow, Index n = N)
{
  return typename NRowsBlockXpr<N>::Type(derived(), startRow, 0, n, cols());
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename ConstNRowsBlockXpr<N>::Type middleRows(Index startRow, Index n = N) const
{
  return typename ConstNRowsBlockXpr<N>::Type(derived(), startRow, 0, n, cols());
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(column-major)

template<typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename NColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
#else
typename NColsBlockXpr<...>::Type
#endif
leftCols(NColsType n)
{
  return typename NColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
            (derived(), 0, 0, rows(), internal::get_runtime_value(n));
}

template<typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstNColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
#else
const typename ConstNColsBlockXpr<...>::Type
#endif
leftCols(NColsType n) const
{
  return typename ConstNColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
            (derived(), 0, 0, rows(), internal::get_runtime_value(n));
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(column-major)

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename NColsBlockXpr<N>::Type leftCols(Index n = N)
{
  return typename NColsBlockXpr<N>::Type(derived(), 0, 0, rows(), n);
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename ConstNColsBlockXpr<N>::Type leftCols(Index n = N) const
{
  return typename ConstNColsBlockXpr<N>::Type(derived(), 0, 0, rows(), n);
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(column-major)

template<typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename NColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
#else
typename NColsBlockXpr<...>::Type
#endif
rightCols(NColsType n)
{
  return typename NColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
            (derived(), 0, cols() - internal::get_runtime_value(n), rows(), internal::get_runtime_value(n));
}

template<typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstNColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
#else
const typename ConstNColsBlockXpr<...>::Type
#endif
rightCols(NColsType n) const
{
  return typename ConstNColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
            (derived(), 0, cols() - internal::get_runtime_value(n), rows(), internal::get_runtime_value(n));
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(column-major)

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename NColsBlockXpr<N>::Type rightCols(Index n = N)
{
  return typename NColsBlockXpr<N>::Type(derived(), 0, cols() - n, rows(), n);
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename ConstNColsBlockXpr<N>::Type rightCols(Index n = N) const
{
  return typename ConstNColsBlockXpr<N>::Type(derived(), 0, cols() - n, rows(), n);
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(column-major)

template<typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename NColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
#else
typename NColsBlockXpr<...>::Type
#endif
middleCols(Index startCol, NColsType numCols)
{
  return typename NColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
            (derived(), 0, startCol, rows(), internal::get_runtime_value(numCols));
}

template<typename NColsType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstNColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
#else
const typename ConstNColsBlockXpr<...>::Type
#endif
middleCols(Index startCol, NColsType numCols) const
{
  return typename ConstNColsBlockXpr<internal::get_fixed_value<NColsType>::value>::Type
            (derived(), 0, startCol, rows(), internal::get_runtime_value(numCols));
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(column-major)

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename NColsBlockXpr<N>::Type middleCols(Index startCol, Index n = N)
{
  return typename NColsBlockXpr<N>::Type(derived(), 0, startCol, rows(), n);
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename ConstNColsBlockXpr<N>::Type middleCols(Index startCol, Index n = N) const
{
  return typename ConstNColsBlockXpr<N>::Type(derived(), 0, startCol, rows(), n);
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<int NRows, int NCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedBlockXpr<NRows,NCols>::Type block(Index startRow, Index startCol)
{
  return typename FixedBlockXpr<NRows,NCols>::Type(derived(), startRow, startCol);
}

template<int NRows, int NCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const typename ConstFixedBlockXpr<NRows,NCols>::Type block(Index startRow, Index startCol) const
{
  return typename ConstFixedBlockXpr<NRows,NCols>::Type(derived(), startRow, startCol);
}

EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL

template<int NRows, int NCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedBlockXpr<NRows,NCols>::Type block(Index startRow, Index startCol,
                                                  Index blockRows, Index blockCols)
{
  return typename FixedBlockXpr<NRows,NCols>::Type(derived(), startRow, startCol, blockRows, blockCols);
}

template<int NRows, int NCols>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const typename ConstFixedBlockXpr<NRows,NCols>::Type block(Index startRow, Index startCol,
                                                              Index blockRows, Index blockCols) const
{
  return typename ConstFixedBlockXpr<NRows,NCols>::Type(derived(), startRow, startCol, blockRows, blockCols);
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(column-major)

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
ColXpr col(Index i)
{
  return ColXpr(derived(), i);
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
ConstColXpr col(Index i) const
{
  return ConstColXpr(derived(), i);
}

EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(row-major)

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
RowXpr row(Index i)
{
  return RowXpr(derived(), i);
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
ConstRowXpr row(Index i) const
{
  return ConstRowXpr(derived(), i);
}

template<typename NType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename FixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
#else
typename FixedSegmentReturnType<...>::Type
#endif
segment(Index start, NType n)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename FixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
            (derived(), start, internal::get_runtime_value(n));
}

template<typename NType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstFixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
#else
const typename ConstFixedSegmentReturnType<...>::Type
#endif
segment(Index start, NType n) const
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename ConstFixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
            (derived(), start, internal::get_runtime_value(n));
}

template<typename NType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename FixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
#else
typename FixedSegmentReturnType<...>::Type
#endif
head(NType n)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename FixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
              (derived(), 0, internal::get_runtime_value(n));
}

template<typename NType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstFixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
#else
const typename ConstFixedSegmentReturnType<...>::Type
#endif
head(NType n) const
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename ConstFixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
            (derived(), 0, internal::get_runtime_value(n));
}

template<typename NType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
typename FixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
#else
typename FixedSegmentReturnType<...>::Type
#endif
tail(NType n)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename FixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
            (derived(), this->size() - internal::get_runtime_value(n), internal::get_runtime_value(n));
}

template<typename NType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
#ifndef EIGEN_PARSED_BY_DOXYGEN
const typename ConstFixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
#else
const typename ConstFixedSegmentReturnType<...>::Type
#endif
tail(NType n) const
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename ConstFixedSegmentReturnType<internal::get_fixed_value<NType>::value>::Type
            (derived(), this->size() - internal::get_runtime_value(n), internal::get_runtime_value(n));
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedSegmentReturnType<N>::Type segment(Index start, Index n = N)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename FixedSegmentReturnType<N>::Type(derived(), start, n);
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename ConstFixedSegmentReturnType<N>::Type segment(Index start, Index n = N) const
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename ConstFixedSegmentReturnType<N>::Type(derived(), start, n);
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedSegmentReturnType<N>::Type head(Index n = N)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename FixedSegmentReturnType<N>::Type(derived(), 0, n);
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename ConstFixedSegmentReturnType<N>::Type head(Index n = N) const
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename ConstFixedSegmentReturnType<N>::Type(derived(), 0, n);
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename FixedSegmentReturnType<N>::Type tail(Index n = N)
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename FixedSegmentReturnType<N>::Type(derived(), size() - n);
}

template<int N>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename ConstFixedSegmentReturnType<N>::Type tail(Index n = N) const
{
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived)
  return typename ConstFixedSegmentReturnType<N>::Type(derived(), size() - n);
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
InnerVectorReturnType innerVector(Index outer)
{ return InnerVectorReturnType(derived(), outer); }

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const ConstInnerVectorReturnType innerVector(Index outer) const
{ return ConstInnerVectorReturnType(derived(), outer); }

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
InnerVectorsReturnType
innerVectors(Index outerStart, Index outerSize)
{
  return Block<Derived,Dynamic,Dynamic,true>(derived(),
                                             IsRowMajor ? outerStart : 0, IsRowMajor ? 0 : outerStart,
                                             IsRowMajor ? outerSize : rows(), IsRowMajor ? cols() : outerSize);

}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const ConstInnerVectorsReturnType
innerVectors(Index outerStart, Index outerSize) const
{
  return Block<const Derived,Dynamic,Dynamic,true>(derived(),
                                                  IsRowMajor ? outerStart : 0, IsRowMajor ? 0 : outerStart,
                                                  IsRowMajor ? outerSize : rows(), IsRowMajor ? cols() : outerSize);

}

template<DirectionType Direction>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename internal::conditional<Direction==Vertical,ColXpr,RowXpr>::type
subVector(Index i)
{
  return typename internal::conditional<Direction==Vertical,ColXpr,RowXpr>::type(derived(),i);
}

template<DirectionType Direction>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typename internal::conditional<Direction==Vertical,ConstColXpr,ConstRowXpr>::type
subVector(Index i) const
{
  return typename internal::conditional<Direction==Vertical,ConstColXpr,ConstRowXpr>::type(derived(),i);
}

template<DirectionType Direction>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
Index subVectors() const
{ return (Direction==Vertical)?cols():rows(); }
