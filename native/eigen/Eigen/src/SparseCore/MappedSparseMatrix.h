#ifndef EIGEN_MAPPED_SPARSEMATRIX_H
#define EIGEN_MAPPED_SPARSEMATRIX_H

namespace Eigen {

namespace internal {
template<typename Scalar_, int _Flags, typename StorageIndex_>
struct traits<MappedSparseMatrix<Scalar_, _Flags, StorageIndex_> > : traits<SparseMatrix<Scalar_, _Flags, StorageIndex_> >
{};
} 

template<typename Scalar_, int _Flags, typename StorageIndex_>
class MappedSparseMatrix
  : public Map<SparseMatrix<Scalar_, _Flags, StorageIndex_> >
{
    typedef Map<SparseMatrix<Scalar_, _Flags, StorageIndex_> > Base;

  public:
    
    typedef typename Base::StorageIndex StorageIndex;
    typedef typename Base::Scalar Scalar;

    inline MappedSparseMatrix(Index rows, Index cols, Index nnz, StorageIndex* outerIndexPtr, StorageIndex* innerIndexPtr, Scalar* valuePtr, StorageIndex* innerNonZeroPtr = 0)
      : Base(rows, cols, nnz, outerIndexPtr, innerIndexPtr, valuePtr, innerNonZeroPtr)
    {}

    inline ~MappedSparseMatrix() {}
};

namespace internal {

template<typename Scalar_, int Options_, typename StorageIndex_>
struct evaluator<MappedSparseMatrix<Scalar_,Options_,StorageIndex_> >
  : evaluator<SparseCompressedBase<MappedSparseMatrix<Scalar_,Options_,StorageIndex_> > >
{
  typedef MappedSparseMatrix<Scalar_,Options_,StorageIndex_> XprType;
  typedef evaluator<SparseCompressedBase<XprType> > Base;
  
  evaluator() : Base() {}
  explicit evaluator(const XprType &mat) : Base(mat) {}
};

}

} 

#endif 
