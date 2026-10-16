#ifndef EIGEN_ORDERING_H
#define EIGEN_ORDERING_H

namespace Eigen {
  
#include "Eigen_Colamd.h"

namespace internal {

template<typename MatrixType> 
void ordering_helper_at_plus_a(const MatrixType& A, MatrixType& symmat)
{
  MatrixType C;
  C = A.transpose(); 
  for (int i = 0; i < C.rows(); i++) 
  {
      for (typename MatrixType::InnerIterator it(C, i); it; ++it)
        it.valueRef() = typename MatrixType::Scalar(0);
  }
  symmat = C + A;
}
    
}

template <typename StorageIndex>
class AMDOrdering
{
  public:
    typedef PermutationMatrix<Dynamic, Dynamic, StorageIndex> PermutationType;

    template <typename MatrixType>
    void operator()(const MatrixType& mat, PermutationType& perm)
    {
      
      SparseMatrix<typename MatrixType::Scalar, ColMajor, StorageIndex> symm;
      internal::ordering_helper_at_plus_a(mat,symm); 

      internal::minimum_degree_ordering(symm, perm);
    }

    template <typename SrcType, unsigned int SrcUpLo> 
    void operator()(const SparseSelfAdjointView<SrcType, SrcUpLo>& mat, PermutationType& perm)
    { 
      SparseMatrix<typename SrcType::Scalar, ColMajor, StorageIndex> C; C = mat;

      internal::minimum_degree_ordering(C, perm);
    }
};

template <typename StorageIndex>
class NaturalOrdering
{
  public:
    typedef PermutationMatrix<Dynamic, Dynamic, StorageIndex> PermutationType;

    template <typename MatrixType>
    void operator()(const MatrixType& , PermutationType& perm)
    {
      perm.resize(0); 
    }
    
};

template<typename StorageIndex>
class COLAMDOrdering
{
  public:
    typedef PermutationMatrix<Dynamic, Dynamic, StorageIndex> PermutationType; 
    typedef Matrix<StorageIndex, Dynamic, 1> IndexVector;

    template <typename MatrixType>
    void operator() (const MatrixType& mat, PermutationType& perm)
    {
      eigen_assert(mat.isCompressed() && "COLAMDOrdering requires a sparse matrix in compressed mode. Call .makeCompressed() before passing it to COLAMDOrdering");
      
      StorageIndex m = StorageIndex(mat.rows());
      StorageIndex n = StorageIndex(mat.cols());
      StorageIndex nnz = StorageIndex(mat.nonZeros());
      
      StorageIndex Alen = internal::Colamd::recommended(nnz, m, n); 
      
      double knobs [internal::Colamd::NKnobs]; 
      StorageIndex stats [internal::Colamd::NStats];
      internal::Colamd::set_defaults(knobs);
      
      IndexVector p(n+1), A(Alen); 
      for(StorageIndex i=0; i <= n; i++)   p(i) = mat.outerIndexPtr()[i];
      for(StorageIndex i=0; i < nnz; i++)  A(i) = mat.innerIndexPtr()[i];
      
      StorageIndex info = internal::Colamd::compute_ordering(m, n, Alen, A.data(), p.data(), knobs, stats); 
      EIGEN_UNUSED_VARIABLE(info);
      eigen_assert( info && "COLAMD failed " );
      
      perm.resize(n);
      for (StorageIndex i = 0; i < n; i++) perm.indices()(p(i)) = i;
    }
};

} 

#endif
