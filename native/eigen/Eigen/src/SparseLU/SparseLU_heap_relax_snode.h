#ifndef SPARSELU_HEAP_RELAX_SNODE_H
#define SPARSELU_HEAP_RELAX_SNODE_H

namespace Eigen {
namespace internal {

template <typename Scalar, typename StorageIndex>
void SparseLUImpl<Scalar,StorageIndex>::heap_relax_snode (const Index n, IndexVector& et, const Index relax_columns, IndexVector& descendants, IndexVector& relax_end)
{

  IndexVector post;
  internal::treePostorder(StorageIndex(n), et, post); 
  IndexVector inv_post(n+1); 
  for (StorageIndex i = 0; i < n+1; ++i) inv_post(post(i)) = i; 

  IndexVector iwork(n);
  IndexVector et_save(n+1);
  for (Index i = 0; i < n; ++i)
  {
    iwork(post(i)) = post(et(i));
  }
  et_save = et; 
  et = iwork; 

  relax_end.setConstant(emptyIdxLU);
  Index j, parent; 
  descendants.setZero();
  for (j = 0; j < n; j++) 
  {
    parent = et(j);
    if (parent != n) 
      descendants(parent) += descendants(j) + 1;
  }
  
  Index snode_start; 
  StorageIndex k;
  Index nsuper_et_post = 0; 
  Index nsuper_et = 0; 
  StorageIndex l; 
  for (j = 0; j < n; )
  {
    parent = et(j);
    snode_start = j; 
    while ( parent != n && descendants(parent) < relax_columns ) 
    {
      j = parent; 
      parent = et(j);
    }
    
    ++nsuper_et_post;
    k = StorageIndex(n);
    for (Index i = snode_start; i <= j; ++i)
      k = (std::min)(k, inv_post(i));
    l = inv_post(j);
    if ( (l - k) == (j - snode_start) )  
    {
      
      relax_end(k) = l; 
      ++nsuper_et; 
    }
    else 
    {
      for (Index i = snode_start; i <= j; ++i) 
      {
        l = inv_post(i);
        if (descendants(i) == 0) 
        {
          relax_end(l) = l;
          ++nsuper_et;
        }
      }
    }
    j++;
    
    while (descendants(j) != 0 && j < n) j++;
  } 

  et = et_save; 
}

} 

} 
#endif 
