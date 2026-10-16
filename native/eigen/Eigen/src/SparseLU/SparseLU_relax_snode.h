#ifndef SPARSELU_RELAX_SNODE_H
#define SPARSELU_RELAX_SNODE_H

namespace Eigen {

namespace internal {

template <typename Scalar, typename StorageIndex>
void SparseLUImpl<Scalar,StorageIndex>::relax_snode (const Index n, IndexVector& et, const Index relax_columns, IndexVector& descendants, IndexVector& relax_end)
{

  Index parent; 
  relax_end.setConstant(emptyIdxLU);
  descendants.setZero();
  for (Index j = 0; j < n; j++) 
  {
    parent = et(j);
    if (parent != n) 
      descendants(parent) += descendants(j) + 1;
  }
  
  Index snode_start; 
  for (Index j = 0; j < n; )
  {
    parent = et(j);
    snode_start = j; 
    while ( parent != n && descendants(parent) < relax_columns ) 
    {
      j = parent; 
      parent = et(j);
    }
    
    relax_end(snode_start) = StorageIndex(j); 
    j++;
    
    while (descendants(j) != 0 && j < n) j++;
  } 
  
}

} 

} 
#endif
