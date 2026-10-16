#ifndef EIGEN_SPARSELU_UTILS_H
#define EIGEN_SPARSELU_UTILS_H

namespace Eigen {
namespace internal {

template <typename Scalar, typename StorageIndex>
void SparseLUImpl<Scalar,StorageIndex>::countnz(const Index n, Index& nnzL, Index& nnzU, GlobalLU_t& glu)
{
 nnzL = 0; 
 nnzU = (glu.xusub)(n); 
 Index nsuper = (glu.supno)(n); 
 Index jlen; 
 Index i, j, fsupc;
 if (n <= 0 ) return; 
 
 for (i = 0; i <= nsuper; i++)
 {
   fsupc = glu.xsup(i); 
   jlen = glu.xlsub(fsupc+1) - glu.xlsub(fsupc); 
   
   for (j = fsupc; j < glu.xsup(i+1); j++)
   {
     nnzL += jlen; 
     nnzU += j - fsupc + 1; 
     jlen--; 
   }
 }
}

template <typename Scalar, typename StorageIndex>
void SparseLUImpl<Scalar,StorageIndex>::fixupL(const Index n, const IndexVector& perm_r, GlobalLU_t& glu)
{
  Index fsupc, i, j, k, jstart; 
  
  StorageIndex nextl = 0; 
  Index nsuper = (glu.supno)(n); 

  for (i = 0; i <= nsuper; i++)
  {
    fsupc = glu.xsup(i); 
    jstart = glu.xlsub(fsupc); 
    glu.xlsub(fsupc) = nextl; 
    for (j = jstart; j < glu.xlsub(fsupc + 1); j++)
    {
      glu.lsub(nextl) = perm_r(glu.lsub(j)); 
      nextl++;
    }
    for (k = fsupc+1; k < glu.xsup(i+1); k++)
      glu.xlsub(k) = nextl; 
  }
  
  glu.xlsub(n) = nextl; 
}

} 

} 
#endif 
