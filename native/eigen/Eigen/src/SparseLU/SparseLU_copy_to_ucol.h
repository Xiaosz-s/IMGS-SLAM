#ifndef SPARSELU_COPY_TO_UCOL_H
#define SPARSELU_COPY_TO_UCOL_H

namespace Eigen {
namespace internal {

template <typename Scalar, typename StorageIndex>
Index SparseLUImpl<Scalar,StorageIndex>::copy_to_ucol(const Index jcol, const Index nseg, IndexVector& segrep,
                                                      BlockIndexVector repfnz ,IndexVector& perm_r, BlockScalarVector dense, GlobalLU_t& glu)
{  
  Index ksub, krep, ksupno; 
    
  Index jsupno = glu.supno(jcol);

  Index k = nseg - 1, i; 
  StorageIndex nextu = glu.xusub(jcol); 
  Index kfnz, isub, segsize; 
  Index new_next,irow; 
  Index fsupc, mem; 
  for (ksub = 0; ksub < nseg; ksub++)
  {
    krep = segrep(k); k--; 
    ksupno = glu.supno(krep); 
    if (jsupno != ksupno ) 
    {
      kfnz = repfnz(krep); 
      if (kfnz != emptyIdxLU)
      { 
        fsupc = glu.xsup(ksupno); 
        isub = glu.xlsub(fsupc) + kfnz - fsupc; 
        segsize = krep - kfnz + 1; 
        new_next = nextu + segsize; 
        while (new_next > glu.nzumax) 
        {
          mem = memXpand<ScalarVector>(glu.ucol, glu.nzumax, nextu, UCOL, glu.num_expansions); 
          if (mem) return mem; 
          mem = memXpand<IndexVector>(glu.usub, glu.nzumax, nextu, USUB, glu.num_expansions); 
          if (mem) return mem; 
          
        }
        
        for (i = 0; i < segsize; i++)
        {
          irow = glu.lsub(isub); 
          glu.usub(nextu) = perm_r(irow); 
          glu.ucol(nextu) = dense(irow); 
          dense(irow) = Scalar(0.0); 
          nextu++;
          isub++;
        }
        
      } 
      
    } 
    
  } 
  glu.xusub(jcol + 1) = nextu; 
  return 0; 
}

} 
} 

#endif 
