#ifndef SPARSELU_PRUNEL_H
#define SPARSELU_PRUNEL_H

namespace Eigen {
namespace internal {

template <typename Scalar, typename StorageIndex>
void SparseLUImpl<Scalar,StorageIndex>::pruneL(const Index jcol, const IndexVector& perm_r, const Index pivrow, const Index nseg,
                                               const IndexVector& segrep, BlockIndexVector repfnz, IndexVector& xprune, GlobalLU_t& glu)
{
  
  Index jsupno = glu.supno(jcol); 
  Index i,irep,irep1; 
  bool movnum, do_prune = false; 
  Index kmin = 0, kmax = 0, minloc, maxloc,krow; 
  for (i = 0; i < nseg; i++)
  {
    irep = segrep(i); 
    irep1 = irep + 1; 
    do_prune = false; 

    if (repfnz(irep) == emptyIdxLU) continue; 

    if (glu.supno(irep) == glu.supno(irep1) ) continue; 

    if (glu.supno(irep) != jsupno )
    {
      if ( xprune (irep) >= glu.xlsub(irep1) )
      {
        kmin = glu.xlsub(irep);
        kmax = glu.xlsub(irep1) - 1; 
        for (krow = kmin; krow <= kmax; krow++)
        {
          if (glu.lsub(krow) == pivrow) 
          {
            do_prune = true; 
            break; 
          }
        }
      }
      
      if (do_prune) 
      {

        movnum = false; 
        if (irep == glu.xsup(glu.supno(irep)) ) 
          movnum = true; 
        
        while (kmin <= kmax)
        {
          if (perm_r(glu.lsub(kmax)) == emptyIdxLU)
            kmax--; 
          else if ( perm_r(glu.lsub(kmin)) != emptyIdxLU)
            kmin++;
          else 
          {

            std::swap(glu.lsub(kmin), glu.lsub(kmax)); 

            if (movnum) 
            {
              minloc = glu.xlusup(irep) + ( kmin - glu.xlsub(irep) ); 
              maxloc = glu.xlusup(irep) + ( kmax - glu.xlsub(irep) ); 
              std::swap(glu.lusup(minloc), glu.lusup(maxloc)); 
            }
            kmin++;
            kmax--;
          }
        } 
        
        xprune(irep) = StorageIndex(kmin);  
      } 
    } 
  } 
}

} 
} 

#endif 
