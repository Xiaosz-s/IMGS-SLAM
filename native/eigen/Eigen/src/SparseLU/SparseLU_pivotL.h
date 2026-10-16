#ifndef SPARSELU_PIVOTL_H
#define SPARSELU_PIVOTL_H

namespace Eigen {
namespace internal {

template <typename Scalar, typename StorageIndex>
Index SparseLUImpl<Scalar,StorageIndex>::pivotL(const Index jcol, const RealScalar& diagpivotthresh, IndexVector& perm_r, IndexVector& iperm_c, Index& pivrow, GlobalLU_t& glu)
{
  
  Index fsupc = (glu.xsup)((glu.supno)(jcol)); 
  Index nsupc = jcol - fsupc; 
  Index lptr = glu.xlsub(fsupc); 
  Index nsupr = glu.xlsub(fsupc+1) - lptr; 
  Index lda = glu.xlusup(fsupc+1) - glu.xlusup(fsupc); 
  Scalar* lu_sup_ptr = &(glu.lusup.data()[glu.xlusup(fsupc)]); 
  Scalar* lu_col_ptr = &(glu.lusup.data()[glu.xlusup(jcol)]); 
  StorageIndex* lsub_ptr = &(glu.lsub.data()[lptr]); 

  Index diagind = iperm_c(jcol); 
  RealScalar pivmax(-1.0);
  Index pivptr = nsupc; 
  Index diag = emptyIdxLU; 
  RealScalar rtemp;
  Index isub, icol, itemp, k; 
  for (isub = nsupc; isub < nsupr; ++isub) {
    using std::abs;
    rtemp = abs(lu_col_ptr[isub]);
    if (rtemp > pivmax) {
      pivmax = rtemp; 
      pivptr = isub;
    } 
    if (lsub_ptr[isub] == diagind) diag = isub;
  }

  if ( pivmax <= RealScalar(0.0) ) {
    
    pivrow = pivmax < RealScalar(0.0) ? diagind : lsub_ptr[pivptr];
    perm_r(pivrow) = StorageIndex(jcol);
    return (jcol+1);
  }
  
  RealScalar thresh = diagpivotthresh * pivmax; 

  {
    
    if (diag >= 0 ) 
    {
      
      using std::abs;
      rtemp = abs(lu_col_ptr[diag]);
      if (rtemp != RealScalar(0.0) && rtemp >= thresh) pivptr = diag;
    }
    pivrow = lsub_ptr[pivptr];
  }

  perm_r(pivrow) = StorageIndex(jcol);
  
  if (pivptr != nsupc )
  {
    std::swap( lsub_ptr[pivptr], lsub_ptr[nsupc] );

    for (icol = 0; icol <= nsupc; icol++)
    {
      itemp = pivptr + icol * lda; 
      std::swap(lu_sup_ptr[itemp], lu_sup_ptr[nsupc + icol * lda]);
    }
  }
  
  Scalar temp = Scalar(1.0) / lu_col_ptr[nsupc];
  for (k = nsupc+1; k < nsupr; k++)
    lu_col_ptr[k] *= temp; 
  return 0;
}

} 
} 

#endif 
