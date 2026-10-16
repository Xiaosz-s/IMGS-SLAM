#ifndef SPARSELU_COLUMN_DFS_H
#define SPARSELU_COLUMN_DFS_H

template <typename Scalar, typename StorageIndex> class SparseLUImpl;
namespace Eigen {

namespace internal {

template<typename IndexVector, typename ScalarVector>
struct column_dfs_traits : no_assignment_operator
{
  typedef typename ScalarVector::Scalar Scalar;
  typedef typename IndexVector::Scalar StorageIndex;
  column_dfs_traits(Index jcol, Index& jsuper, typename SparseLUImpl<Scalar, StorageIndex>::GlobalLU_t& glu, SparseLUImpl<Scalar, StorageIndex>& luImpl)
   : m_jcol(jcol), m_jsuper_ref(jsuper), m_glu(glu), m_luImpl(luImpl)
 {}
  bool update_segrep(Index , Index )
  {
    return true;
  }
  void mem_expand(IndexVector& lsub, Index& nextl, Index chmark)
  {
    if (nextl >= m_glu.nzlmax)
      m_luImpl.memXpand(lsub, m_glu.nzlmax, nextl, LSUB, m_glu.num_expansions); 
    if (chmark != (m_jcol-1)) m_jsuper_ref = emptyIdxLU;
  }
  enum { ExpandMem = true };
  
  Index m_jcol;
  Index& m_jsuper_ref;
  typename SparseLUImpl<Scalar, StorageIndex>::GlobalLU_t& m_glu;
  SparseLUImpl<Scalar, StorageIndex>& m_luImpl;
};

template <typename Scalar, typename StorageIndex>
Index SparseLUImpl<Scalar,StorageIndex>::column_dfs(const Index m, const Index jcol, IndexVector& perm_r, Index maxsuper, Index& nseg,
                                                    BlockIndexVector lsub_col, IndexVector& segrep, BlockIndexVector repfnz, IndexVector& xprune,
                                                    IndexVector& marker, IndexVector& parent, IndexVector& xplore, GlobalLU_t& glu)
{
  
  Index jsuper = glu.supno(jcol); 
  Index nextl = glu.xlsub(jcol); 
  VectorBlock<IndexVector> marker2(marker, 2*m, m); 

  column_dfs_traits<IndexVector, ScalarVector> traits(jcol, jsuper, glu, *this);

  for (Index k = 0; ((k < m) ? lsub_col[k] != emptyIdxLU : false) ; k++)
  {
    Index krow = lsub_col(k); 
    lsub_col(k) = emptyIdxLU; 
    Index kmark = marker2(krow); 

    if (kmark == jcol) continue;
    
    dfs_kernel(StorageIndex(jcol), perm_r, nseg, glu.lsub, segrep, repfnz, xprune, marker2, parent,
                   xplore, glu, nextl, krow, traits);
  } 
  
  Index fsupc;
  StorageIndex nsuper = glu.supno(jcol);
  StorageIndex jcolp1 = StorageIndex(jcol) + 1;
  Index jcolm1 = jcol - 1;

  if ( jcol == 0 )
  { 
    nsuper = glu.supno(0) = 0 ;
  }
  else 
  {
    fsupc = glu.xsup(nsuper); 
    StorageIndex jptr = glu.xlsub(jcol); 
    StorageIndex jm1ptr = glu.xlsub(jcolm1); 

    if ( (nextl-jptr != jptr-jm1ptr-1) ) jsuper = emptyIdxLU;

    if ( (jcol - fsupc) >= maxsuper) jsuper = emptyIdxLU; 

    if (jsuper == emptyIdxLU)
    { 
      if ( (fsupc < jcolm1-1) ) 
      { 
        StorageIndex ito = glu.xlsub(fsupc+1);
        glu.xlsub(jcolm1) = ito; 
        StorageIndex istop = ito + jptr - jm1ptr; 
        xprune(jcolm1) = istop; 
        glu.xlsub(jcol) = istop; 
        
        for (StorageIndex ifrom = jm1ptr; ifrom < nextl; ++ifrom, ++ito)
          glu.lsub(ito) = glu.lsub(ifrom); 
        nextl = ito;  
      }
      nsuper++; 
      glu.supno(jcol) = nsuper; 
    } 
  } 

  glu.xsup(nsuper+1) = jcolp1; 
  glu.supno(jcolp1) = nsuper; 
  xprune(jcol) = StorageIndex(nextl);  
  glu.xlsub(jcolp1) = StorageIndex(nextl); 
  
  return 0; 
}

} 

} 

#endif
