#ifndef SPARSELU_PANEL_DFS_H
#define SPARSELU_PANEL_DFS_H

namespace Eigen {

namespace internal {
  
template<typename IndexVector>
struct panel_dfs_traits
{
  typedef typename IndexVector::Scalar StorageIndex;
  panel_dfs_traits(Index jcol, StorageIndex* marker)
    : m_jcol(jcol), m_marker(marker)
  {}
  bool update_segrep(Index krep, StorageIndex jj)
  {
    if(m_marker[krep]<m_jcol)
    {
      m_marker[krep] = jj; 
      return true;
    }
    return false;
  }
  void mem_expand(IndexVector& , Index , Index ) {}
  enum { ExpandMem = false };
  Index m_jcol;
  StorageIndex* m_marker;
};

template <typename Scalar, typename StorageIndex>
template <typename Traits>
void SparseLUImpl<Scalar,StorageIndex>::dfs_kernel(const StorageIndex jj, IndexVector& perm_r,
                   Index& nseg, IndexVector& panel_lsub, IndexVector& segrep,
                   Ref<IndexVector> repfnz_col, IndexVector& xprune, Ref<IndexVector> marker, IndexVector& parent,
                   IndexVector& xplore, GlobalLU_t& glu,
                   Index& nextl_col, Index krow, Traits& traits
                  )
{
  
  StorageIndex kmark = marker(krow);

  marker(krow) = jj; 
  StorageIndex kperm = perm_r(krow); 
  if (kperm == emptyIdxLU ) {
    
    panel_lsub(nextl_col++) = StorageIndex(krow);  
    
    traits.mem_expand(panel_lsub, nextl_col, kmark);
  }
  else 
  {

    StorageIndex krep = glu.xsup(glu.supno(kperm)+1) - 1; 
    
    StorageIndex myfnz = repfnz_col(krep); 
    
    if (myfnz != emptyIdxLU )
    {
      
      if (myfnz > kperm ) repfnz_col(krep) = kperm; 
      
    }
    else 
    {
      
      StorageIndex oldrep = emptyIdxLU; 
      parent(krep) = oldrep; 
      repfnz_col(krep) = kperm; 
      StorageIndex xdfs =  glu.xlsub(krep); 
      Index maxdfs = xprune(krep); 
      
      StorageIndex kpar;
      do 
      {
        
        while (xdfs < maxdfs) 
        {
          StorageIndex kchild = glu.lsub(xdfs); 
          xdfs++; 
          StorageIndex chmark = marker(kchild); 
          
          if (chmark != jj ) 
          {
            marker(kchild) = jj; 
            StorageIndex chperm = perm_r(kchild); 
            
            if (chperm == emptyIdxLU) 
            {
              
              panel_lsub(nextl_col++) = kchild;
              traits.mem_expand(panel_lsub, nextl_col, chmark);
            }
            else
            {

              StorageIndex chrep = glu.xsup(glu.supno(chperm)+1) - 1; 
              myfnz = repfnz_col(chrep); 
              
              if (myfnz != emptyIdxLU) 
              { 
                if (myfnz > chperm) 
                  repfnz_col(chrep) = chperm; 
              }
              else 
              { 
                xplore(krep) = xdfs; 
                oldrep = krep; 
                krep = chrep; 
                parent(krep) = oldrep; 
                repfnz_col(krep) = chperm; 
                xdfs = glu.xlsub(krep); 
                maxdfs = xprune(krep); 
                
              } 
            } 
                
          } 
        } 

        if(traits.update_segrep(krep,jj))
        
        {
          segrep(nseg) = krep; 
          ++nseg; 
          
        }
        
        kpar = parent(krep); 
        if (kpar == emptyIdxLU) 
          break; 
        krep = kpar; 
        xdfs = xplore(krep); 
        maxdfs = xprune(krep); 

      } while (kpar != emptyIdxLU); 
      
    } 

  } 
}

template <typename Scalar, typename StorageIndex>
void SparseLUImpl<Scalar,StorageIndex>::panel_dfs(const Index m, const Index w, const Index jcol, MatrixType& A, IndexVector& perm_r, Index& nseg, ScalarVector& dense, IndexVector& panel_lsub, IndexVector& segrep, IndexVector& repfnz, IndexVector& xprune, IndexVector& marker, IndexVector& parent, IndexVector& xplore, GlobalLU_t& glu)
{
  Index nextl_col; 

  VectorBlock<IndexVector> marker1(marker, m, m); 
  nseg = 0; 
  
  panel_dfs_traits<IndexVector> traits(jcol, marker1.data());

  for (StorageIndex jj = StorageIndex(jcol); jj < jcol + w; jj++) 
  {
    nextl_col = (jj - jcol) * m; 
    
    VectorBlock<IndexVector> repfnz_col(repfnz, nextl_col, m); 
    VectorBlock<ScalarVector> dense_col(dense,nextl_col, m); 

    for (typename MatrixType::InnerIterator it(A, jj); it; ++it)
    {
      Index krow = it.row(); 
      dense_col(krow) = it.value();
      
      StorageIndex kmark = marker(krow); 
      if (kmark == jj) 
        continue; 
      
      dfs_kernel(jj, perm_r, nseg, panel_lsub, segrep, repfnz_col, xprune, marker, parent,
                   xplore, glu, nextl_col, krow, traits);
    }
    
  } 
}

} 
} 

#endif 
