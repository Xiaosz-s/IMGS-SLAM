#ifndef SPARSE_COLETREE_H
#define SPARSE_COLETREE_H

namespace Eigen {

namespace internal {

template<typename Index, typename IndexVector>
Index etree_find (Index i, IndexVector& pp)
{
  Index p = pp(i); 
  Index gp = pp(p); 
  while (gp != p) 
  {
    pp(i) = gp; 
    i = gp; 
    p = pp(i);
    gp = pp(p);
  }
  return p; 
}

template <typename MatrixType, typename IndexVector>
int coletree(const MatrixType& mat, IndexVector& parent, IndexVector& firstRowElt, typename MatrixType::StorageIndex *perm=0)
{
  typedef typename MatrixType::StorageIndex StorageIndex;
  StorageIndex nc = convert_index<StorageIndex>(mat.cols()); 
  StorageIndex m = convert_index<StorageIndex>(mat.rows());
  StorageIndex diagSize = (std::min)(nc,m);
  IndexVector root(nc); 
  root.setZero();
  IndexVector pp(nc); 
  pp.setZero(); 
  parent.resize(mat.cols());
  
  firstRowElt.resize(m);
  firstRowElt.setConstant(nc);
  firstRowElt.segment(0, diagSize).setLinSpaced(diagSize, 0, diagSize-1);
  bool found_diag;
  for (StorageIndex col = 0; col < nc; col++)
  {
    StorageIndex pcol = col;
    if(perm) pcol  = perm[col];
    for (typename MatrixType::InnerIterator it(mat, pcol); it; ++it)
    { 
      Index row = it.row();
      firstRowElt(row) = (std::min)(firstRowElt(row), col);
    }
  }

  StorageIndex rset, cset, rroot;
  for (StorageIndex col = 0; col < nc; col++) 
  {
    found_diag = col>=m;
    pp(col) = col; 
    cset = col; 
    root(cset) = col; 
    parent(col) = nc; 

    StorageIndex pcol = col;
    if(perm) pcol  = perm[col];
    for (typename MatrixType::InnerIterator it(mat, pcol); it||!found_diag; ++it)
    { 
      Index i = col;
      if(it) i = it.index();
      if (i == col) found_diag = true;
      
      StorageIndex row = firstRowElt(i);
      if (row >= col) continue; 
      rset = internal::etree_find(row, pp); 
      rroot = root(rset);
      if (rroot != col) 
      {
        parent(rroot) = col; 
        pp(cset) = rset; 
        cset = rset; 
        root(cset) = col; 
      }
    }
  }
  return 0;  
}

template <typename IndexVector>
void nr_etdfs (typename IndexVector::Scalar n, IndexVector& parent, IndexVector& first_kid, IndexVector& next_kid, IndexVector& post, typename IndexVector::Scalar postnum)
{
  typedef typename IndexVector::Scalar StorageIndex;
  StorageIndex current = n, first, next;
  while (postnum != n) 
  {
    
    first = first_kid(current);

    if (first == -1) 
    {
      
      post(current) = postnum++;

      next = next_kid(current); 
      while (next == -1) 
      {
        
        current = parent(current); 
        
        post(current) = postnum++;

        next = next_kid(current); 
      }
      
      if (postnum == n+1) return; 

      current = next; 
    }
    else 
    {
      current = first; 
    }
  }
}

template <typename IndexVector>
void treePostorder(typename IndexVector::Scalar n, IndexVector& parent, IndexVector& post)
{
  typedef typename IndexVector::Scalar StorageIndex;
  IndexVector first_kid, next_kid; 
  StorageIndex postnum; 
  
  first_kid.resize(n+1); 
  next_kid.setZero(n+1);
  post.setZero(n+1);

  first_kid.setConstant(-1); 
  for (StorageIndex v = n-1; v >= 0; v--) 
  {
    StorageIndex dad = parent(v);
    next_kid(v) = first_kid(dad); 
    first_kid(dad) = v; 
  }

  postnum = 0; 
  internal::nr_etdfs(n, parent, first_kid, next_kid, post, postnum);
}

} 

} 

#endif 
