#ifndef METIS_SUPPORT_H
#define METIS_SUPPORT_H

namespace Eigen {

template <typename StorageIndex>
class MetisOrdering
{
public:
  typedef PermutationMatrix<Dynamic,Dynamic,StorageIndex> PermutationType;
  typedef Matrix<StorageIndex,Dynamic,1> IndexVector; 
  
  template <typename MatrixType>
  void get_symmetrized_graph(const MatrixType& A)
  {
    Index m = A.cols(); 
    eigen_assert((A.rows() == A.cols()) && "ONLY FOR SQUARED MATRICES");
    
    MatrixType At = A.transpose(); 
    
    Index TotNz = 0; 
    IndexVector visited(m); 
    visited.setConstant(-1); 
    for (StorageIndex j = 0; j < m; j++)
    {
      
      visited(j) = j; 
      
      for (typename MatrixType::InnerIterator it(A, j); it; ++it)
      {
        Index idx = it.index(); 
        if (visited(idx) != j ) 
        {
          visited(idx) = j; 
          ++TotNz; 
        }
      }
      
      for (typename MatrixType::InnerIterator it(At, j); it; ++it)
      {
        Index idx = it.index(); 
        if(visited(idx) != j)
        {
          visited(idx) = j; 
          ++TotNz; 
        }
      }
    }
    
    m_indexPtr.resize(m+1);
    m_innerIndices.resize(TotNz); 

    visited.setConstant(-1); 
    StorageIndex CurNz = 0; 
    for (StorageIndex j = 0; j < m; j++)
    {
      m_indexPtr(j) = CurNz; 
      
      visited(j) = j; 
      
      for (typename MatrixType::InnerIterator it(A,j); it; ++it)
      {
        StorageIndex idx = it.index(); 
        if (visited(idx) != j ) 
        {
          visited(idx) = j; 
          m_innerIndices(CurNz) = idx; 
          CurNz++; 
        }
      }
      
      for (typename MatrixType::InnerIterator it(At, j); it; ++it)
      {
        StorageIndex idx = it.index(); 
        if(visited(idx) != j)
        {
          visited(idx) = j; 
          m_innerIndices(CurNz) = idx; 
          ++CurNz; 
        }
      }
    }
    m_indexPtr(m) = CurNz;    
  }
  
  template <typename MatrixType>
  void operator() (const MatrixType& A, PermutationType& matperm)
  {
     StorageIndex m = internal::convert_index<StorageIndex>(A.cols()); 
     IndexVector perm(m),iperm(m); 
    
     get_symmetrized_graph(A); 
     int output_error;

     output_error = METIS_NodeND(&m, m_indexPtr.data(), m_innerIndices.data(), NULL, NULL, perm.data(), iperm.data());
     
    if(output_error != METIS_OK) 
    {
      
     std::cerr << "ERROR WHILE CALLING THE METIS PACKAGE \n"; 
     return; 
    }

     matperm.resize(m);
     for (int j = 0; j < m; j++)
       matperm.indices()(iperm(j)) = j;
   
  }
  
  protected:
    IndexVector m_indexPtr; 
    IndexVector m_innerIndices; 
};

}
#endif
