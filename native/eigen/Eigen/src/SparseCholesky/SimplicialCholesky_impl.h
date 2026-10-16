#ifndef EIGEN_SIMPLICIAL_CHOLESKY_IMPL_H
#define EIGEN_SIMPLICIAL_CHOLESKY_IMPL_H

namespace Eigen {

template<typename Derived>
void SimplicialCholeskyBase<Derived>::analyzePattern_preordered(const CholMatrixType& ap, bool doLDLT)
{
  const StorageIndex size = StorageIndex(ap.rows());
  m_matrix.resize(size, size);
  m_parent.resize(size);
  m_nonZerosPerCol.resize(size);

  ei_declare_aligned_stack_constructed_variable(StorageIndex, tags, size, 0);

  for(StorageIndex k = 0; k < size; ++k)
  {
    
    m_parent[k] = -1;             
    tags[k] = k;                  
    m_nonZerosPerCol[k] = 0;      
    for(typename CholMatrixType::InnerIterator it(ap,k); it; ++it)
    {
      StorageIndex i = it.index();
      if(i < k)
      {
        
        for(; tags[i] != k; i = m_parent[i])
        {
          
          if (m_parent[i] == -1)
            m_parent[i] = k;
          m_nonZerosPerCol[i]++;        
          tags[i] = k;                  
        }
      }
    }
  }

  StorageIndex* Lp = m_matrix.outerIndexPtr();
  Lp[0] = 0;
  for(StorageIndex k = 0; k < size; ++k)
    Lp[k+1] = Lp[k] + m_nonZerosPerCol[k] + (doLDLT ? 0 : 1);

  m_matrix.resizeNonZeros(Lp[size]);

  m_isInitialized     = true;
  m_info              = Success;
  m_analysisIsOk      = true;
  m_factorizationIsOk = false;
}

template<typename Derived>
template<bool DoLDLT>
void SimplicialCholeskyBase<Derived>::factorize_preordered(const CholMatrixType& ap)
{
  using std::sqrt;

  eigen_assert(m_analysisIsOk && "You must first call analyzePattern()");
  eigen_assert(ap.rows()==ap.cols());
  eigen_assert(m_parent.size()==ap.rows());
  eigen_assert(m_nonZerosPerCol.size()==ap.rows());

  const StorageIndex size = StorageIndex(ap.rows());
  const StorageIndex* Lp = m_matrix.outerIndexPtr();
  StorageIndex* Li = m_matrix.innerIndexPtr();
  Scalar* Lx = m_matrix.valuePtr();

  ei_declare_aligned_stack_constructed_variable(Scalar, y, size, 0);
  ei_declare_aligned_stack_constructed_variable(StorageIndex,  pattern, size, 0);
  ei_declare_aligned_stack_constructed_variable(StorageIndex,  tags, size, 0);

  bool ok = true;
  m_diag.resize(DoLDLT ? size : 0);

  for(StorageIndex k = 0; k < size; ++k)
  {
    
    y[k] = Scalar(0);                     
    StorageIndex top = size;               
    tags[k] = k;                    
    m_nonZerosPerCol[k] = 0;        
    for(typename CholMatrixType::InnerIterator it(ap,k); it; ++it)
    {
      StorageIndex i = it.index();
      if(i <= k)
      {
        y[i] += numext::conj(it.value());            
        Index len;
        for(len = 0; tags[i] != k; i = m_parent[i])
        {
          pattern[len++] = i;     
          tags[i] = k;            
        }
        while(len > 0)
          pattern[--top] = pattern[--len];
      }
    }

    RealScalar d = numext::real(y[k]) * m_shiftScale + m_shiftOffset;    
    y[k] = Scalar(0);
    for(; top < size; ++top)
    {
      Index i = pattern[top];       
      Scalar yi = y[i];             
      y[i] = Scalar(0);

      Scalar l_ki;
      if(DoLDLT)
        l_ki = yi / numext::real(m_diag[i]);
      else
        yi = l_ki = yi / Lx[Lp[i]];

      Index p2 = Lp[i] + m_nonZerosPerCol[i];
      Index p;
      for(p = Lp[i] + (DoLDLT ? 0 : 1); p < p2; ++p)
        y[Li[p]] -= numext::conj(Lx[p]) * yi;
      d -= numext::real(l_ki * numext::conj(yi));
      Li[p] = k;                          
      Lx[p] = l_ki;
      ++m_nonZerosPerCol[i];              
    }
    if(DoLDLT)
    {
      m_diag[k] = d;
      if(d == RealScalar(0))
      {
        ok = false;                         
        break;
      }
    }
    else
    {
      Index p = Lp[k] + m_nonZerosPerCol[k]++;
      Li[p] = k ;                
      if(d <= RealScalar(0)) {
        ok = false;              
        break;
      }
      Lx[p] = sqrt(d) ;
    }
  }

  m_info = ok ? Success : NumericalIssue;
  m_factorizationIsOk = true;
}

} 

#endif 
