#ifndef EIGEN_INCOMPLETE_CHOlESKY_H
#define EIGEN_INCOMPLETE_CHOlESKY_H

#include <vector>
#include <list>

namespace Eigen {

template <typename Scalar, int UpLo_ = Lower, typename OrderingType_ = AMDOrdering<int> >
class IncompleteCholesky : public SparseSolverBase<IncompleteCholesky<Scalar,UpLo_,OrderingType_> >
{
  protected:
    typedef SparseSolverBase<IncompleteCholesky<Scalar,UpLo_,OrderingType_> > Base;
    using Base::m_isInitialized;
  public:
    typedef typename NumTraits<Scalar>::Real RealScalar;
    typedef OrderingType_ OrderingType;
    typedef typename OrderingType::PermutationType PermutationType;
    typedef typename PermutationType::StorageIndex StorageIndex;
    typedef SparseMatrix<Scalar,ColMajor,StorageIndex> FactorType;
    typedef Matrix<Scalar,Dynamic,1> VectorSx;
    typedef Matrix<RealScalar,Dynamic,1> VectorRx;
    typedef Matrix<StorageIndex,Dynamic, 1> VectorIx;
    typedef std::vector<std::list<StorageIndex> > VectorList;
    enum { UpLo = UpLo_ };
    enum {
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };
  public:

    IncompleteCholesky() : m_initialShift(1e-3),m_analysisIsOk(false),m_factorizationIsOk(false) {}

    template<typename MatrixType>
    IncompleteCholesky(const MatrixType& matrix) : m_initialShift(1e-3),m_analysisIsOk(false),m_factorizationIsOk(false)
    {
      compute(matrix);
    }

    EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return m_L.rows(); }

    EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return m_L.cols(); }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "IncompleteCholesky is not initialized.");
      return m_info;
    }

    void setInitialShift(RealScalar shift) { m_initialShift = shift; }

    template<typename MatrixType>
    void analyzePattern(const MatrixType& mat)
    {
      OrderingType ord;
      PermutationType pinv;
      ord(mat.template selfadjointView<UpLo>(), pinv);
      if(pinv.size()>0) m_perm = pinv.inverse();
      else              m_perm.resize(0);
      m_L.resize(mat.rows(), mat.cols());
      m_analysisIsOk = true;
      m_isInitialized = true;
      m_info = Success;
    }

    template<typename MatrixType>
    void factorize(const MatrixType& mat);

    template<typename MatrixType>
    void compute(const MatrixType& mat)
    {
      analyzePattern(mat);
      factorize(mat);
    }

    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
      eigen_assert(m_factorizationIsOk && "factorize() should be called first");
      if (m_perm.rows() == b.rows())  x = m_perm * b;
      else                            x = b;
      x = m_scale.asDiagonal() * x;
      x = m_L.template triangularView<Lower>().solve(x);
      x = m_L.adjoint().template triangularView<Upper>().solve(x);
      x = m_scale.asDiagonal() * x;
      if (m_perm.rows() == b.rows())
        x = m_perm.inverse() * x;
    }

    const FactorType& matrixL() const { eigen_assert("m_factorizationIsOk"); return m_L; }

    const VectorRx& scalingS() const { eigen_assert("m_factorizationIsOk"); return m_scale; }

    const PermutationType& permutationP() const { eigen_assert("m_analysisIsOk"); return m_perm; }

  protected:
    FactorType m_L;              
    VectorRx m_scale;            
    RealScalar m_initialShift;   
    bool m_analysisIsOk;
    bool m_factorizationIsOk;
    ComputationInfo m_info;
    PermutationType m_perm;

  private:
    inline void updateList(Ref<const VectorIx> colPtr, Ref<VectorIx> rowIdx, Ref<VectorSx> vals, const Index& col, const Index& jk, VectorIx& firstElt, VectorList& listCol);
};

template<typename Scalar, int UpLo_, typename OrderingType>
template<typename MatrixType_>
void IncompleteCholesky<Scalar,UpLo_, OrderingType>::factorize(const MatrixType_& mat)
{
  using std::sqrt;
  eigen_assert(m_analysisIsOk && "analyzePattern() should be called first");

  if (m_perm.rows() == mat.rows() ) 
  {
    
    FactorType tmp(mat.rows(), mat.cols());
    tmp = mat.template selfadjointView<UpLo_>().twistedBy(m_perm);
    m_L.template selfadjointView<Lower>() = tmp.template selfadjointView<Lower>();
  }
  else
  {
    m_L.template selfadjointView<Lower>() = mat.template selfadjointView<UpLo_>();
  }

  Index n = m_L.cols();
  Index nnz = m_L.nonZeros();
  Map<VectorSx> vals(m_L.valuePtr(), nnz);         
  Map<VectorIx> rowIdx(m_L.innerIndexPtr(), nnz);  
  Map<VectorIx> colPtr( m_L.outerIndexPtr(), n+1); 
  VectorIx firstElt(n-1); 
  VectorList listCol(n);  
  VectorSx col_vals(n);   
  VectorIx col_irow(n);   
  VectorIx col_pattern(n);
  col_pattern.fill(-1);
  StorageIndex col_nnz;

  m_scale.resize(n);
  m_scale.setZero();
  for (Index j = 0; j < n; j++)
    for (Index k = colPtr[j]; k < colPtr[j+1]; k++)
    {
      m_scale(j) += numext::abs2(vals(k));
      if(rowIdx[k]!=j)
        m_scale(rowIdx[k]) += numext::abs2(vals(k));
    }

  m_scale = m_scale.cwiseSqrt().cwiseSqrt();

  for (Index j = 0; j < n; ++j)
    if(m_scale(j)>(std::numeric_limits<RealScalar>::min)())
      m_scale(j) = RealScalar(1)/m_scale(j);
    else
      m_scale(j) = 1;

  RealScalar mindiag = NumTraits<RealScalar>::highest();
  for (Index j = 0; j < n; j++)
  {
    for (Index k = colPtr[j]; k < colPtr[j+1]; k++)
      vals[k] *= (m_scale(j)*m_scale(rowIdx[k]));
    eigen_internal_assert(rowIdx[colPtr[j]]==j && "IncompleteCholesky: only the lower triangular part must be stored");
    mindiag = numext::mini(numext::real(vals[colPtr[j]]), mindiag);
  }

  FactorType L_save = m_L;

  RealScalar shift = 0;
  if(mindiag <= RealScalar(0.))
    shift = m_initialShift - mindiag;

  m_info = NumericalIssue;

  int iter = 0;
  do
  {
    
    for (Index j = 0; j < n; j++)
      vals[colPtr[j]] += shift;

    Index j=0;
    for (; j < n; ++j)
    {

      Scalar diag = vals[colPtr[j]];  
      col_nnz = 0;
      for (Index i = colPtr[j] + 1; i < colPtr[j+1]; i++)
      {
        StorageIndex l = rowIdx[i];
        col_vals(col_nnz) = vals[i];
        col_irow(col_nnz) = l;
        col_pattern(l) = col_nnz;
        col_nnz++;
      }
      {
        typename std::list<StorageIndex>::iterator k;
        
        for(k = listCol[j].begin(); k != listCol[j].end(); k++)
        {
          Index jk = firstElt(*k); 
          eigen_internal_assert(rowIdx[jk]==j);
          Scalar v_j_jk = numext::conj(vals[jk]);

          jk += 1;
          for (Index i = jk; i < colPtr[*k+1]; i++)
          {
            StorageIndex l = rowIdx[i];
            if(col_pattern[l]<0)
            {
              col_vals(col_nnz) = vals[i] * v_j_jk;
              col_irow[col_nnz] = l;
              col_pattern(l) = col_nnz;
              col_nnz++;
            }
            else
              col_vals(col_pattern[l]) -= vals[i] * v_j_jk;
          }
          updateList(colPtr,rowIdx,vals, *k, jk, firstElt, listCol);
        }
      }

      if(numext::real(diag) <= 0)
      {
        if(++iter>=10)
          return;

        shift = numext::maxi(m_initialShift,RealScalar(2)*shift);
        
        vals = Map<const VectorSx>(L_save.valuePtr(), nnz);
        rowIdx = Map<const VectorIx>(L_save.innerIndexPtr(), nnz);
        colPtr = Map<const VectorIx>(L_save.outerIndexPtr(), n+1);
        col_pattern.fill(-1);
        for(Index i=0; i<n; ++i)
          listCol[i].clear();

        break;
      }

      RealScalar rdiag = sqrt(numext::real(diag));
      vals[colPtr[j]] = rdiag;
      for (Index k = 0; k<col_nnz; ++k)
      {
        Index i = col_irow[k];
        
        col_vals(k) /= rdiag;
        
        vals[colPtr[i]] -= numext::abs2(col_vals(k));
      }

      Index p = colPtr[j+1] - colPtr[j] - 1 ;
      Ref<VectorSx> cvals = col_vals.head(col_nnz);
      Ref<VectorIx> cirow = col_irow.head(col_nnz);
      internal::QuickSplit(cvals,cirow, p);
      
      Index cpt = 0;
      for (Index i = colPtr[j]+1; i < colPtr[j+1]; i++)
      {
        vals[i] = col_vals(cpt);
        rowIdx[i] = col_irow(cpt);
        
        col_pattern(col_irow(cpt)) = -1;
        cpt++;
      }
      
      Index jk = colPtr(j)+1;
      updateList(colPtr,rowIdx,vals,j,jk,firstElt,listCol);
    }

    if(j==n)
    {
      m_factorizationIsOk = true;
      m_info = Success;
    }
  } while(m_info!=Success);
}

template<typename Scalar, int UpLo_, typename OrderingType>
inline void IncompleteCholesky<Scalar,UpLo_, OrderingType>::updateList(Ref<const VectorIx> colPtr, Ref<VectorIx> rowIdx, Ref<VectorSx> vals, const Index& col, const Index& jk, VectorIx& firstElt, VectorList& listCol)
{
  if (jk < colPtr(col+1) )
  {
    Index p = colPtr(col+1) - jk;
    Index minpos;
    rowIdx.segment(jk,p).minCoeff(&minpos);
    minpos += jk;
    if (rowIdx(minpos) != rowIdx(jk))
    {
      
      std::swap(rowIdx(jk),rowIdx(minpos));
      std::swap(vals(jk),vals(minpos));
    }
    firstElt(col) = internal::convert_index<StorageIndex,Index>(jk);
    listCol[rowIdx(jk)].push_back(internal::convert_index<StorageIndex,Index>(col));
  }
}

} 

#endif
