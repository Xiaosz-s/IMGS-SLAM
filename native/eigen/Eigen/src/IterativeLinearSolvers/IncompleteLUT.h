#ifndef EIGEN_INCOMPLETE_LUT_H
#define EIGEN_INCOMPLETE_LUT_H

namespace Eigen {

namespace internal {

template <typename VectorV, typename VectorI>
Index QuickSplit(VectorV &row, VectorI &ind, Index ncut)
{
  typedef typename VectorV::RealScalar RealScalar;
  using std::swap;
  using std::abs;
  Index mid;
  Index n = row.size(); 
  Index first, last ;

  ncut--; 
  first = 0;
  last = n-1;
  if (ncut < first || ncut > last ) return 0;

  do {
    mid = first;
    RealScalar abskey = abs(row(mid));
    for (Index j = first + 1; j <= last; j++) {
      if ( abs(row(j)) > abskey) {
        ++mid;
        swap(row(mid), row(j));
        swap(ind(mid), ind(j));
      }
    }
    
    swap(row(mid), row(first));
    swap(ind(mid), ind(first));

    if (mid > ncut) last = mid - 1;
    else if (mid < ncut ) first = mid + 1;
  } while (mid != ncut );

  return 0; 
}

}

template <typename Scalar_, typename StorageIndex_ = int>
class IncompleteLUT : public SparseSolverBase<IncompleteLUT<Scalar_, StorageIndex_> >
{
  protected:
    typedef SparseSolverBase<IncompleteLUT> Base;
    using Base::m_isInitialized;
  public:
    typedef Scalar_ Scalar;
    typedef StorageIndex_ StorageIndex;
    typedef typename NumTraits<Scalar>::Real RealScalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    typedef Matrix<StorageIndex,Dynamic,1> VectorI;
    typedef SparseMatrix<Scalar,RowMajor,StorageIndex> FactorType;

    enum {
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };

  public:

    IncompleteLUT()
      : m_droptol(NumTraits<Scalar>::dummy_precision()), m_fillfactor(10),
        m_analysisIsOk(false), m_factorizationIsOk(false)
    {}

    template<typename MatrixType>
    explicit IncompleteLUT(const MatrixType& mat, const RealScalar& droptol=NumTraits<Scalar>::dummy_precision(), int fillfactor = 10)
      : m_droptol(droptol),m_fillfactor(fillfactor),
        m_analysisIsOk(false),m_factorizationIsOk(false)
    {
      eigen_assert(fillfactor != 0);
      compute(mat);
    }

    EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return m_lu.rows(); }

    EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return m_lu.cols(); }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "IncompleteLUT is not initialized.");
      return m_info;
    }

    template<typename MatrixType>
    void analyzePattern(const MatrixType& amat);

    template<typename MatrixType>
    void factorize(const MatrixType& amat);

    template<typename MatrixType>
    IncompleteLUT& compute(const MatrixType& amat)
    {
      analyzePattern(amat);
      factorize(amat);
      return *this;
    }

    void setDroptol(const RealScalar& droptol);
    void setFillfactor(int fillfactor);

    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
      x = m_Pinv * b;
      x = m_lu.template triangularView<UnitLower>().solve(x);
      x = m_lu.template triangularView<Upper>().solve(x);
      x = m_P * x;
    }

protected:

    struct keep_diag {
      inline bool operator() (const Index& row, const Index& col, const Scalar&) const
      {
        return row!=col;
      }
    };

protected:

    FactorType m_lu;
    RealScalar m_droptol;
    int m_fillfactor;
    bool m_analysisIsOk;
    bool m_factorizationIsOk;
    ComputationInfo m_info;
    PermutationMatrix<Dynamic,Dynamic,StorageIndex> m_P;     
    PermutationMatrix<Dynamic,Dynamic,StorageIndex> m_Pinv;  
};

template<typename Scalar, typename StorageIndex>
void IncompleteLUT<Scalar,StorageIndex>::setDroptol(const RealScalar& droptol)
{
  this->m_droptol = droptol;
}

template<typename Scalar, typename StorageIndex>
void IncompleteLUT<Scalar,StorageIndex>::setFillfactor(int fillfactor)
{
  this->m_fillfactor = fillfactor;
}

template <typename Scalar, typename StorageIndex>
template<typename MatrixType_>
void IncompleteLUT<Scalar,StorageIndex>::analyzePattern(const MatrixType_& amat)
{

  SparseMatrix<Scalar,ColMajor, StorageIndex> mat1 = amat;
  SparseMatrix<Scalar,ColMajor, StorageIndex> mat2 = amat.transpose();

  SparseMatrix<Scalar,ColMajor, StorageIndex> AtA = mat2 + mat1;
  AMDOrdering<StorageIndex> ordering;
  ordering(AtA,m_P);
  m_Pinv  = m_P.inverse(); 
  m_analysisIsOk = true;
  m_factorizationIsOk = false;
  m_isInitialized = true;
}

template <typename Scalar, typename StorageIndex>
template<typename MatrixType_>
void IncompleteLUT<Scalar,StorageIndex>::factorize(const MatrixType_& amat)
{
  using std::sqrt;
  using std::swap;
  using std::abs;
  using internal::convert_index;

  eigen_assert((amat.rows() == amat.cols()) && "The factorization should be done on a square matrix");
  Index n = amat.cols();  
  m_lu.resize(n,n);
  
  Vector u(n) ;     
  VectorI ju(n);   
  VectorI jr(n);   

  eigen_assert(m_analysisIsOk && "You must first call analyzePattern()");
  SparseMatrix<Scalar,RowMajor, StorageIndex> mat;
  mat = amat.twistedBy(m_Pinv);

  jr.fill(-1);
  ju.fill(0);
  u.fill(0);

  Index fill_in = (amat.nonZeros()*m_fillfactor)/n + 1;
  if (fill_in > n) fill_in = n;

  Index nnzL = fill_in/2;
  Index nnzU = nnzL;
  m_lu.reserve(n * (nnzL + nnzU + 1));

  for (Index ii = 0; ii < n; ii++)
  {

    Index sizeu = 1; 
    Index sizel = 0; 
    ju(ii)    = convert_index<StorageIndex>(ii);
    u(ii)     = 0;
    jr(ii)    = convert_index<StorageIndex>(ii);
    RealScalar rownorm = 0;

    typename FactorType::InnerIterator j_it(mat, ii); 
    for (; j_it; ++j_it)
    {
      Index k = j_it.index();
      if (k < ii)
      {
        
        ju(sizel) = convert_index<StorageIndex>(k);
        u(sizel) = j_it.value();
        jr(k) = convert_index<StorageIndex>(sizel);
        ++sizel;
      }
      else if (k == ii)
      {
        u(ii) = j_it.value();
      }
      else
      {
        
        Index jpos = ii + sizeu;
        ju(jpos) = convert_index<StorageIndex>(k);
        u(jpos) = j_it.value();
        jr(k) = convert_index<StorageIndex>(jpos);
        ++sizeu;
      }
      rownorm += numext::abs2(j_it.value());
    }

    if(rownorm==0)
    {
      m_info = NumericalIssue;
      return;
    }
    
    rownorm = sqrt(rownorm);

    Index jj = 0;
    Index len = 0;
    while (jj < sizel)
    {

      Index k;
      Index minrow = ju.segment(jj,sizel-jj).minCoeff(&k); 
      k += jj;
      if (minrow != ju(jj))
      {
        
        Index j = ju(jj);
        swap(ju(jj), ju(k));
        jr(minrow) = convert_index<StorageIndex>(jj);
        jr(j) = convert_index<StorageIndex>(k);
        swap(u(jj), u(k));
      }
      
      jr(minrow) = -1;

      typename FactorType::InnerIterator ki_it(m_lu, minrow);
      while (ki_it && ki_it.index() < minrow) ++ki_it;
      eigen_internal_assert(ki_it && ki_it.col()==minrow);
      Scalar fact = u(jj) / ki_it.value();

      if(abs(fact) <= m_droptol)
      {
        jj++;
        continue;
      }

      ++ki_it;
      for (; ki_it; ++ki_it)
      {
        Scalar prod = fact * ki_it.value();
        Index j     = ki_it.index();
        Index jpos  = jr(j);
        if (jpos == -1) 
        {
          Index newpos;
          if (j >= ii) 
          {
            newpos = ii + sizeu;
            sizeu++;
            eigen_internal_assert(sizeu<=n);
          }
          else 
          {
            newpos = sizel;
            sizel++;
            eigen_internal_assert(sizel<=ii);
          }
          ju(newpos) = convert_index<StorageIndex>(j);
          u(newpos) = -prod;
          jr(j) = convert_index<StorageIndex>(newpos);
        }
        else
          u(jpos) -= prod;
      }
      
      u(len)  = fact;
      ju(len) = convert_index<StorageIndex>(minrow);
      ++len;

      jj++;
    } 

    for(Index k = 0; k <sizeu; k++) jr(ju(ii+k)) = -1;

    sizel = len;
    len = (std::min)(sizel, nnzL);
    typename Vector::SegmentReturnType ul(u.segment(0, sizel));
    typename VectorI::SegmentReturnType jul(ju.segment(0, sizel));
    internal::QuickSplit(ul, jul, len);

    m_lu.startVec(ii);
    for(Index k = 0; k < len; k++)
      m_lu.insertBackByOuterInnerUnordered(ii,ju(k)) = u(k);

    if (u(ii) == Scalar(0))
      u(ii) = sqrt(m_droptol) * rownorm;
    m_lu.insertBackByOuterInnerUnordered(ii, ii) = u(ii);

    len = 0;
    for(Index k = 1; k < sizeu; k++)
    {
      if(abs(u(ii+k)) > m_droptol * rownorm )
      {
        ++len;
        u(ii + len)  = u(ii + k);
        ju(ii + len) = ju(ii + k);
      }
    }
    sizeu = len + 1; 
    len = (std::min)(sizeu, nnzU);
    typename Vector::SegmentReturnType uu(u.segment(ii+1, sizeu-1));
    typename VectorI::SegmentReturnType juu(ju.segment(ii+1, sizeu-1));
    internal::QuickSplit(uu, juu, len);

    for(Index k = ii + 1; k < ii + len; k++)
      m_lu.insertBackByOuterInnerUnordered(ii,ju(k)) = u(k);
  }
  m_lu.finalize();
  m_lu.makeCompressed();

  m_factorizationIsOk = true;
  m_info = Success;
}

} 

#endif 
