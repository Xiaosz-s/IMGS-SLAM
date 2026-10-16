#ifndef EIGEN_SPARSEMATRIX_H
#define EIGEN_SPARSEMATRIX_H

namespace Eigen { 

namespace internal {
template<typename Scalar_, int Options_, typename StorageIndex_>
struct traits<SparseMatrix<Scalar_, Options_, StorageIndex_> >
{
  typedef Scalar_ Scalar;
  typedef StorageIndex_ StorageIndex;
  typedef Sparse StorageKind;
  typedef MatrixXpr XprKind;
  enum {
    RowsAtCompileTime = Dynamic,
    ColsAtCompileTime = Dynamic,
    MaxRowsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
    Flags = Options_ | NestByRefBit | LvalueBit | CompressedAccessBit,
    SupportedAccessPatterns = InnerRandomAccessPattern
  };
};

template<typename Scalar_, int Options_, typename StorageIndex_, int DiagIndex>
struct traits<Diagonal<SparseMatrix<Scalar_, Options_, StorageIndex_>, DiagIndex> >
{
  typedef SparseMatrix<Scalar_, Options_, StorageIndex_> MatrixType;
  typedef typename ref_selector<MatrixType>::type MatrixTypeNested;
  typedef typename remove_reference<MatrixTypeNested>::type _MatrixTypeNested;

  typedef Scalar_ Scalar;
  typedef Dense StorageKind;
  typedef StorageIndex_ StorageIndex;
  typedef MatrixXpr XprKind;

  enum {
    RowsAtCompileTime = Dynamic,
    ColsAtCompileTime = 1,
    MaxRowsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = 1,
    Flags = LvalueBit
  };
};

template<typename Scalar_, int Options_, typename StorageIndex_, int DiagIndex>
struct traits<Diagonal<const SparseMatrix<Scalar_, Options_, StorageIndex_>, DiagIndex> >
 : public traits<Diagonal<SparseMatrix<Scalar_, Options_, StorageIndex_>, DiagIndex> >
{
  enum {
    Flags = 0
  };
};

} 

template<typename Scalar_, int Options_, typename StorageIndex_>
class SparseMatrix
  : public SparseCompressedBase<SparseMatrix<Scalar_, Options_, StorageIndex_> >
{
    typedef SparseCompressedBase<SparseMatrix> Base;
    using Base::convert_index;
    friend class SparseVector<Scalar_,0,StorageIndex_>;
    template<typename, typename, typename, typename, typename>
    friend struct internal::Assignment;
  public:
    using Base::isCompressed;
    using Base::nonZeros;
    EIGEN_SPARSE_PUBLIC_INTERFACE(SparseMatrix)
    using Base::operator+=;
    using Base::operator-=;

    typedef MappedSparseMatrix<Scalar,Flags> Map;
    typedef Diagonal<SparseMatrix> DiagonalReturnType;
    typedef Diagonal<const SparseMatrix> ConstDiagonalReturnType;
    typedef typename Base::InnerIterator InnerIterator;
    typedef typename Base::ReverseInnerIterator ReverseInnerIterator;

    using Base::IsRowMajor;
    typedef internal::CompressedStorage<Scalar,StorageIndex> Storage;
    enum {
      Options = Options_
    };

    typedef typename Base::IndexVector IndexVector;
    typedef typename Base::ScalarVector ScalarVector;
  protected:
    typedef SparseMatrix<Scalar,(Flags&~RowMajorBit)|(IsRowMajor?RowMajorBit:0)> TransposedSparseMatrix;

    Index m_outerSize;
    Index m_innerSize;
    StorageIndex* m_outerIndex;
    StorageIndex* m_innerNonZeros;     
    Storage m_data;

  public:

    inline Index rows() const { return IsRowMajor ? m_outerSize : m_innerSize; }
    
    inline Index cols() const { return IsRowMajor ? m_innerSize : m_outerSize; }

    inline Index innerSize() const { return m_innerSize; }
    
    inline Index outerSize() const { return m_outerSize; }

    inline const Scalar* valuePtr() const { return m_data.valuePtr(); }

    inline Scalar* valuePtr() { return m_data.valuePtr(); }

    inline const StorageIndex* innerIndexPtr() const { return m_data.indexPtr(); }

    inline StorageIndex* innerIndexPtr() { return m_data.indexPtr(); }

    inline const StorageIndex* outerIndexPtr() const { return m_outerIndex; }

    inline StorageIndex* outerIndexPtr() { return m_outerIndex; }

    inline const StorageIndex* innerNonZeroPtr() const { return m_innerNonZeros; }

    inline StorageIndex* innerNonZeroPtr() { return m_innerNonZeros; }

    inline Storage& data() { return m_data; }
    
    inline const Storage& data() const { return m_data; }

    inline Scalar coeff(Index row, Index col) const
    {
      eigen_assert(row>=0 && row<rows() && col>=0 && col<cols());
      
      const Index outer = IsRowMajor ? row : col;
      const Index inner = IsRowMajor ? col : row;
      Index end = m_innerNonZeros ? m_outerIndex[outer] + m_innerNonZeros[outer] : m_outerIndex[outer+1];
      return m_data.atInRange(m_outerIndex[outer], end, StorageIndex(inner));
    }

    inline Scalar& coeffRef(Index row, Index col)
    {
      eigen_assert(row>=0 && row<rows() && col>=0 && col<cols());
      
      const Index outer = IsRowMajor ? row : col;
      const Index inner = IsRowMajor ? col : row;

      Index start = m_outerIndex[outer];
      Index end = m_innerNonZeros ? m_outerIndex[outer] + m_innerNonZeros[outer] : m_outerIndex[outer+1];
      eigen_assert(end>=start && "you probably called coeffRef on a non finalized matrix");
      if(end<=start)
        return insert(row,col);
      const Index p = m_data.searchLowerIndex(start,end-1,StorageIndex(inner));
      if((p<end) && (m_data.index(p)==inner))
        return m_data.value(p);
      else
        return insert(row,col);
    }

    Scalar& insert(Index row, Index col);

  public:

    inline void setZero()
    {
      m_data.clear();
      std::fill_n(m_outerIndex, m_outerSize + 1, StorageIndex(0));
      if(m_innerNonZeros) {
        std::fill_n(m_innerNonZeros, m_outerSize, StorageIndex(0));
      }
    }

    inline void reserve(Index reserveSize)
    {
      eigen_assert(isCompressed() && "This function does not make sense in non compressed mode.");
      m_data.reserve(reserveSize);
    }
    
    #ifdef EIGEN_PARSED_BY_DOXYGEN

    template<class SizesType>
    inline void reserve(const SizesType& reserveSizes);
    #else
    template<class SizesType>
    inline void reserve(const SizesType& reserveSizes, const typename SizesType::value_type& enableif =
    #if (!EIGEN_COMP_MSVC) || (EIGEN_COMP_MSVC>=1500) 
        typename
    #endif
        SizesType::value_type())
    {
      EIGEN_UNUSED_VARIABLE(enableif);
      reserveInnerVectors(reserveSizes);
    }
    #endif 
  protected:
    template<class SizesType>
    inline void reserveInnerVectors(const SizesType& reserveSizes)
    {
      if(isCompressed())
      {
        Index totalReserveSize = 0;
        
        m_innerNonZeros = static_cast<StorageIndex*>(std::malloc(m_outerSize * sizeof(StorageIndex)));
        if (!m_innerNonZeros) internal::throw_std_bad_alloc();

        StorageIndex* newOuterIndex = m_innerNonZeros;
        
        StorageIndex count = 0;
        for(Index j=0; j<m_outerSize; ++j)
        {
          newOuterIndex[j] = count;
          count += reserveSizes[j] + (m_outerIndex[j+1]-m_outerIndex[j]);
          totalReserveSize += reserveSizes[j];
        }
        m_data.reserve(totalReserveSize);
        StorageIndex previousOuterIndex = m_outerIndex[m_outerSize];
        for(Index j=m_outerSize-1; j>=0; --j)
        {
          StorageIndex innerNNZ = previousOuterIndex - m_outerIndex[j];
          for(Index i=innerNNZ-1; i>=0; --i)
          {
            m_data.index(newOuterIndex[j]+i) = m_data.index(m_outerIndex[j]+i);
            m_data.value(newOuterIndex[j]+i) = m_data.value(m_outerIndex[j]+i);
          }
          previousOuterIndex = m_outerIndex[j];
          m_outerIndex[j] = newOuterIndex[j];
          m_innerNonZeros[j] = innerNNZ;
        }
        if(m_outerSize>0)
          m_outerIndex[m_outerSize] = m_outerIndex[m_outerSize-1] + m_innerNonZeros[m_outerSize-1] + reserveSizes[m_outerSize-1];
        
        m_data.resize(m_outerIndex[m_outerSize]);
      }
      else
      {
        StorageIndex* newOuterIndex = static_cast<StorageIndex*>(std::malloc((m_outerSize+1)*sizeof(StorageIndex)));
        if (!newOuterIndex) internal::throw_std_bad_alloc();
        
        StorageIndex count = 0;
        for(Index j=0; j<m_outerSize; ++j)
        {
          newOuterIndex[j] = count;
          StorageIndex alreadyReserved = (m_outerIndex[j+1]-m_outerIndex[j]) - m_innerNonZeros[j];
          StorageIndex toReserve = std::max<StorageIndex>(reserveSizes[j], alreadyReserved);
          count += toReserve + m_innerNonZeros[j];
        }
        newOuterIndex[m_outerSize] = count;
        
        m_data.resize(count);
        for(Index j=m_outerSize-1; j>=0; --j)
        {
          Index offset = newOuterIndex[j] - m_outerIndex[j];
          if(offset>0)
          {
            StorageIndex innerNNZ = m_innerNonZeros[j];
            for(Index i=innerNNZ-1; i>=0; --i)
            {
              m_data.index(newOuterIndex[j]+i) = m_data.index(m_outerIndex[j]+i);
              m_data.value(newOuterIndex[j]+i) = m_data.value(m_outerIndex[j]+i);
            }
          }
        }
        
        std::swap(m_outerIndex, newOuterIndex);
        std::free(newOuterIndex);
      }
      
    }
  public:

    inline Scalar& insertBack(Index row, Index col)
    {
      return insertBackByOuterInner(IsRowMajor?row:col, IsRowMajor?col:row);
    }

    inline Scalar& insertBackByOuterInner(Index outer, Index inner)
    {
      eigen_assert(Index(m_outerIndex[outer+1]) == m_data.size() && "Invalid ordered insertion (invalid outer index)");
      eigen_assert( (m_outerIndex[outer+1]-m_outerIndex[outer]==0 || m_data.index(m_data.size()-1)<inner) && "Invalid ordered insertion (invalid inner index)");
      Index p = m_outerIndex[outer+1];
      ++m_outerIndex[outer+1];
      m_data.append(Scalar(0), inner);
      return m_data.value(p);
    }

    inline Scalar& insertBackByOuterInnerUnordered(Index outer, Index inner)
    {
      Index p = m_outerIndex[outer+1];
      ++m_outerIndex[outer+1];
      m_data.append(Scalar(0), inner);
      return m_data.value(p);
    }

    inline void startVec(Index outer)
    {
      eigen_assert(m_outerIndex[outer]==Index(m_data.size()) && "You must call startVec for each inner vector sequentially");
      eigen_assert(m_outerIndex[outer+1]==0 && "You must call startVec for each inner vector sequentially");
      m_outerIndex[outer+1] = m_outerIndex[outer];
    }

    inline void finalize()
    {
      if(isCompressed())
      {
        StorageIndex size = internal::convert_index<StorageIndex>(m_data.size());
        Index i = m_outerSize;
        
        while (i>=0 && m_outerIndex[i]==0)
          --i;
        ++i;
        while (i<=m_outerSize)
        {
          m_outerIndex[i] = size;
          ++i;
        }
      }
    }

    template<typename InputIterators>
    void setFromTriplets(const InputIterators& begin, const InputIterators& end);

    template<typename InputIterators,typename DupFunctor>
    void setFromTriplets(const InputIterators& begin, const InputIterators& end, DupFunctor dup_func);

    void sumupDuplicates() { collapseDuplicates(internal::scalar_sum_op<Scalar,Scalar>()); }

    template<typename DupFunctor>
    void collapseDuplicates(DupFunctor dup_func = DupFunctor());

    Scalar& insertByOuterInner(Index j, Index i)
    {
      return insert(IsRowMajor ? j : i, IsRowMajor ? i : j);
    }

    void makeCompressed()
    {
      if(isCompressed())
        return;
      
      eigen_internal_assert(m_outerIndex!=0 && m_outerSize>0);
      
      Index oldStart = m_outerIndex[1];
      m_outerIndex[1] = m_innerNonZeros[0];
      for(Index j=1; j<m_outerSize; ++j)
      {
        Index nextOldStart = m_outerIndex[j+1];
        Index offset = oldStart - m_outerIndex[j];
        if(offset>0)
        {
          for(Index k=0; k<m_innerNonZeros[j]; ++k)
          {
            m_data.index(m_outerIndex[j]+k) = m_data.index(oldStart+k);
            m_data.value(m_outerIndex[j]+k) = m_data.value(oldStart+k);
          }
        }
        m_outerIndex[j+1] = m_outerIndex[j] + m_innerNonZeros[j];
        oldStart = nextOldStart;
      }
      std::free(m_innerNonZeros);
      m_innerNonZeros = 0;
      m_data.resize(m_outerIndex[m_outerSize]);
      m_data.squeeze();
    }

    void uncompress()
    {
      if(m_innerNonZeros != 0)
        return; 
      m_innerNonZeros = static_cast<StorageIndex*>(std::malloc(m_outerSize * sizeof(StorageIndex)));
      for (Index i = 0; i < m_outerSize; i++)
      {
        m_innerNonZeros[i] = m_outerIndex[i+1] - m_outerIndex[i]; 
      }
    }

    void prune(const Scalar& reference, const RealScalar& epsilon = NumTraits<RealScalar>::dummy_precision())
    {
      prune(default_prunning_func(reference,epsilon));
    }

    template<typename KeepFunc>
    void prune(const KeepFunc& keep = KeepFunc())
    {
      
      makeCompressed();

      StorageIndex k = 0;
      for(Index j=0; j<m_outerSize; ++j)
      {
        Index previousStart = m_outerIndex[j];
        m_outerIndex[j] = k;
        Index end = m_outerIndex[j+1];
        for(Index i=previousStart; i<end; ++i)
        {
          if(keep(IsRowMajor?j:m_data.index(i), IsRowMajor?m_data.index(i):j, m_data.value(i)))
          {
            m_data.value(k) = m_data.value(i);
            m_data.index(k) = m_data.index(i);
            ++k;
          }
        }
      }
      m_outerIndex[m_outerSize] = k;
      m_data.resize(k,0);
    }

    void conservativeResize(Index rows, Index cols) 
    {
      
      if (this->rows() == rows && this->cols() == cols) return;

      if(rows==0 || cols==0) return resize(rows,cols);

      Index innerChange = IsRowMajor ? cols - this->cols() : rows - this->rows();
      Index outerChange = IsRowMajor ? rows - this->rows() : cols - this->cols();
      StorageIndex newInnerSize = convert_index(IsRowMajor ? cols : rows);

      if (m_innerNonZeros)
      {
        
        StorageIndex *newInnerNonZeros = static_cast<StorageIndex*>(std::realloc(m_innerNonZeros, (m_outerSize + outerChange) * sizeof(StorageIndex)));
        if (!newInnerNonZeros) internal::throw_std_bad_alloc();
        m_innerNonZeros = newInnerNonZeros;
        
        for(Index i=m_outerSize; i<m_outerSize+outerChange; i++)          
          m_innerNonZeros[i] = 0;
      } 
      else if (innerChange < 0) 
      {
        
        m_innerNonZeros = static_cast<StorageIndex*>(std::malloc((m_outerSize + outerChange) * sizeof(StorageIndex)));
        if (!m_innerNonZeros) internal::throw_std_bad_alloc();
        for(Index i = 0; i < m_outerSize + (std::min)(outerChange, Index(0)); i++)
          m_innerNonZeros[i] = m_outerIndex[i+1] - m_outerIndex[i];
        for(Index i = m_outerSize; i < m_outerSize + outerChange; i++)
          m_innerNonZeros[i] = 0;
      }

      if (m_innerNonZeros && innerChange < 0)
      {
        for(Index i = 0; i < m_outerSize + (std::min)(outerChange, Index(0)); i++)
        {
          StorageIndex &n = m_innerNonZeros[i];
          StorageIndex start = m_outerIndex[i];
          while (n > 0 && m_data.index(start+n-1) >= newInnerSize) --n; 
        }
      }
      
      m_innerSize = newInnerSize;

      if (outerChange == 0)
        return;
          
      StorageIndex *newOuterIndex = static_cast<StorageIndex*>(std::realloc(m_outerIndex, (m_outerSize + outerChange + 1) * sizeof(StorageIndex)));
      if (!newOuterIndex) internal::throw_std_bad_alloc();
      m_outerIndex = newOuterIndex;
      if (outerChange > 0)
      {
        StorageIndex lastIdx = m_outerSize == 0 ? 0 : m_outerIndex[m_outerSize];
        for(Index i=m_outerSize; i<m_outerSize+outerChange+1; i++)          
          m_outerIndex[i] = lastIdx; 
      }
      m_outerSize += outerChange;
    }

    void resize(Index rows, Index cols)
    {
      const Index outerSize = IsRowMajor ? rows : cols;
      m_innerSize = IsRowMajor ? cols : rows;
      m_data.clear();
      if (m_outerSize != outerSize || m_outerSize==0)
      {
        std::free(m_outerIndex);
        m_outerIndex = static_cast<StorageIndex*>(std::malloc((outerSize + 1) * sizeof(StorageIndex)));
        if (!m_outerIndex) internal::throw_std_bad_alloc();
        
        m_outerSize = outerSize;
      }
      if(m_innerNonZeros)
      {
        std::free(m_innerNonZeros);
        m_innerNonZeros = 0;
      }
      std::fill_n(m_outerIndex, m_outerSize + 1, StorageIndex(0));
    }

    void resizeNonZeros(Index size)
    {
      m_data.resize(size);
    }

    const ConstDiagonalReturnType diagonal() const { return ConstDiagonalReturnType(*this); }

    DiagonalReturnType diagonal() { return DiagonalReturnType(*this); }

    inline SparseMatrix()
      : m_outerSize(-1), m_innerSize(0), m_outerIndex(0), m_innerNonZeros(0)
    {
      check_template_parameters();
      resize(0, 0);
    }

    inline SparseMatrix(Index rows, Index cols)
      : m_outerSize(0), m_innerSize(0), m_outerIndex(0), m_innerNonZeros(0)
    {
      check_template_parameters();
      resize(rows, cols);
    }

    template<typename OtherDerived>
    inline SparseMatrix(const SparseMatrixBase<OtherDerived>& other)
      : m_outerSize(0), m_innerSize(0), m_outerIndex(0), m_innerNonZeros(0)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<Scalar, typename OtherDerived::Scalar>::value),
        YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)
      check_template_parameters();
      const bool needToTranspose = (Flags & RowMajorBit) != (internal::evaluator<OtherDerived>::Flags & RowMajorBit);
      if (needToTranspose)
        *this = other.derived();
      else
      {
        #ifdef EIGEN_SPARSE_CREATE_TEMPORARY_PLUGIN
          EIGEN_SPARSE_CREATE_TEMPORARY_PLUGIN
        #endif
        internal::call_assignment_no_alias(*this, other.derived());
      }
    }

    template<typename OtherDerived, unsigned int UpLo>
    inline SparseMatrix(const SparseSelfAdjointView<OtherDerived, UpLo>& other)
      : m_outerSize(0), m_innerSize(0), m_outerIndex(0), m_innerNonZeros(0)
    {
      check_template_parameters();
      Base::operator=(other);
    }

    inline SparseMatrix(const SparseMatrix& other)
      : Base(), m_outerSize(0), m_innerSize(0), m_outerIndex(0), m_innerNonZeros(0)
    {
      check_template_parameters();
      *this = other.derived();
    }

    template<typename OtherDerived>
    SparseMatrix(const ReturnByValue<OtherDerived>& other)
      : Base(), m_outerSize(0), m_innerSize(0), m_outerIndex(0), m_innerNonZeros(0)
    {
      check_template_parameters();
      initAssignment(other);
      other.evalTo(*this);
    }

    template<typename OtherDerived>
    explicit SparseMatrix(const DiagonalBase<OtherDerived>& other)
      : Base(), m_outerSize(0), m_innerSize(0), m_outerIndex(0), m_innerNonZeros(0)
    {
      check_template_parameters();
      *this = other.derived();
    }

    inline void swap(SparseMatrix& other)
    {
      
      std::swap(m_outerIndex, other.m_outerIndex);
      std::swap(m_innerSize, other.m_innerSize);
      std::swap(m_outerSize, other.m_outerSize);
      std::swap(m_innerNonZeros, other.m_innerNonZeros);
      m_data.swap(other.m_data);
    }

    inline void setIdentity()
    {
      eigen_assert(rows() == cols() && "ONLY FOR SQUARED MATRICES");
      this->m_data.resize(rows());
      Eigen::Map<IndexVector>(this->m_data.indexPtr(), rows()).setLinSpaced(0, StorageIndex(rows()-1));
      Eigen::Map<ScalarVector>(this->m_data.valuePtr(), rows()).setOnes();
      Eigen::Map<IndexVector>(this->m_outerIndex, rows()+1).setLinSpaced(0, StorageIndex(rows()));
      std::free(m_innerNonZeros);
      m_innerNonZeros = 0;
    }
    inline SparseMatrix& operator=(const SparseMatrix& other)
    {
      if (other.isRValue())
      {
        swap(other.const_cast_derived());
      }
      else if(this!=&other)
      {
        #ifdef EIGEN_SPARSE_CREATE_TEMPORARY_PLUGIN
          EIGEN_SPARSE_CREATE_TEMPORARY_PLUGIN
        #endif
        initAssignment(other);
        if(other.isCompressed())
        {
          internal::smart_copy(other.m_outerIndex, other.m_outerIndex + m_outerSize + 1, m_outerIndex);
          m_data = other.m_data;
        }
        else
        {
          Base::operator=(other);
        }
      }
      return *this;
    }

#ifndef EIGEN_PARSED_BY_DOXYGEN
    template<typename OtherDerived>
    inline SparseMatrix& operator=(const EigenBase<OtherDerived>& other)
    { return Base::operator=(other.derived()); }

    template<typename Lhs, typename Rhs>
    inline SparseMatrix& operator=(const Product<Lhs,Rhs,AliasFreeProduct>& other);
#endif 

    template<typename OtherDerived>
    EIGEN_DONT_INLINE SparseMatrix& operator=(const SparseMatrixBase<OtherDerived>& other);

    friend std::ostream & operator << (std::ostream & s, const SparseMatrix& m)
    {
      EIGEN_DBG_SPARSE(
        s << "Nonzero entries:\n";
        if(m.isCompressed())
        {
          for (Index i=0; i<m.nonZeros(); ++i)
            s << "(" << m.m_data.value(i) << "," << m.m_data.index(i) << ") ";
        }
        else
        {
          for (Index i=0; i<m.outerSize(); ++i)
          {
            Index p = m.m_outerIndex[i];
            Index pe = m.m_outerIndex[i]+m.m_innerNonZeros[i];
            Index k=p;
            for (; k<pe; ++k) {
              s << "(" << m.m_data.value(k) << "," << m.m_data.index(k) << ") ";
            }
            for (; k<m.m_outerIndex[i+1]; ++k) {
              s << "(_,_) ";
            }
          }
        }
        s << std::endl;
        s << std::endl;
        s << "Outer pointers:\n";
        for (Index i=0; i<m.outerSize(); ++i) {
          s << m.m_outerIndex[i] << " ";
        }
        s << " $" << std::endl;
        if(!m.isCompressed())
        {
          s << "Inner non zeros:\n";
          for (Index i=0; i<m.outerSize(); ++i) {
            s << m.m_innerNonZeros[i] << " ";
          }
          s << " $" << std::endl;
        }
        s << std::endl;
      );
      s << static_cast<const SparseMatrixBase<SparseMatrix>&>(m);
      return s;
    }

    inline ~SparseMatrix()
    {
      std::free(m_outerIndex);
      std::free(m_innerNonZeros);
    }

    Scalar sum() const;
    
#   ifdef EIGEN_SPARSEMATRIX_PLUGIN
#     include EIGEN_SPARSEMATRIX_PLUGIN
#   endif

protected:

    template<typename Other>
    void initAssignment(const Other& other)
    {
      resize(other.rows(), other.cols());
      if(m_innerNonZeros)
      {
        std::free(m_innerNonZeros);
        m_innerNonZeros = 0;
      }
    }

    EIGEN_DONT_INLINE Scalar& insertCompressed(Index row, Index col);

    class SingletonVector
    {
        StorageIndex m_index;
        StorageIndex m_value;
      public:
        typedef StorageIndex value_type;
        SingletonVector(Index i, Index v)
          : m_index(convert_index(i)), m_value(convert_index(v))
        {}

        StorageIndex operator[](Index i) const { return i==m_index ? m_value : 0; }
    };

    EIGEN_DONT_INLINE Scalar& insertUncompressed(Index row, Index col);

public:

    EIGEN_STRONG_INLINE Scalar& insertBackUncompressed(Index row, Index col)
    {
      const Index outer = IsRowMajor ? row : col;
      const Index inner = IsRowMajor ? col : row;

      eigen_assert(!isCompressed());
      eigen_assert(m_innerNonZeros[outer]<=(m_outerIndex[outer+1] - m_outerIndex[outer]));

      Index p = m_outerIndex[outer] + m_innerNonZeros[outer]++;
      m_data.index(p) = convert_index(inner);
      return (m_data.value(p) = Scalar(0));
    }
protected:
    struct IndexPosPair {
      IndexPosPair(Index a_i, Index a_p) : i(a_i), p(a_p) {}
      Index i;
      Index p;
    };

    template<typename DiagXpr, typename Func>
    void assignDiagonal(const DiagXpr diagXpr, const Func& assignFunc)
    {
      Index n = diagXpr.size();

      const bool overwrite = internal::is_same<Func, internal::assign_op<Scalar,Scalar> >::value;
      if(overwrite)
      {
        if((this->rows()!=n) || (this->cols()!=n))
          this->resize(n, n);
      }

      if(m_data.size()==0 || overwrite)
      {
        typedef Array<StorageIndex,Dynamic,1> ArrayXI;  
        this->makeCompressed();
        this->resizeNonZeros(n);
        Eigen::Map<ArrayXI>(this->innerIndexPtr(), n).setLinSpaced(0,StorageIndex(n)-1);
        Eigen::Map<ArrayXI>(this->outerIndexPtr(), n+1).setLinSpaced(0,StorageIndex(n));
        Eigen::Map<Array<Scalar,Dynamic,1> > values = this->coeffs();
        values.setZero();
        internal::call_assignment_no_alias(values, diagXpr, assignFunc);
      }
      else
      {
        bool isComp = isCompressed();
        internal::evaluator<DiagXpr> diaEval(diagXpr);
        std::vector<IndexPosPair> newEntries;

        for(Index i = 0; i<n; ++i)
        {
          internal::LowerBoundIndex lb = this->lower_bound(i,i);
          Index p = lb.value;
          if(lb.found)
          {
            
            assignFunc.assignCoeff(m_data.value(p), diaEval.coeff(i));
          }
          else if((!isComp) && m_innerNonZeros[i] < (m_outerIndex[i+1]-m_outerIndex[i]))
          {
            
            m_data.moveChunk(p, p+1, m_outerIndex[i]+m_innerNonZeros[i]-p);
            m_innerNonZeros[i]++;
            m_data.value(p) = Scalar(0);
            m_data.index(p) = StorageIndex(i);
            assignFunc.assignCoeff(m_data.value(p), diaEval.coeff(i));
          }
          else
          {
            
            newEntries.push_back(IndexPosPair(i,p));
          }
        }
        
        Index n_entries = Index(newEntries.size());
        if(n_entries>0)
        {
          Storage newData(m_data.size()+n_entries);
          Index prev_p = 0;
          Index prev_i = 0;
          for(Index k=0; k<n_entries;++k)
          {
            Index i = newEntries[k].i;
            Index p = newEntries[k].p;
            internal::smart_copy(m_data.valuePtr()+prev_p, m_data.valuePtr()+p, newData.valuePtr()+prev_p+k);
            internal::smart_copy(m_data.indexPtr()+prev_p, m_data.indexPtr()+p, newData.indexPtr()+prev_p+k);
            for(Index j=prev_i;j<i;++j)
              m_outerIndex[j+1] += k;
            if(!isComp)
              m_innerNonZeros[i]++;
            prev_p = p;
            prev_i = i;
            newData.value(p+k) = Scalar(0);
            newData.index(p+k) = StorageIndex(i);
            assignFunc.assignCoeff(newData.value(p+k), diaEval.coeff(i));
          }
          {
            internal::smart_copy(m_data.valuePtr()+prev_p, m_data.valuePtr()+m_data.size(), newData.valuePtr()+prev_p+n_entries);
            internal::smart_copy(m_data.indexPtr()+prev_p, m_data.indexPtr()+m_data.size(), newData.indexPtr()+prev_p+n_entries);
            for(Index j=prev_i+1;j<=m_outerSize;++j)
              m_outerIndex[j] += n_entries;
          }
          m_data.swap(newData);
        }
      }
    }

private:
  static void check_template_parameters()
  {
    EIGEN_STATIC_ASSERT(NumTraits<StorageIndex>::IsSigned,THE_INDEX_TYPE_MUST_BE_A_SIGNED_TYPE);
    EIGEN_STATIC_ASSERT((Options&(ColMajor|RowMajor))==Options,INVALID_MATRIX_TEMPLATE_PARAMETERS);
  }

  struct default_prunning_func {
    default_prunning_func(const Scalar& ref, const RealScalar& eps) : reference(ref), epsilon(eps) {}
    inline bool operator() (const Index&, const Index&, const Scalar& value) const
    {
      return !internal::isMuchSmallerThan(value, reference, epsilon);
    }
    Scalar reference;
    RealScalar epsilon;
  };
};

namespace internal {

template<typename InputIterator, typename SparseMatrixType, typename DupFunctor>
void set_from_triplets(const InputIterator& begin, const InputIterator& end, SparseMatrixType& mat, DupFunctor dup_func)
{
  enum { IsRowMajor = SparseMatrixType::IsRowMajor };
  typedef typename SparseMatrixType::Scalar Scalar;
  typedef typename SparseMatrixType::StorageIndex StorageIndex;
  SparseMatrix<Scalar,IsRowMajor?ColMajor:RowMajor,StorageIndex> trMat(mat.rows(),mat.cols());

  if(begin!=end)
  {
    
    typename SparseMatrixType::IndexVector wi(trMat.outerSize());
    wi.setZero();
    for(InputIterator it(begin); it!=end; ++it)
    {
      eigen_assert(it->row()>=0 && it->row()<mat.rows() && it->col()>=0 && it->col()<mat.cols());
      wi(IsRowMajor ? it->col() : it->row())++;
    }

    trMat.reserve(wi);
    for(InputIterator it(begin); it!=end; ++it)
      trMat.insertBackUncompressed(it->row(),it->col()) = it->value();

    trMat.collapseDuplicates(dup_func);
  }

  mat = trMat;
}

}

template<typename Scalar, int Options_, typename StorageIndex_>
template<typename InputIterators>
void SparseMatrix<Scalar,Options_,StorageIndex_>::setFromTriplets(const InputIterators& begin, const InputIterators& end)
{
  internal::set_from_triplets<InputIterators, SparseMatrix<Scalar,Options_,StorageIndex_> >(begin, end, *this, internal::scalar_sum_op<Scalar,Scalar>());
}

template<typename Scalar, int Options_, typename StorageIndex_>
template<typename InputIterators,typename DupFunctor>
void SparseMatrix<Scalar,Options_,StorageIndex_>::setFromTriplets(const InputIterators& begin, const InputIterators& end, DupFunctor dup_func)
{
  internal::set_from_triplets<InputIterators, SparseMatrix<Scalar,Options_,StorageIndex_>, DupFunctor>(begin, end, *this, dup_func);
}

template<typename Scalar, int Options_, typename StorageIndex_>
template<typename DupFunctor>
void SparseMatrix<Scalar,Options_,StorageIndex_>::collapseDuplicates(DupFunctor dup_func)
{
  eigen_assert(!isCompressed());
  
  IndexVector wi(innerSize());
  wi.fill(-1);
  StorageIndex count = 0;
  
  for(Index j=0; j<outerSize(); ++j)
  {
    StorageIndex start   = count;
    Index oldEnd  = m_outerIndex[j]+m_innerNonZeros[j];
    for(Index k=m_outerIndex[j]; k<oldEnd; ++k)
    {
      Index i = m_data.index(k);
      if(wi(i)>=start)
      {
        
        m_data.value(wi(i)) = dup_func(m_data.value(wi(i)), m_data.value(k));
      }
      else
      {
        m_data.value(count) = m_data.value(k);
        m_data.index(count) = m_data.index(k);
        wi(i) = count;
        ++count;
      }
    }
    m_outerIndex[j] = start;
  }
  m_outerIndex[m_outerSize] = count;

  std::free(m_innerNonZeros);
  m_innerNonZeros = 0;
  m_data.resize(m_outerIndex[m_outerSize]);
}

template<typename Scalar, int Options_, typename StorageIndex_>
template<typename OtherDerived>
EIGEN_DONT_INLINE SparseMatrix<Scalar,Options_,StorageIndex_>& SparseMatrix<Scalar,Options_,StorageIndex_>::operator=(const SparseMatrixBase<OtherDerived>& other)
{
  EIGEN_STATIC_ASSERT((internal::is_same<Scalar, typename OtherDerived::Scalar>::value),
        YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)

  #ifdef EIGEN_SPARSE_CREATE_TEMPORARY_PLUGIN
    EIGEN_SPARSE_CREATE_TEMPORARY_PLUGIN
  #endif
      
  const bool needToTranspose = (Flags & RowMajorBit) != (internal::evaluator<OtherDerived>::Flags & RowMajorBit);
  if (needToTranspose)
  {
    #ifdef EIGEN_SPARSE_TRANSPOSED_COPY_PLUGIN
      EIGEN_SPARSE_TRANSPOSED_COPY_PLUGIN
    #endif

    typedef typename internal::nested_eval<OtherDerived,2,typename internal::plain_matrix_type<OtherDerived>::type >::type OtherCopy;
    typedef typename internal::remove_all<OtherCopy>::type _OtherCopy;
    typedef internal::evaluator<_OtherCopy> OtherCopyEval;
    OtherCopy otherCopy(other.derived());
    OtherCopyEval otherCopyEval(otherCopy);

    SparseMatrix dest(other.rows(),other.cols());
    Eigen::Map<IndexVector> (dest.m_outerIndex,dest.outerSize()).setZero();

    for (Index j=0; j<otherCopy.outerSize(); ++j)
      for (typename OtherCopyEval::InnerIterator it(otherCopyEval, j); it; ++it)
        ++dest.m_outerIndex[it.index()];

    StorageIndex count = 0;
    IndexVector positions(dest.outerSize());
    for (Index j=0; j<dest.outerSize(); ++j)
    {
      StorageIndex tmp = dest.m_outerIndex[j];
      dest.m_outerIndex[j] = count;
      positions[j] = count;
      count += tmp;
    }
    dest.m_outerIndex[dest.outerSize()] = count;
    
    dest.m_data.resize(count);
    
    for (StorageIndex j=0; j<otherCopy.outerSize(); ++j)
    {
      for (typename OtherCopyEval::InnerIterator it(otherCopyEval, j); it; ++it)
      {
        Index pos = positions[it.index()]++;
        dest.m_data.index(pos) = j;
        dest.m_data.value(pos) = it.value();
      }
    }
    this->swap(dest);
    return *this;
  }
  else
  {
    if(other.isRValue())
    {
      initAssignment(other.derived());
    }
    
    return Base::operator=(other.derived());
  }
}

template<typename Scalar_, int Options_, typename StorageIndex_>
typename SparseMatrix<Scalar_,Options_,StorageIndex_>::Scalar& SparseMatrix<Scalar_,Options_,StorageIndex_>::insert(Index row, Index col)
{
  eigen_assert(row>=0 && row<rows() && col>=0 && col<cols());
  
  const Index outer = IsRowMajor ? row : col;
  const Index inner = IsRowMajor ? col : row;
  
  if(isCompressed())
  {
    if(nonZeros()==0)
    {
      
      if(m_data.allocatedSize()==0)
        m_data.reserve(2*m_innerSize);

      m_innerNonZeros = static_cast<StorageIndex*>(std::malloc(m_outerSize * sizeof(StorageIndex)));
      if(!m_innerNonZeros) internal::throw_std_bad_alloc();
      
      std::fill(m_innerNonZeros, m_innerNonZeros + m_outerSize, StorageIndex(0));

      StorageIndex end = convert_index(m_data.allocatedSize());
      for(Index j=1; j<=m_outerSize; ++j)
        m_outerIndex[j] = end;
    }
    else
    {
      
      m_innerNonZeros = static_cast<StorageIndex*>(std::malloc(m_outerSize * sizeof(StorageIndex)));
      if(!m_innerNonZeros) internal::throw_std_bad_alloc();
      for(Index j=0; j<m_outerSize; ++j)
        m_innerNonZeros[j] = m_outerIndex[j+1]-m_outerIndex[j];
    }
  }

  Index data_end = m_data.allocatedSize();

  if(m_outerIndex[outer]==data_end)
  {
    eigen_internal_assert(m_innerNonZeros[outer]==0);

    StorageIndex p = convert_index(m_data.size());
    Index j = outer;
    while(j>=0 && m_innerNonZeros[j]==0)
      m_outerIndex[j--] = p;

    ++m_innerNonZeros[outer];
    m_data.append(Scalar(0), inner);

    if(data_end != m_data.allocatedSize())
    {

      eigen_internal_assert(data_end < m_data.allocatedSize());
      StorageIndex new_end = convert_index(m_data.allocatedSize());
      for(Index k=outer+1; k<=m_outerSize; ++k)
        if(m_outerIndex[k]==data_end)
          m_outerIndex[k] = new_end;
    }
    return m_data.value(p);
  }

  if(m_outerIndex[outer+1]==data_end && m_outerIndex[outer]+m_innerNonZeros[outer]==m_data.size())
  {
    eigen_internal_assert(outer+1==m_outerSize || m_innerNonZeros[outer+1]==0);

    ++m_innerNonZeros[outer];
    m_data.resize(m_data.size()+1);

    if(data_end != m_data.allocatedSize())
    {

      eigen_internal_assert(data_end < m_data.allocatedSize());
      StorageIndex new_end = convert_index(m_data.allocatedSize());
      for(Index k=outer+1; k<=m_outerSize; ++k)
        if(m_outerIndex[k]==data_end)
          m_outerIndex[k] = new_end;
    }

    Index startId = m_outerIndex[outer];
    Index p = m_outerIndex[outer]+m_innerNonZeros[outer]-1;
    while ( (p > startId) && (m_data.index(p-1) > inner) )
    {
      m_data.index(p) = m_data.index(p-1);
      m_data.value(p) = m_data.value(p-1);
      --p;
    }
    
    m_data.index(p) = convert_index(inner);
    return (m_data.value(p) = Scalar(0));
  }
  
  if(m_data.size() != m_data.allocatedSize())
  {
    
    m_data.resize(m_data.allocatedSize());
    this->reserveInnerVectors(Array<StorageIndex,Dynamic,1>::Constant(m_outerSize, 2));
  }
  
  return insertUncompressed(row,col);
}
    
template<typename Scalar_, int Options_, typename StorageIndex_>
EIGEN_DONT_INLINE typename SparseMatrix<Scalar_,Options_,StorageIndex_>::Scalar& SparseMatrix<Scalar_,Options_,StorageIndex_>::insertUncompressed(Index row, Index col)
{
  eigen_assert(!isCompressed());

  const Index outer = IsRowMajor ? row : col;
  const StorageIndex inner = convert_index(IsRowMajor ? col : row);

  Index room = m_outerIndex[outer+1] - m_outerIndex[outer];
  StorageIndex innerNNZ = m_innerNonZeros[outer];
  if(innerNNZ>=room)
  {
    
    reserve(SingletonVector(outer,std::max<StorageIndex>(2,innerNNZ)));
  }

  Index startId = m_outerIndex[outer];
  Index p = startId + m_innerNonZeros[outer];
  while ( (p > startId) && (m_data.index(p-1) > inner) )
  {
    m_data.index(p) = m_data.index(p-1);
    m_data.value(p) = m_data.value(p-1);
    --p;
  }
  eigen_assert((p<=startId || m_data.index(p-1)!=inner) && "you cannot insert an element that already exists, you must call coeffRef to this end");

  m_innerNonZeros[outer]++;

  m_data.index(p) = inner;
  return (m_data.value(p) = Scalar(0));
}

template<typename Scalar_, int Options_, typename StorageIndex_>
EIGEN_DONT_INLINE typename SparseMatrix<Scalar_,Options_,StorageIndex_>::Scalar& SparseMatrix<Scalar_,Options_,StorageIndex_>::insertCompressed(Index row, Index col)
{
  eigen_assert(isCompressed());

  const Index outer = IsRowMajor ? row : col;
  const Index inner = IsRowMajor ? col : row;

  Index previousOuter = outer;
  if (m_outerIndex[outer+1]==0)
  {
    
    while (previousOuter>=0 && m_outerIndex[previousOuter]==0)
    {
      m_outerIndex[previousOuter] = convert_index(m_data.size());
      --previousOuter;
    }
    m_outerIndex[outer+1] = m_outerIndex[outer];
  }

  bool isLastVec = (!(previousOuter==-1 && m_data.size()!=0))
                && (std::size_t(m_outerIndex[outer+1]) == m_data.size());

  std::size_t startId = m_outerIndex[outer];
  
  std::size_t p = m_outerIndex[outer+1];
  ++m_outerIndex[outer+1];

  double reallocRatio = 1;
  if (m_data.allocatedSize()<=m_data.size())
  {
    
    if (m_data.size()==0)
    {
      m_data.reserve(32);
    }
    else
    {

      double nnzEstimate = double(m_outerIndex[outer])*double(m_outerSize)/double(outer+1);
      reallocRatio = (nnzEstimate-double(m_data.size()))/double(m_data.size());

      reallocRatio = (std::min)((std::max)(reallocRatio,1.5),8.);
    }
  }
  m_data.resize(m_data.size()+1,reallocRatio);

  if (!isLastVec)
  {
    if (previousOuter==-1)
    {

      for (Index k=0; k<=(outer+1); ++k)
        m_outerIndex[k] = 0;
      Index k=outer+1;
      while(m_outerIndex[k]==0)
        m_outerIndex[k++] = 1;
      while (k<=m_outerSize && m_outerIndex[k]!=0)
        m_outerIndex[k++]++;
      p = 0;
      --k;
      k = m_outerIndex[k]-1;
      while (k>0)
      {
        m_data.index(k) = m_data.index(k-1);
        m_data.value(k) = m_data.value(k-1);
        k--;
      }
    }
    else
    {

      Index j = outer+2;
      while (j<=m_outerSize && m_outerIndex[j]!=0)
        m_outerIndex[j++]++;
      --j;
      
      Index k = m_outerIndex[j]-1;
      while (k>=Index(p))
      {
        m_data.index(k) = m_data.index(k-1);
        m_data.value(k) = m_data.value(k-1);
        k--;
      }
    }
  }

  while ( (p > startId) && (m_data.index(p-1) > inner) )
  {
    m_data.index(p) = m_data.index(p-1);
    m_data.value(p) = m_data.value(p-1);
    --p;
  }

  m_data.index(p) = inner;
  return (m_data.value(p) = Scalar(0));
}

namespace internal {

template<typename Scalar_, int Options_, typename StorageIndex_>
struct evaluator<SparseMatrix<Scalar_,Options_,StorageIndex_> >
  : evaluator<SparseCompressedBase<SparseMatrix<Scalar_,Options_,StorageIndex_> > >
{
  typedef evaluator<SparseCompressedBase<SparseMatrix<Scalar_,Options_,StorageIndex_> > > Base;
  typedef SparseMatrix<Scalar_,Options_,StorageIndex_> SparseMatrixType;
  evaluator() : Base() {}
  explicit evaluator(const SparseMatrixType &mat) : Base(mat) {}
};

}

} 

#endif 
