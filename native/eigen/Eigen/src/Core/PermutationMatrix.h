#ifndef EIGEN_PERMUTATIONMATRIX_H
#define EIGEN_PERMUTATIONMATRIX_H

namespace Eigen { 

namespace internal {

enum PermPermProduct_t {PermPermProduct};

} 

template<typename Derived>
class PermutationBase : public EigenBase<Derived>
{
    typedef internal::traits<Derived> Traits;
    typedef EigenBase<Derived> Base;
  public:

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    typedef typename Traits::IndicesType IndicesType;
    enum {
      Flags = Traits::Flags,
      RowsAtCompileTime = Traits::RowsAtCompileTime,
      ColsAtCompileTime = Traits::ColsAtCompileTime,
      MaxRowsAtCompileTime = Traits::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = Traits::MaxColsAtCompileTime
    };
    typedef typename Traits::StorageIndex StorageIndex;
    typedef Matrix<StorageIndex,RowsAtCompileTime,ColsAtCompileTime,0,MaxRowsAtCompileTime,MaxColsAtCompileTime>
            DenseMatrixType;
    typedef PermutationMatrix<IndicesType::SizeAtCompileTime,IndicesType::MaxSizeAtCompileTime,StorageIndex>
            PlainPermutationType;
    typedef PlainPermutationType PlainObject;
    using Base::derived;
    typedef Inverse<Derived> InverseReturnType;
    typedef void Scalar;
    #endif

    template<typename OtherDerived>
    Derived& operator=(const PermutationBase<OtherDerived>& other)
    {
      indices() = other.indices();
      return derived();
    }

    template<typename OtherDerived>
    Derived& operator=(const TranspositionsBase<OtherDerived>& tr)
    {
      setIdentity(tr.size());
      for(Index k=size()-1; k>=0; --k)
        applyTranspositionOnTheRight(k,tr.coeff(k));
      return derived();
    }

    inline EIGEN_DEVICE_FUNC Index rows() const { return Index(indices().size()); }

    inline EIGEN_DEVICE_FUNC Index cols() const { return Index(indices().size()); }

    inline EIGEN_DEVICE_FUNC Index size() const { return Index(indices().size()); }

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    template<typename DenseDerived>
    void evalTo(MatrixBase<DenseDerived>& other) const
    {
      other.setZero();
      for (Index i=0; i<rows(); ++i)
        other.coeffRef(indices().coeff(i),i) = typename DenseDerived::Scalar(1);
    }
    #endif

    DenseMatrixType toDenseMatrix() const
    {
      return derived();
    }

    const IndicesType& indices() const { return derived().indices(); }
    
    IndicesType& indices() { return derived().indices(); }

    inline void resize(Index newSize)
    {
      indices().resize(newSize);
    }

    void setIdentity()
    {
      StorageIndex n = StorageIndex(size());
      for(StorageIndex i = 0; i < n; ++i)
        indices().coeffRef(i) = i;
    }

    void setIdentity(Index newSize)
    {
      resize(newSize);
      setIdentity();
    }

    Derived& applyTranspositionOnTheLeft(Index i, Index j)
    {
      eigen_assert(i>=0 && j>=0 && i<size() && j<size());
      for(Index k = 0; k < size(); ++k)
      {
        if(indices().coeff(k) == i) indices().coeffRef(k) = StorageIndex(j);
        else if(indices().coeff(k) == j) indices().coeffRef(k) = StorageIndex(i);
      }
      return derived();
    }

    Derived& applyTranspositionOnTheRight(Index i, Index j)
    {
      eigen_assert(i>=0 && j>=0 && i<size() && j<size());
      std::swap(indices().coeffRef(i), indices().coeffRef(j));
      return derived();
    }

    inline InverseReturnType inverse() const
    { return InverseReturnType(derived()); }

    inline InverseReturnType transpose() const
    { return InverseReturnType(derived()); }

#ifndef EIGEN_PARSED_BY_DOXYGEN
  protected:
    template<typename OtherDerived>
    void assignTranspose(const PermutationBase<OtherDerived>& other)
    {
      for (Index i=0; i<rows();++i) indices().coeffRef(other.indices().coeff(i)) = i;
    }
    template<typename Lhs,typename Rhs>
    void assignProduct(const Lhs& lhs, const Rhs& rhs)
    {
      eigen_assert(lhs.cols() == rhs.rows());
      for (Index i=0; i<rows();++i) indices().coeffRef(i) = lhs.indices().coeff(rhs.indices().coeff(i));
    }
#endif

  public:

    template<typename Other>
    inline PlainPermutationType operator*(const PermutationBase<Other>& other) const
    { return PlainPermutationType(internal::PermPermProduct, derived(), other.derived()); }

    template<typename Other>
    inline PlainPermutationType operator*(const InverseImpl<Other,PermutationStorage>& other) const
    { return PlainPermutationType(internal::PermPermProduct, *this, other.eval()); }

    template<typename Other> friend
    inline PlainPermutationType operator*(const InverseImpl<Other, PermutationStorage>& other, const PermutationBase& perm)
    { return PlainPermutationType(internal::PermPermProduct, other.eval(), perm); }

    Index determinant() const
    {
      Index res = 1;
      Index n = size();
      Matrix<bool,RowsAtCompileTime,1,0,MaxRowsAtCompileTime> mask(n);
      mask.fill(false);
      Index r = 0;
      while(r < n)
      {
        
        while(r<n && mask[r]) r++;
        if(r>=n)
          break;
        
        Index k0 = r++;
        mask.coeffRef(k0) = true;
        for(Index k=indices().coeff(k0); k!=k0; k=indices().coeff(k))
        {
          mask.coeffRef(k) = true;
          res = -res;
        }
      }
      return res;
    }

  protected:

};

namespace internal {
template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename StorageIndex_>
struct traits<PermutationMatrix<SizeAtCompileTime, MaxSizeAtCompileTime, StorageIndex_> >
 : traits<Matrix<StorageIndex_,SizeAtCompileTime,SizeAtCompileTime,0,MaxSizeAtCompileTime,MaxSizeAtCompileTime> >
{
  typedef PermutationStorage StorageKind;
  typedef Matrix<StorageIndex_, SizeAtCompileTime, 1, 0, MaxSizeAtCompileTime, 1> IndicesType;
  typedef StorageIndex_ StorageIndex;
  typedef void Scalar;
};
}

template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename StorageIndex_>
class PermutationMatrix : public PermutationBase<PermutationMatrix<SizeAtCompileTime, MaxSizeAtCompileTime, StorageIndex_> >
{
    typedef PermutationBase<PermutationMatrix> Base;
    typedef internal::traits<PermutationMatrix> Traits;
  public:

    typedef const PermutationMatrix& Nested;

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    typedef typename Traits::IndicesType IndicesType;
    typedef typename Traits::StorageIndex StorageIndex;
    #endif

    inline PermutationMatrix()
    {}

    explicit inline PermutationMatrix(Index size) : m_indices(size)
    {
      eigen_internal_assert(size <= NumTraits<StorageIndex>::highest());
    }

    template<typename OtherDerived>
    inline PermutationMatrix(const PermutationBase<OtherDerived>& other)
      : m_indices(other.indices()) {}

    template<typename Other>
    explicit inline PermutationMatrix(const MatrixBase<Other>& indices) : m_indices(indices)
    {}

    template<typename Other>
    explicit PermutationMatrix(const TranspositionsBase<Other>& tr)
      : m_indices(tr.size())
    {
      *this = tr;
    }

    template<typename Other>
    PermutationMatrix& operator=(const PermutationBase<Other>& other)
    {
      m_indices = other.indices();
      return *this;
    }

    template<typename Other>
    PermutationMatrix& operator=(const TranspositionsBase<Other>& tr)
    {
      return Base::operator=(tr.derived());
    }

    const IndicesType& indices() const { return m_indices; }
    
    IndicesType& indices() { return m_indices; }

#ifndef EIGEN_PARSED_BY_DOXYGEN
    template<typename Other>
    PermutationMatrix(const InverseImpl<Other,PermutationStorage>& other)
      : m_indices(other.derived().nestedExpression().size())
    {
      eigen_internal_assert(m_indices.size() <= NumTraits<StorageIndex>::highest());
      StorageIndex end = StorageIndex(m_indices.size());
      for (StorageIndex i=0; i<end;++i)
        m_indices.coeffRef(other.derived().nestedExpression().indices().coeff(i)) = i;
    }
    template<typename Lhs,typename Rhs>
    PermutationMatrix(internal::PermPermProduct_t, const Lhs& lhs, const Rhs& rhs)
      : m_indices(lhs.indices().size())
    {
      Base::assignProduct(lhs,rhs);
    }
#endif

  protected:

    IndicesType m_indices;
};

namespace internal {
template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename StorageIndex_, int _PacketAccess>
struct traits<Map<PermutationMatrix<SizeAtCompileTime, MaxSizeAtCompileTime, StorageIndex_>,_PacketAccess> >
 : traits<Matrix<StorageIndex_,SizeAtCompileTime,SizeAtCompileTime,0,MaxSizeAtCompileTime,MaxSizeAtCompileTime> >
{
  typedef PermutationStorage StorageKind;
  typedef Map<const Matrix<StorageIndex_, SizeAtCompileTime, 1, 0, MaxSizeAtCompileTime, 1>, _PacketAccess> IndicesType;
  typedef StorageIndex_ StorageIndex;
  typedef void Scalar;
};
}

template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename StorageIndex_, int _PacketAccess>
class Map<PermutationMatrix<SizeAtCompileTime, MaxSizeAtCompileTime, StorageIndex_>,_PacketAccess>
  : public PermutationBase<Map<PermutationMatrix<SizeAtCompileTime, MaxSizeAtCompileTime, StorageIndex_>,_PacketAccess> >
{
    typedef PermutationBase<Map> Base;
    typedef internal::traits<Map> Traits;
  public:

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    typedef typename Traits::IndicesType IndicesType;
    typedef typename IndicesType::Scalar StorageIndex;
    #endif

    inline Map(const StorageIndex* indicesPtr)
      : m_indices(indicesPtr)
    {}

    inline Map(const StorageIndex* indicesPtr, Index size)
      : m_indices(indicesPtr,size)
    {}

    template<typename Other>
    Map& operator=(const PermutationBase<Other>& other)
    { return Base::operator=(other.derived()); }

    template<typename Other>
    Map& operator=(const TranspositionsBase<Other>& tr)
    { return Base::operator=(tr.derived()); }

    #ifndef EIGEN_PARSED_BY_DOXYGEN

    Map& operator=(const Map& other)
    {
      m_indices = other.m_indices;
      return *this;
    }
    #endif

    const IndicesType& indices() const { return m_indices; }
    
    IndicesType& indices() { return m_indices; }

  protected:

    IndicesType m_indices;
};

template<typename IndicesType_> class TranspositionsWrapper;
namespace internal {
template<typename IndicesType_>
struct traits<PermutationWrapper<IndicesType_> >
{
  typedef PermutationStorage StorageKind;
  typedef void Scalar;
  typedef typename IndicesType_::Scalar StorageIndex;
  typedef IndicesType_ IndicesType;
  enum {
    RowsAtCompileTime = IndicesType_::SizeAtCompileTime,
    ColsAtCompileTime = IndicesType_::SizeAtCompileTime,
    MaxRowsAtCompileTime = IndicesType::MaxSizeAtCompileTime,
    MaxColsAtCompileTime = IndicesType::MaxSizeAtCompileTime,
    Flags = 0
  };
};
}

template<typename IndicesType_>
class PermutationWrapper : public PermutationBase<PermutationWrapper<IndicesType_> >
{
    typedef PermutationBase<PermutationWrapper> Base;
    typedef internal::traits<PermutationWrapper> Traits;
  public:

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    typedef typename Traits::IndicesType IndicesType;
    #endif

    inline PermutationWrapper(const IndicesType& indices)
      : m_indices(indices)
    {}

    const typename internal::remove_all<typename IndicesType::Nested>::type&
    indices() const { return m_indices; }

  protected:

    typename IndicesType::Nested m_indices;
};

template<typename MatrixDerived, typename PermutationDerived>
EIGEN_DEVICE_FUNC
const Product<MatrixDerived, PermutationDerived, AliasFreeProduct>
operator*(const MatrixBase<MatrixDerived> &matrix,
          const PermutationBase<PermutationDerived>& permutation)
{
  return Product<MatrixDerived, PermutationDerived, AliasFreeProduct>
            (matrix.derived(), permutation.derived());
}

template<typename PermutationDerived, typename MatrixDerived>
EIGEN_DEVICE_FUNC
const Product<PermutationDerived, MatrixDerived, AliasFreeProduct>
operator*(const PermutationBase<PermutationDerived> &permutation,
          const MatrixBase<MatrixDerived>& matrix)
{
  return Product<PermutationDerived, MatrixDerived, AliasFreeProduct>
            (permutation.derived(), matrix.derived());
}

template<typename PermutationType>
class InverseImpl<PermutationType, PermutationStorage>
  : public EigenBase<Inverse<PermutationType> >
{
    typedef typename PermutationType::PlainPermutationType PlainPermutationType;
    typedef internal::traits<PermutationType> PermTraits;
  protected:
    InverseImpl() {}
  public:
    typedef Inverse<PermutationType> InverseType;
    using EigenBase<Inverse<PermutationType> >::derived;

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    typedef typename PermutationType::DenseMatrixType DenseMatrixType;
    enum {
      RowsAtCompileTime = PermTraits::RowsAtCompileTime,
      ColsAtCompileTime = PermTraits::ColsAtCompileTime,
      MaxRowsAtCompileTime = PermTraits::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = PermTraits::MaxColsAtCompileTime
    };
    #endif

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    template<typename DenseDerived>
    void evalTo(MatrixBase<DenseDerived>& other) const
    {
      other.setZero();
      for (Index i=0; i<derived().rows();++i)
        other.coeffRef(i, derived().nestedExpression().indices().coeff(i)) = typename DenseDerived::Scalar(1);
    }
    #endif

    PlainPermutationType eval() const { return derived(); }

    DenseMatrixType toDenseMatrix() const { return derived(); }

    template<typename OtherDerived> friend
    const Product<OtherDerived, InverseType, AliasFreeProduct>
    operator*(const MatrixBase<OtherDerived>& matrix, const InverseType& trPerm)
    {
      return Product<OtherDerived, InverseType, AliasFreeProduct>(matrix.derived(), trPerm.derived());
    }

    template<typename OtherDerived>
    const Product<InverseType, OtherDerived, AliasFreeProduct>
    operator*(const MatrixBase<OtherDerived>& matrix) const
    {
      return Product<InverseType, OtherDerived, AliasFreeProduct>(derived(), matrix.derived());
    }
};

template<typename Derived>
const PermutationWrapper<const Derived> MatrixBase<Derived>::asPermutation() const
{
  return derived();
}

namespace internal {

template<> struct AssignmentKind<DenseShape,PermutationShape> { typedef EigenBase2EigenBase Kind; };

} 

} 

#endif 
