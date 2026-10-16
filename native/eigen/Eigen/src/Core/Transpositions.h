#ifndef EIGEN_TRANSPOSITIONS_H
#define EIGEN_TRANSPOSITIONS_H

namespace Eigen {

template<typename Derived>
class TranspositionsBase
{
    typedef internal::traits<Derived> Traits;

  public:

    typedef typename Traits::IndicesType IndicesType;
    typedef typename IndicesType::Scalar StorageIndex;
    typedef Eigen::Index Index; 

    EIGEN_DEVICE_FUNC
    Derived& derived() { return *static_cast<Derived*>(this); }
    EIGEN_DEVICE_FUNC
    const Derived& derived() const { return *static_cast<const Derived*>(this); }

    template<typename OtherDerived>
    Derived& operator=(const TranspositionsBase<OtherDerived>& other)
    {
      indices() = other.indices();
      return derived();
    }

    EIGEN_DEVICE_FUNC
    Index size() const { return indices().size(); }
    
    EIGEN_DEVICE_FUNC
    Index rows() const { return indices().size(); }
    
    EIGEN_DEVICE_FUNC
    Index cols() const { return indices().size(); }

    EIGEN_DEVICE_FUNC
    inline const StorageIndex& coeff(Index i) const { return indices().coeff(i); }
    
    inline StorageIndex& coeffRef(Index i) { return indices().coeffRef(i); }
    
    inline const StorageIndex& operator()(Index i) const { return indices()(i); }
    
    inline StorageIndex& operator()(Index i) { return indices()(i); }
    
    inline const StorageIndex& operator[](Index i) const { return indices()(i); }
    
    inline StorageIndex& operator[](Index i) { return indices()(i); }

    EIGEN_DEVICE_FUNC
    const IndicesType& indices() const { return derived().indices(); }
    
    EIGEN_DEVICE_FUNC
    IndicesType& indices() { return derived().indices(); }

    inline void resize(Index newSize)
    {
      indices().resize(newSize);
    }

    void setIdentity()
    {
      for(StorageIndex i = 0; i < indices().size(); ++i)
        coeffRef(i) = i;
    }

    inline Transpose<TranspositionsBase> inverse() const
    { return Transpose<TranspositionsBase>(derived()); }

    inline Transpose<TranspositionsBase> transpose() const
    { return Transpose<TranspositionsBase>(derived()); }

  protected:
};

namespace internal {
template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename StorageIndex_>
struct traits<Transpositions<SizeAtCompileTime,MaxSizeAtCompileTime,StorageIndex_> >
 : traits<PermutationMatrix<SizeAtCompileTime,MaxSizeAtCompileTime,StorageIndex_> >
{
  typedef Matrix<StorageIndex_, SizeAtCompileTime, 1, 0, MaxSizeAtCompileTime, 1> IndicesType;
  typedef TranspositionsStorage StorageKind;
};
}

template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename StorageIndex_>
class Transpositions : public TranspositionsBase<Transpositions<SizeAtCompileTime,MaxSizeAtCompileTime,StorageIndex_> >
{
    typedef internal::traits<Transpositions> Traits;
  public:

    typedef TranspositionsBase<Transpositions> Base;
    typedef typename Traits::IndicesType IndicesType;
    typedef typename IndicesType::Scalar StorageIndex;

    inline Transpositions() {}

    template<typename OtherDerived>
    inline Transpositions(const TranspositionsBase<OtherDerived>& other)
      : m_indices(other.indices()) {}

    template<typename Other>
    explicit inline Transpositions(const MatrixBase<Other>& indices) : m_indices(indices)
    {}

    template<typename OtherDerived>
    Transpositions& operator=(const TranspositionsBase<OtherDerived>& other)
    {
      return Base::operator=(other);
    }

    inline Transpositions(Index size) : m_indices(size)
    {}

    EIGEN_DEVICE_FUNC
    const IndicesType& indices() const { return m_indices; }
    
    EIGEN_DEVICE_FUNC
    IndicesType& indices() { return m_indices; }

  protected:

    IndicesType m_indices;
};

namespace internal {
template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename StorageIndex_, int _PacketAccess>
struct traits<Map<Transpositions<SizeAtCompileTime,MaxSizeAtCompileTime,StorageIndex_>,_PacketAccess> >
 : traits<PermutationMatrix<SizeAtCompileTime,MaxSizeAtCompileTime,StorageIndex_> >
{
  typedef Map<const Matrix<StorageIndex_,SizeAtCompileTime,1,0,MaxSizeAtCompileTime,1>, _PacketAccess> IndicesType;
  typedef StorageIndex_ StorageIndex;
  typedef TranspositionsStorage StorageKind;
};
}

template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename StorageIndex_, int PacketAccess>
class Map<Transpositions<SizeAtCompileTime,MaxSizeAtCompileTime,StorageIndex_>,PacketAccess>
 : public TranspositionsBase<Map<Transpositions<SizeAtCompileTime,MaxSizeAtCompileTime,StorageIndex_>,PacketAccess> >
{
    typedef internal::traits<Map> Traits;
  public:

    typedef TranspositionsBase<Map> Base;
    typedef typename Traits::IndicesType IndicesType;
    typedef typename IndicesType::Scalar StorageIndex;

    explicit inline Map(const StorageIndex* indicesPtr)
      : m_indices(indicesPtr)
    {}

    inline Map(const StorageIndex* indicesPtr, Index size)
      : m_indices(indicesPtr,size)
    {}

    template<typename OtherDerived>
    Map& operator=(const TranspositionsBase<OtherDerived>& other)
    {
      return Base::operator=(other);
    }

    #ifndef EIGEN_PARSED_BY_DOXYGEN

    Map& operator=(const Map& other)
    {
      m_indices = other.m_indices;
      return *this;
    }
    #endif

    EIGEN_DEVICE_FUNC
    const IndicesType& indices() const { return m_indices; }

    EIGEN_DEVICE_FUNC
    IndicesType& indices() { return m_indices; }

  protected:

    IndicesType m_indices;
};

namespace internal {
template<typename IndicesType_>
struct traits<TranspositionsWrapper<IndicesType_> >
 : traits<PermutationWrapper<IndicesType_> >
{
  typedef TranspositionsStorage StorageKind;
};
}

template<typename IndicesType_>
class TranspositionsWrapper
 : public TranspositionsBase<TranspositionsWrapper<IndicesType_> >
{
    typedef internal::traits<TranspositionsWrapper> Traits;
  public:

    typedef TranspositionsBase<TranspositionsWrapper> Base;
    typedef typename Traits::IndicesType IndicesType;
    typedef typename IndicesType::Scalar StorageIndex;

    explicit inline TranspositionsWrapper(IndicesType& indices)
      : m_indices(indices)
    {}

    template<typename OtherDerived>
    TranspositionsWrapper& operator=(const TranspositionsBase<OtherDerived>& other)
    {
      return Base::operator=(other);
    }

    EIGEN_DEVICE_FUNC
    const IndicesType& indices() const { return m_indices; }

    EIGEN_DEVICE_FUNC
    IndicesType& indices() { return m_indices; }

  protected:

    typename IndicesType::Nested m_indices;
};

template<typename MatrixDerived, typename TranspositionsDerived>
EIGEN_DEVICE_FUNC
const Product<MatrixDerived, TranspositionsDerived, AliasFreeProduct>
operator*(const MatrixBase<MatrixDerived> &matrix,
          const TranspositionsBase<TranspositionsDerived>& transpositions)
{
  return Product<MatrixDerived, TranspositionsDerived, AliasFreeProduct>
            (matrix.derived(), transpositions.derived());
}

template<typename TranspositionsDerived, typename MatrixDerived>
EIGEN_DEVICE_FUNC
const Product<TranspositionsDerived, MatrixDerived, AliasFreeProduct>
operator*(const TranspositionsBase<TranspositionsDerived> &transpositions,
          const MatrixBase<MatrixDerived>& matrix)
{
  return Product<TranspositionsDerived, MatrixDerived, AliasFreeProduct>
            (transpositions.derived(), matrix.derived());
}

namespace internal {

template<typename Derived>
struct traits<Transpose<TranspositionsBase<Derived> > >
 : traits<Derived>
{};

} 

template<typename TranspositionsDerived>
class Transpose<TranspositionsBase<TranspositionsDerived> >
{
    typedef TranspositionsDerived TranspositionType;
    typedef typename TranspositionType::IndicesType IndicesType;
  public:

    explicit Transpose(const TranspositionType& t) : m_transpositions(t) {}

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    Index size() const EIGEN_NOEXCEPT { return m_transpositions.size(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    Index rows() const EIGEN_NOEXCEPT { return m_transpositions.size(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    Index cols() const EIGEN_NOEXCEPT { return m_transpositions.size(); }

    template<typename OtherDerived> friend
    const Product<OtherDerived, Transpose, AliasFreeProduct>
    operator*(const MatrixBase<OtherDerived>& matrix, const Transpose& trt)
    {
      return Product<OtherDerived, Transpose, AliasFreeProduct>(matrix.derived(), trt);
    }

    template<typename OtherDerived>
    const Product<Transpose, OtherDerived, AliasFreeProduct>
    operator*(const MatrixBase<OtherDerived>& matrix) const
    {
      return Product<Transpose, OtherDerived, AliasFreeProduct>(*this, matrix.derived());
    }

    EIGEN_DEVICE_FUNC
    const TranspositionType& nestedExpression() const { return m_transpositions; }

  protected:
    const TranspositionType& m_transpositions;
};

} 

#endif 
