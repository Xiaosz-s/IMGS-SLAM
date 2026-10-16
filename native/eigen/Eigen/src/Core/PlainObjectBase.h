#ifndef EIGEN_DENSESTORAGEBASE_H
#define EIGEN_DENSESTORAGEBASE_H

#if defined(EIGEN_INITIALIZE_MATRICES_BY_ZERO)
# define EIGEN_INITIALIZE_COEFFS
# define EIGEN_INITIALIZE_COEFFS_IF_THAT_OPTION_IS_ENABLED for(Index i=0;i<base().size();++i) coeffRef(i)=Scalar(0);
#elif defined(EIGEN_INITIALIZE_MATRICES_BY_NAN)
# define EIGEN_INITIALIZE_COEFFS
# define EIGEN_INITIALIZE_COEFFS_IF_THAT_OPTION_IS_ENABLED for(Index i=0;i<base().size();++i) coeffRef(i)=std::numeric_limits<Scalar>::quiet_NaN();
#else
# undef EIGEN_INITIALIZE_COEFFS
# define EIGEN_INITIALIZE_COEFFS_IF_THAT_OPTION_IS_ENABLED
#endif

namespace Eigen {

namespace internal {

template<int MaxSizeAtCompileTime> struct check_rows_cols_for_overflow {
  template<typename Index>
  EIGEN_DEVICE_FUNC
  static EIGEN_ALWAYS_INLINE void run(Index, Index)
  {
  }
};

template<> struct check_rows_cols_for_overflow<Dynamic> {
  template<typename Index>
  EIGEN_DEVICE_FUNC
  static EIGEN_ALWAYS_INLINE void run(Index rows, Index cols)
  {

    Index max_index = (std::size_t(1) << (8 * sizeof(Index) - 1)) - 1; 
    bool error = (rows == 0 || cols == 0) ? false
               : (rows > max_index / cols);
    if (error)
      throw_std_bad_alloc();
  }
};

template <typename Derived,
          typename OtherDerived = Derived,
          bool IsVector = bool(Derived::IsVectorAtCompileTime) && bool(OtherDerived::IsVectorAtCompileTime)>
struct conservative_resize_like_impl;

template<typename MatrixTypeA, typename MatrixTypeB, bool SwapPointers> struct matrix_swap_impl;

} 

#ifdef EIGEN_PARSED_BY_DOXYGEN
namespace doxygen {

template<typename Derived> struct dense_xpr_base_dispatcher;

template<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
struct dense_xpr_base_dispatcher<Matrix<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_> >
    : public MatrixBase {};

template<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
struct dense_xpr_base_dispatcher<Array<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_> >
    : public ArrayBase {};

} 

template<typename Derived>
class PlainObjectBase : public doxygen::dense_xpr_base_dispatcher<Derived>
#else
template<typename Derived>
class PlainObjectBase : public internal::dense_xpr_base<Derived>::type
#endif
{
  public:
    enum { Options = internal::traits<Derived>::Options };
    typedef typename internal::dense_xpr_base<Derived>::type Base;

    typedef typename internal::traits<Derived>::StorageKind StorageKind;
    typedef typename internal::traits<Derived>::Scalar Scalar;

    typedef typename internal::packet_traits<Scalar>::type PacketScalar;
    typedef typename NumTraits<Scalar>::Real RealScalar;
    typedef Derived DenseType;

    using Base::RowsAtCompileTime;
    using Base::ColsAtCompileTime;
    using Base::SizeAtCompileTime;
    using Base::MaxRowsAtCompileTime;
    using Base::MaxColsAtCompileTime;
    using Base::MaxSizeAtCompileTime;
    using Base::IsVectorAtCompileTime;
    using Base::Flags;

    typedef Eigen::Map<Derived, Unaligned>  MapType;
    typedef const Eigen::Map<const Derived, Unaligned> ConstMapType;
    typedef Eigen::Map<Derived, AlignedMax> AlignedMapType;
    typedef const Eigen::Map<const Derived, AlignedMax> ConstAlignedMapType;
    template<typename StrideType> struct StridedMapType { typedef Eigen::Map<Derived, Unaligned, StrideType> type; };
    template<typename StrideType> struct StridedConstMapType { typedef Eigen::Map<const Derived, Unaligned, StrideType> type; };
    template<typename StrideType> struct StridedAlignedMapType { typedef Eigen::Map<Derived, AlignedMax, StrideType> type; };
    template<typename StrideType> struct StridedConstAlignedMapType { typedef Eigen::Map<const Derived, AlignedMax, StrideType> type; };

  protected:
    DenseStorage<Scalar, Base::MaxSizeAtCompileTime, Base::RowsAtCompileTime, Base::ColsAtCompileTime, Options> m_storage;

  public:
    enum { NeedsToAlign = (SizeAtCompileTime != Dynamic) && (internal::traits<Derived>::Alignment>0) };
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

    EIGEN_DEVICE_FUNC
    Base& base() { return *static_cast<Base*>(this); }
    EIGEN_DEVICE_FUNC
    const Base& base() const { return *static_cast<const Base*>(this); }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index rows() const EIGEN_NOEXCEPT { return m_storage.rows(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    Index cols() const EIGEN_NOEXCEPT { return m_storage.cols(); }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE const Scalar& coeff(Index rowId, Index colId) const
    {
      if(Flags & RowMajorBit)
        return m_storage.data()[colId + rowId * m_storage.cols()];
      else 
        return m_storage.data()[rowId + colId * m_storage.rows()];
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE const Scalar& coeff(Index index) const
    {
      return m_storage.data()[index];
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar& coeffRef(Index rowId, Index colId)
    {
      if(Flags & RowMajorBit)
        return m_storage.data()[colId + rowId * m_storage.cols()];
      else 
        return m_storage.data()[rowId + colId * m_storage.rows()];
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar& coeffRef(Index index)
    {
      return m_storage.data()[index];
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE const Scalar& coeffRef(Index rowId, Index colId) const
    {
      if(Flags & RowMajorBit)
        return m_storage.data()[colId + rowId * m_storage.cols()];
      else 
        return m_storage.data()[rowId + colId * m_storage.rows()];
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE const Scalar& coeffRef(Index index) const
    {
      return m_storage.data()[index];
    }

    template<int LoadMode>
    EIGEN_STRONG_INLINE PacketScalar packet(Index rowId, Index colId) const
    {
      return internal::ploadt<PacketScalar, LoadMode>
               (m_storage.data() + (Flags & RowMajorBit
                                   ? colId + rowId * m_storage.cols()
                                   : rowId + colId * m_storage.rows()));
    }

    template<int LoadMode>
    EIGEN_STRONG_INLINE PacketScalar packet(Index index) const
    {
      return internal::ploadt<PacketScalar, LoadMode>(m_storage.data() + index);
    }

    template<int StoreMode>
    EIGEN_STRONG_INLINE void writePacket(Index rowId, Index colId, const PacketScalar& val)
    {
      internal::pstoret<Scalar, PacketScalar, StoreMode>
              (m_storage.data() + (Flags & RowMajorBit
                                   ? colId + rowId * m_storage.cols()
                                   : rowId + colId * m_storage.rows()), val);
    }

    template<int StoreMode>
    EIGEN_STRONG_INLINE void writePacket(Index index, const PacketScalar& val)
    {
      internal::pstoret<Scalar, PacketScalar, StoreMode>(m_storage.data() + index, val);
    }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar *data() const
    { return m_storage.data(); }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Scalar *data()
    { return m_storage.data(); }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void resize(Index rows, Index cols)
    {
      eigen_assert(   EIGEN_IMPLIES(RowsAtCompileTime!=Dynamic,rows==RowsAtCompileTime)
                   && EIGEN_IMPLIES(ColsAtCompileTime!=Dynamic,cols==ColsAtCompileTime)
                   && EIGEN_IMPLIES(RowsAtCompileTime==Dynamic && MaxRowsAtCompileTime!=Dynamic,rows<=MaxRowsAtCompileTime)
                   && EIGEN_IMPLIES(ColsAtCompileTime==Dynamic && MaxColsAtCompileTime!=Dynamic,cols<=MaxColsAtCompileTime)
                   && rows>=0 && cols>=0 && "Invalid sizes when resizing a matrix or array.");
      internal::check_rows_cols_for_overflow<MaxSizeAtCompileTime>::run(rows, cols);
      #ifdef EIGEN_INITIALIZE_COEFFS
        Index size = rows*cols;
        bool size_changed = size != this->size();
        m_storage.resize(size, rows, cols);
        if(size_changed) EIGEN_INITIALIZE_COEFFS_IF_THAT_OPTION_IS_ENABLED
      #else
        m_storage.resize(rows*cols, rows, cols);
      #endif
    }

    EIGEN_DEVICE_FUNC
    inline void resize(Index size)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(PlainObjectBase)
      eigen_assert(((SizeAtCompileTime == Dynamic && (MaxSizeAtCompileTime==Dynamic || size<=MaxSizeAtCompileTime)) || SizeAtCompileTime == size) && size>=0);
      #ifdef EIGEN_INITIALIZE_COEFFS
        bool size_changed = size != this->size();
      #endif
      if(RowsAtCompileTime == 1)
        m_storage.resize(size, 1, size);
      else
        m_storage.resize(size, size, 1);
      #ifdef EIGEN_INITIALIZE_COEFFS
        if(size_changed) EIGEN_INITIALIZE_COEFFS_IF_THAT_OPTION_IS_ENABLED
      #endif
    }

    EIGEN_DEVICE_FUNC
    inline void resize(NoChange_t, Index cols)
    {
      resize(rows(), cols);
    }

    EIGEN_DEVICE_FUNC
    inline void resize(Index rows, NoChange_t)
    {
      resize(rows, cols());
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void resizeLike(const EigenBase<OtherDerived>& _other)
    {
      const OtherDerived& other = _other.derived();
      internal::check_rows_cols_for_overflow<MaxSizeAtCompileTime>::run(other.rows(), other.cols());
      const Index othersize = other.rows()*other.cols();
      if(RowsAtCompileTime == 1)
      {
        eigen_assert(other.rows() == 1 || other.cols() == 1);
        resize(1, othersize);
      }
      else if(ColsAtCompileTime == 1)
      {
        eigen_assert(other.rows() == 1 || other.cols() == 1);
        resize(othersize, 1);
      }
      else resize(other.rows(), other.cols());
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void conservativeResize(Index rows, Index cols)
    {
      internal::conservative_resize_like_impl<Derived>::run(*this, rows, cols);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void conservativeResize(Index rows, NoChange_t)
    {
      
      conservativeResize(rows, cols());
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void conservativeResize(NoChange_t, Index cols)
    {
      
      conservativeResize(rows(), cols);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void conservativeResize(Index size)
    {
      internal::conservative_resize_like_impl<Derived>::run(*this, size);
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void conservativeResizeLike(const DenseBase<OtherDerived>& other)
    {
      internal::conservative_resize_like_impl<Derived,OtherDerived>::run(*this, other);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Derived& operator=(const PlainObjectBase& other)
    {
      return _set(other);
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Derived& lazyAssign(const DenseBase<OtherDerived>& other)
    {
      _resize_to_match(other);
      return Base::lazyAssign(other.derived());
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Derived& operator=(const ReturnByValue<OtherDerived>& func)
    {
      resize(func.rows(), func.cols());
      return Base::operator=(func);
    }

  protected:

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE PlainObjectBase() : m_storage()
    {

    }

#ifndef EIGEN_PARSED_BY_DOXYGEN

    EIGEN_DEVICE_FUNC
    explicit PlainObjectBase(internal::constructor_without_unaligned_array_assert)
      : m_storage(internal::constructor_without_unaligned_array_assert())
    {

    }
#endif

#if EIGEN_HAS_RVALUE_REFERENCES
    EIGEN_DEVICE_FUNC
    PlainObjectBase(PlainObjectBase&& other) EIGEN_NOEXCEPT
      : m_storage( std::move(other.m_storage) )
    {
    }

    EIGEN_DEVICE_FUNC
    PlainObjectBase& operator=(PlainObjectBase&& other) EIGEN_NOEXCEPT
    {
      _check_template_params();
      m_storage = std::move(other.m_storage);
      return *this;
    }
#endif

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE PlainObjectBase(const PlainObjectBase& other)
      : Base(), m_storage(other.m_storage) { }
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE PlainObjectBase(Index size, Index rows, Index cols)
      : m_storage(size, rows, cols)
    {

    }

    #if EIGEN_HAS_CXX11

    template <typename... ArgTypes>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    PlainObjectBase(const Scalar& a0, const Scalar& a1, const Scalar& a2,  const Scalar& a3, const ArgTypes&... args)
      : m_storage()
    {
      _check_template_params();
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(PlainObjectBase, sizeof...(args) + 4);
      m_storage.data()[0] = a0;
      m_storage.data()[1] = a1;
      m_storage.data()[2] = a2;
      m_storage.data()[3] = a3;
      Index i = 4;
      auto x = {(m_storage.data()[i++] = args, 0)...};
      static_cast<void>(x);
    }

    EIGEN_DEVICE_FUNC
    explicit EIGEN_STRONG_INLINE PlainObjectBase(const std::initializer_list<std::initializer_list<Scalar>>& list)
      : m_storage()
    {
      _check_template_params();

      size_t list_size = 0;
      if (list.begin() != list.end()) {
        list_size = list.begin()->size();
      }

      if (ColsAtCompileTime == 1 && list.size() == 1) {
        eigen_assert(list_size == static_cast<size_t>(RowsAtCompileTime) || RowsAtCompileTime == Dynamic);
        resize(list_size, ColsAtCompileTime);
        std::copy(list.begin()->begin(), list.begin()->end(), m_storage.data());
      } else {
        eigen_assert(list.size() == static_cast<size_t>(RowsAtCompileTime) || RowsAtCompileTime == Dynamic);
        eigen_assert(list_size == static_cast<size_t>(ColsAtCompileTime) || ColsAtCompileTime == Dynamic);
        resize(list.size(), list_size);

        Index row_index = 0;
        for (const std::initializer_list<Scalar>& row : list) {
          eigen_assert(list_size == row.size());
          Index col_index = 0;
          for (const Scalar& e : row) {
            coeffRef(row_index, col_index) = e;
            ++col_index;
          }
          ++row_index;
        }
      }
    }
    #endif  

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE PlainObjectBase(const DenseBase<OtherDerived> &other)
      : m_storage()
    {
      _check_template_params();
      resizeLike(other);
      _set_noalias(other);
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE PlainObjectBase(const EigenBase<OtherDerived> &other)
      : m_storage()
    {
      _check_template_params();
      resizeLike(other);
      *this = other.derived();
    }
    
    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE PlainObjectBase(const ReturnByValue<OtherDerived>& other)
    {
      _check_template_params();
      
      resize(other.rows(), other.cols());
      other.evalTo(this->derived());
    }

  public:

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Derived& operator=(const EigenBase<OtherDerived> &other)
    {
      _resize_to_match(other);
      Base::operator=(other.derived());
      return this->derived();
    }

    static inline ConstMapType Map(const Scalar* data)
    { return ConstMapType(data); }
    static inline MapType Map(Scalar* data)
    { return MapType(data); }
    static inline ConstMapType Map(const Scalar* data, Index size)
    { return ConstMapType(data, size); }
    static inline MapType Map(Scalar* data, Index size)
    { return MapType(data, size); }
    static inline ConstMapType Map(const Scalar* data, Index rows, Index cols)
    { return ConstMapType(data, rows, cols); }
    static inline MapType Map(Scalar* data, Index rows, Index cols)
    { return MapType(data, rows, cols); }

    static inline ConstAlignedMapType MapAligned(const Scalar* data)
    { return ConstAlignedMapType(data); }
    static inline AlignedMapType MapAligned(Scalar* data)
    { return AlignedMapType(data); }
    static inline ConstAlignedMapType MapAligned(const Scalar* data, Index size)
    { return ConstAlignedMapType(data, size); }
    static inline AlignedMapType MapAligned(Scalar* data, Index size)
    { return AlignedMapType(data, size); }
    static inline ConstAlignedMapType MapAligned(const Scalar* data, Index rows, Index cols)
    { return ConstAlignedMapType(data, rows, cols); }
    static inline AlignedMapType MapAligned(Scalar* data, Index rows, Index cols)
    { return AlignedMapType(data, rows, cols); }

    template<int Outer, int Inner>
    static inline typename StridedConstMapType<Stride<Outer, Inner> >::type Map(const Scalar* data, const Stride<Outer, Inner>& stride)
    { return typename StridedConstMapType<Stride<Outer, Inner> >::type(data, stride); }
    template<int Outer, int Inner>
    static inline typename StridedMapType<Stride<Outer, Inner> >::type Map(Scalar* data, const Stride<Outer, Inner>& stride)
    { return typename StridedMapType<Stride<Outer, Inner> >::type(data, stride); }
    template<int Outer, int Inner>
    static inline typename StridedConstMapType<Stride<Outer, Inner> >::type Map(const Scalar* data, Index size, const Stride<Outer, Inner>& stride)
    { return typename StridedConstMapType<Stride<Outer, Inner> >::type(data, size, stride); }
    template<int Outer, int Inner>
    static inline typename StridedMapType<Stride<Outer, Inner> >::type Map(Scalar* data, Index size, const Stride<Outer, Inner>& stride)
    { return typename StridedMapType<Stride<Outer, Inner> >::type(data, size, stride); }
    template<int Outer, int Inner>
    static inline typename StridedConstMapType<Stride<Outer, Inner> >::type Map(const Scalar* data, Index rows, Index cols, const Stride<Outer, Inner>& stride)
    { return typename StridedConstMapType<Stride<Outer, Inner> >::type(data, rows, cols, stride); }
    template<int Outer, int Inner>
    static inline typename StridedMapType<Stride<Outer, Inner> >::type Map(Scalar* data, Index rows, Index cols, const Stride<Outer, Inner>& stride)
    { return typename StridedMapType<Stride<Outer, Inner> >::type(data, rows, cols, stride); }

    template<int Outer, int Inner>
    static inline typename StridedConstAlignedMapType<Stride<Outer, Inner> >::type MapAligned(const Scalar* data, const Stride<Outer, Inner>& stride)
    { return typename StridedConstAlignedMapType<Stride<Outer, Inner> >::type(data, stride); }
    template<int Outer, int Inner>
    static inline typename StridedAlignedMapType<Stride<Outer, Inner> >::type MapAligned(Scalar* data, const Stride<Outer, Inner>& stride)
    { return typename StridedAlignedMapType<Stride<Outer, Inner> >::type(data, stride); }
    template<int Outer, int Inner>
    static inline typename StridedConstAlignedMapType<Stride<Outer, Inner> >::type MapAligned(const Scalar* data, Index size, const Stride<Outer, Inner>& stride)
    { return typename StridedConstAlignedMapType<Stride<Outer, Inner> >::type(data, size, stride); }
    template<int Outer, int Inner>
    static inline typename StridedAlignedMapType<Stride<Outer, Inner> >::type MapAligned(Scalar* data, Index size, const Stride<Outer, Inner>& stride)
    { return typename StridedAlignedMapType<Stride<Outer, Inner> >::type(data, size, stride); }
    template<int Outer, int Inner>
    static inline typename StridedConstAlignedMapType<Stride<Outer, Inner> >::type MapAligned(const Scalar* data, Index rows, Index cols, const Stride<Outer, Inner>& stride)
    { return typename StridedConstAlignedMapType<Stride<Outer, Inner> >::type(data, rows, cols, stride); }
    template<int Outer, int Inner>
    static inline typename StridedAlignedMapType<Stride<Outer, Inner> >::type MapAligned(Scalar* data, Index rows, Index cols, const Stride<Outer, Inner>& stride)
    { return typename StridedAlignedMapType<Stride<Outer, Inner> >::type(data, rows, cols, stride); }

    using Base::setConstant;
    EIGEN_DEVICE_FUNC Derived& setConstant(Index size, const Scalar& val);
    EIGEN_DEVICE_FUNC Derived& setConstant(Index rows, Index cols, const Scalar& val);
    EIGEN_DEVICE_FUNC Derived& setConstant(NoChange_t, Index cols, const Scalar& val);
    EIGEN_DEVICE_FUNC Derived& setConstant(Index rows, NoChange_t, const Scalar& val);

    using Base::setZero;
    EIGEN_DEVICE_FUNC Derived& setZero(Index size);
    EIGEN_DEVICE_FUNC Derived& setZero(Index rows, Index cols);
    EIGEN_DEVICE_FUNC Derived& setZero(NoChange_t, Index cols);
    EIGEN_DEVICE_FUNC Derived& setZero(Index rows, NoChange_t);

    using Base::setOnes;
    EIGEN_DEVICE_FUNC Derived& setOnes(Index size);
    EIGEN_DEVICE_FUNC Derived& setOnes(Index rows, Index cols);
    EIGEN_DEVICE_FUNC Derived& setOnes(NoChange_t, Index cols);
    EIGEN_DEVICE_FUNC Derived& setOnes(Index rows, NoChange_t);

    using Base::setRandom;
    Derived& setRandom(Index size);
    Derived& setRandom(Index rows, Index cols);
    Derived& setRandom(NoChange_t, Index cols);
    Derived& setRandom(Index rows, NoChange_t);

    #ifdef EIGEN_PLAINOBJECTBASE_PLUGIN
    #include EIGEN_PLAINOBJECTBASE_PLUGIN
    #endif

  protected:

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _resize_to_match(const EigenBase<OtherDerived>& other)
    {
      #ifdef EIGEN_NO_AUTOMATIC_RESIZING
      eigen_assert((this->size()==0 || (IsVectorAtCompileTime ? (this->size() == other.size())
                 : (rows() == other.rows() && cols() == other.cols())))
        && "Size mismatch. Automatic resizing is disabled because EIGEN_NO_AUTOMATIC_RESIZING is defined");
      EIGEN_ONLY_USED_FOR_DEBUG(other);
      #else
      resizeLike(other);
      #endif
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Derived& _set(const DenseBase<OtherDerived>& other)
    {
      internal::call_assignment(this->derived(), other.derived());
      return this->derived();
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Derived& _set_noalias(const DenseBase<OtherDerived>& other)
    {

      internal::call_assignment_no_alias(this->derived(), other.derived(), internal::assign_op<Scalar,typename OtherDerived::Scalar>());
      return this->derived();
    }

    template<typename T0, typename T1>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init2(Index rows, Index cols, typename internal::enable_if<Base::SizeAtCompileTime!=2,T0>::type* = 0)
    {
      const bool t0_is_integer_alike = internal::is_valid_index_type<T0>::value;
      const bool t1_is_integer_alike = internal::is_valid_index_type<T1>::value;
      EIGEN_STATIC_ASSERT(t0_is_integer_alike &&
                          t1_is_integer_alike,
                          FLOATING_POINT_ARGUMENT_PASSED__INTEGER_WAS_EXPECTED)
      resize(rows,cols);
    }

    template<typename T0, typename T1>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init2(const T0& val0, const T1& val1, typename internal::enable_if<Base::SizeAtCompileTime==2,T0>::type* = 0)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(PlainObjectBase, 2)
      m_storage.data()[0] = Scalar(val0);
      m_storage.data()[1] = Scalar(val1);
    }

    template<typename T0, typename T1>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init2(const Index& val0, const Index& val1,
                                    typename internal::enable_if<    (!internal::is_same<Index,Scalar>::value)
                                                                  && (internal::is_same<T0,Index>::value)
                                                                  && (internal::is_same<T1,Index>::value)
                                                                  && Base::SizeAtCompileTime==2,T1>::type* = 0)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(PlainObjectBase, 2)
      m_storage.data()[0] = Scalar(val0);
      m_storage.data()[1] = Scalar(val1);
    }

    template<typename T>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init1(Index size, typename internal::enable_if<    (Base::SizeAtCompileTime!=1 || !internal::is_convertible<T, Scalar>::value)
                                                                              && ((!internal::is_same<typename internal::traits<Derived>::XprKind,ArrayXpr>::value || Base::SizeAtCompileTime==Dynamic)),T>::type* = 0)
    {
      
      const bool is_integer_alike = internal::is_valid_index_type<T>::value;
      EIGEN_UNUSED_VARIABLE(is_integer_alike);
      EIGEN_STATIC_ASSERT(is_integer_alike,
                          FLOATING_POINT_ARGUMENT_PASSED__INTEGER_WAS_EXPECTED)
      resize(size);
    }

    template<typename T>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init1(const Scalar& val0, typename internal::enable_if<Base::SizeAtCompileTime==1 && internal::is_convertible<T, Scalar>::value,T>::type* = 0)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(PlainObjectBase, 1)
      m_storage.data()[0] = val0;
    }

    template<typename T>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init1(const Index& val0,
                                    typename internal::enable_if<    (!internal::is_same<Index,Scalar>::value)
                                                                  && (internal::is_same<Index,T>::value)
                                                                  && Base::SizeAtCompileTime==1
                                                                  && internal::is_convertible<T, Scalar>::value,T*>::type* = 0)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(PlainObjectBase, 1)
      m_storage.data()[0] = Scalar(val0);
    }

    template<typename T>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init1(const Scalar* data){
      this->_set_noalias(ConstMapType(data));
    }

    template<typename T, typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init1(const DenseBase<OtherDerived>& other){
      this->_set_noalias(other);
    }

    template<typename T>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init1(const Derived& other){
      this->_set_noalias(other);
    }

    template<typename T, typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init1(const EigenBase<OtherDerived>& other){
      this->derived() = other;
    }

    template<typename T, typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init1(const ReturnByValue<OtherDerived>& other)
    {
      resize(other.rows(), other.cols());
      other.evalTo(this->derived());
    }

    template<typename T, typename OtherDerived, int ColsAtCompileTime>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init1(const RotationBase<OtherDerived,ColsAtCompileTime>& r)
    {
      this->derived() = r;
    }

    template<typename T>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init1(const Scalar& val0,
                                    typename internal::enable_if<    Base::SizeAtCompileTime!=Dynamic
                                                                  && Base::SizeAtCompileTime!=1
                                                                  && internal::is_convertible<T, Scalar>::value
                                                                  && internal::is_same<typename internal::traits<Derived>::XprKind,ArrayXpr>::value,T>::type* = 0)
    {
      Base::setConstant(val0);
    }

    template<typename T>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE void _init1(const Index& val0,
                                    typename internal::enable_if<    (!internal::is_same<Index,Scalar>::value)
                                                                  && (internal::is_same<Index,T>::value)
                                                                  && Base::SizeAtCompileTime!=Dynamic
                                                                  && Base::SizeAtCompileTime!=1
                                                                  && internal::is_convertible<T, Scalar>::value
                                                                  && internal::is_same<typename internal::traits<Derived>::XprKind,ArrayXpr>::value,T*>::type* = 0)
    {
      Base::setConstant(val0);
    }

    template<typename MatrixTypeA, typename MatrixTypeB, bool SwapPointers>
    friend struct internal::matrix_swap_impl;

  public:

#ifndef EIGEN_PARSED_BY_DOXYGEN

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    void swap(DenseBase<OtherDerived> & other)
    {
      enum { SwapPointers = internal::is_same<Derived, OtherDerived>::value && Base::SizeAtCompileTime==Dynamic };
      internal::matrix_swap_impl<Derived, OtherDerived, bool(SwapPointers)>::run(this->derived(), other.derived());
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    void swap(DenseBase<OtherDerived> const & other)
    { Base::swap(other.derived()); }

    EIGEN_DEVICE_FUNC
    static EIGEN_STRONG_INLINE void _check_template_params()
    {
      EIGEN_STATIC_ASSERT((EIGEN_IMPLIES(MaxRowsAtCompileTime==1 && MaxColsAtCompileTime!=1, (int(Options)&RowMajor)==RowMajor)
                        && EIGEN_IMPLIES(MaxColsAtCompileTime==1 && MaxRowsAtCompileTime!=1, (int(Options)&RowMajor)==0)
                        && ((RowsAtCompileTime == Dynamic) || (RowsAtCompileTime >= 0))
                        && ((ColsAtCompileTime == Dynamic) || (ColsAtCompileTime >= 0))
                        && ((MaxRowsAtCompileTime == Dynamic) || (MaxRowsAtCompileTime >= 0))
                        && ((MaxColsAtCompileTime == Dynamic) || (MaxColsAtCompileTime >= 0))
                        && (MaxRowsAtCompileTime == RowsAtCompileTime || RowsAtCompileTime==Dynamic)
                        && (MaxColsAtCompileTime == ColsAtCompileTime || ColsAtCompileTime==Dynamic)
                        && (Options & (DontAlign|RowMajor)) == Options),
        INVALID_MATRIX_TEMPLATE_PARAMETERS)
    }

    enum { IsPlainObjectBase = 1 };
#endif
  public:

    template<typename PlainObjectType, int MapOptions, typename StrideType> friend class Eigen::Map;
    friend class Eigen::Map<Derived, Unaligned>;
    friend class Eigen::Map<const Derived, Unaligned>;
#if EIGEN_MAX_ALIGN_BYTES>0
    
    friend class Eigen::Map<Derived, AlignedMax>;
    friend class Eigen::Map<const Derived, AlignedMax>;
#endif
};

namespace internal {

template <typename Derived, typename OtherDerived, bool IsVector>
struct conservative_resize_like_impl
{
  #if EIGEN_HAS_TYPE_TRAITS
  static const bool IsRelocatable = std::is_trivially_copyable<typename Derived::Scalar>::value;
  #else
  static const bool IsRelocatable = !NumTraits<typename Derived::Scalar>::RequireInitialization;
  #endif
  static void run(DenseBase<Derived>& _this, Index rows, Index cols)
  {
    if (_this.rows() == rows && _this.cols() == cols) return;
    EIGEN_STATIC_ASSERT_DYNAMIC_SIZE(Derived)

    if ( IsRelocatable
          && (( Derived::IsRowMajor && _this.cols() == cols) ||  
              (!Derived::IsRowMajor && _this.rows() == rows) ))  
    {
      internal::check_rows_cols_for_overflow<Derived::MaxSizeAtCompileTime>::run(rows, cols);
      _this.derived().m_storage.conservativeResize(rows*cols,rows,cols);
    }
    else
    {
      
      Derived tmp(rows,cols);
      const Index common_rows = numext::mini(rows, _this.rows());
      const Index common_cols = numext::mini(cols, _this.cols());
      tmp.block(0,0,common_rows,common_cols) = _this.block(0,0,common_rows,common_cols);
      _this.derived().swap(tmp);
    }
  }

  static void run(DenseBase<Derived>& _this, const DenseBase<OtherDerived>& other)
  {
    if (_this.rows() == other.rows() && _this.cols() == other.cols()) return;

    EIGEN_STATIC_ASSERT_DYNAMIC_SIZE(Derived)
    EIGEN_STATIC_ASSERT_DYNAMIC_SIZE(OtherDerived)

    if ( IsRelocatable &&
          (( Derived::IsRowMajor && _this.cols() == other.cols()) ||  
           (!Derived::IsRowMajor && _this.rows() == other.rows()) ))  
    {
      const Index new_rows = other.rows() - _this.rows();
      const Index new_cols = other.cols() - _this.cols();
      _this.derived().m_storage.conservativeResize(other.size(),other.rows(),other.cols());
      if (new_rows>0)
        _this.bottomRightCorner(new_rows, other.cols()) = other.bottomRows(new_rows);
      else if (new_cols>0)
        _this.bottomRightCorner(other.rows(), new_cols) = other.rightCols(new_cols);
    }
    else
    {
      
      Derived tmp(other);
      const Index common_rows = numext::mini(tmp.rows(), _this.rows());
      const Index common_cols = numext::mini(tmp.cols(), _this.cols());
      tmp.block(0,0,common_rows,common_cols) = _this.block(0,0,common_rows,common_cols);
      _this.derived().swap(tmp);
    }
  }
};

template <typename Derived, typename OtherDerived>
struct conservative_resize_like_impl<Derived,OtherDerived,true>
  : conservative_resize_like_impl<Derived,OtherDerived,false>
{
  typedef conservative_resize_like_impl<Derived,OtherDerived,false> Base;
  using Base::run;
  using Base::IsRelocatable;

  static void run(DenseBase<Derived>& _this, Index size)
  {
    const Index new_rows = Derived::RowsAtCompileTime==1 ? 1 : size;
    const Index new_cols = Derived::RowsAtCompileTime==1 ? size : 1;
    if(IsRelocatable)
      _this.derived().m_storage.conservativeResize(size,new_rows,new_cols);
    else
      Base::run(_this.derived(), new_rows, new_cols);
  }

  static void run(DenseBase<Derived>& _this, const DenseBase<OtherDerived>& other)
  {
    if (_this.rows() == other.rows() && _this.cols() == other.cols()) return;

    const Index num_new_elements = other.size() - _this.size();

    const Index new_rows = Derived::RowsAtCompileTime==1 ? 1 : other.rows();
    const Index new_cols = Derived::RowsAtCompileTime==1 ? other.cols() : 1;
    if(IsRelocatable)
      _this.derived().m_storage.conservativeResize(other.size(),new_rows,new_cols);
    else
      Base::run(_this.derived(), new_rows, new_cols);

    if (num_new_elements > 0)
      _this.tail(num_new_elements) = other.tail(num_new_elements);
  }
};

template<typename MatrixTypeA, typename MatrixTypeB, bool SwapPointers>
struct matrix_swap_impl
{
  EIGEN_DEVICE_FUNC
  static EIGEN_STRONG_INLINE void run(MatrixTypeA& a, MatrixTypeB& b)
  {
    a.base().swap(b);
  }
};

template<typename MatrixTypeA, typename MatrixTypeB>
struct matrix_swap_impl<MatrixTypeA, MatrixTypeB, true>
{
  EIGEN_DEVICE_FUNC
  static inline void run(MatrixTypeA& a, MatrixTypeB& b)
  {
    static_cast<typename MatrixTypeA::Base&>(a).m_storage.swap(static_cast<typename MatrixTypeB::Base&>(b).m_storage);
  }
};

} 

} 

#endif 
