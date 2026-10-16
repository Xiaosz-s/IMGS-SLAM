#ifndef EIGEN_ARRAY_H
#define EIGEN_ARRAY_H

namespace Eigen {

namespace internal {
template<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
struct traits<Array<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_> > : traits<Matrix<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_> >
{
  typedef ArrayXpr XprKind;
  typedef ArrayBase<Array<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_> > XprBase;
};
}

template<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
class Array
  : public PlainObjectBase<Array<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_> >
{
  public:

    typedef PlainObjectBase<Array> Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(Array)

    enum { Options = Options_ };
    typedef typename Base::PlainObject PlainObject;

  protected:
    template <typename Derived, typename OtherDerived, bool IsVector>
    friend struct internal::conservative_resize_like_impl;

    using Base::m_storage;

  public:

    using Base::base;
    using Base::coeff;
    using Base::coeffRef;

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Array& operator=(const EigenBase<OtherDerived> &other)
    {
      return Base::operator=(other);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Array& operator=(const Scalar &value)
    {
      Base::setConstant(value);
      return *this;
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Array& operator=(const DenseBase<OtherDerived>& other)
    {
      return Base::_set(other);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Array& operator=(const Array& other)
    {
      return Base::_set(other);
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Array() : Base()
    {
      Base::_check_template_params();
      EIGEN_INITIALIZE_COEFFS_IF_THAT_OPTION_IS_ENABLED
    }

#ifndef EIGEN_PARSED_BY_DOXYGEN

    EIGEN_DEVICE_FUNC
    Array(internal::constructor_without_unaligned_array_assert)
      : Base(internal::constructor_without_unaligned_array_assert())
    {
      Base::_check_template_params();
      EIGEN_INITIALIZE_COEFFS_IF_THAT_OPTION_IS_ENABLED
    }
#endif

#if EIGEN_HAS_RVALUE_REFERENCES
    EIGEN_DEVICE_FUNC
    Array(Array&& other) EIGEN_NOEXCEPT_IF(std::is_nothrow_move_constructible<Scalar>::value)
      : Base(std::move(other))
    {
      Base::_check_template_params();
    }
    EIGEN_DEVICE_FUNC
    Array& operator=(Array&& other) EIGEN_NOEXCEPT_IF(std::is_nothrow_move_assignable<Scalar>::value)
    {
      Base::operator=(std::move(other));
      return *this;
    }
#endif

    #if EIGEN_HAS_CXX11

    template <typename... ArgTypes>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Array(const Scalar& a0, const Scalar& a1, const Scalar& a2, const Scalar& a3, const ArgTypes&... args)
      : Base(a0, a1, a2, a3, args...) {}

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Array(const std::initializer_list<std::initializer_list<Scalar>>& list) : Base(list) {}
    #endif 

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    template<typename T>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE explicit Array(const T& x)
    {
      Base::_check_template_params();
      Base::template _init1<T>(x);
    }

    template<typename T0, typename T1>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Array(const T0& val0, const T1& val1)
    {
      Base::_check_template_params();
      this->template _init2<T0,T1>(val0, val1);
    }

    #else
    
    EIGEN_DEVICE_FUNC explicit Array(const Scalar *data);

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE explicit Array(Index dim);

    Array(const Scalar& value);

    Array(Index rows, Index cols);

    Array(const Scalar& val0, const Scalar& val1);
    #endif  

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Array(const Scalar& val0, const Scalar& val1, const Scalar& val2)
    {
      Base::_check_template_params();
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Array, 3)
      m_storage.data()[0] = val0;
      m_storage.data()[1] = val1;
      m_storage.data()[2] = val2;
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Array(const Scalar& val0, const Scalar& val1, const Scalar& val2, const Scalar& val3)
    {
      Base::_check_template_params();
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Array, 4)
      m_storage.data()[0] = val0;
      m_storage.data()[1] = val1;
      m_storage.data()[2] = val2;
      m_storage.data()[3] = val3;
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Array(const Array& other)
            : Base(other)
    { }

  private:
    struct PrivateType {};
  public:

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Array(const EigenBase<OtherDerived> &other,
                              typename internal::enable_if<internal::is_convertible<typename OtherDerived::Scalar,Scalar>::value,
                                                           PrivateType>::type = PrivateType())
      : Base(other.derived())
    { }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index innerStride() const EIGEN_NOEXCEPT{ return 1; }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index outerStride() const EIGEN_NOEXCEPT { return this->innerSize(); }

    #ifdef EIGEN_ARRAY_PLUGIN
    #include EIGEN_ARRAY_PLUGIN
    #endif

  private:

    template<typename MatrixType, typename OtherDerived, bool SwapPointers>
    friend struct internal::matrix_swap_impl;
};

#define EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, Size, SizeSuffix)   \
                                    \
typedef Array<Type, Size, Size> Array##SizeSuffix##SizeSuffix##TypeSuffix;  \
                                    \
typedef Array<Type, Size, 1>    Array##SizeSuffix##TypeSuffix;

#define EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, Size)         \
                                    \
typedef Array<Type, Size, Dynamic> Array##Size##X##TypeSuffix;  \
                                    \
typedef Array<Type, Dynamic, Size> Array##X##Size##TypeSuffix;

#define EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES(Type, TypeSuffix) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 2, 2) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 3, 3) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 4, 4) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, Dynamic, X) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 2) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 3) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 4)

EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES(int,                  i)
EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES(float,                f)
EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES(double,               d)
EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES(std::complex<float>,  cf)
EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES(std::complex<double>, cd)

#undef EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES
#undef EIGEN_MAKE_ARRAY_TYPEDEFS
#undef EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS

#if EIGEN_HAS_CXX11

#define EIGEN_MAKE_ARRAY_TYPEDEFS(Size, SizeSuffix)               \
                                     \
                                              \
template <typename Type>                                          \
using Array##SizeSuffix##SizeSuffix = Array<Type, Size, Size>;    \
                                     \
                                              \
template <typename Type>                                          \
using Array##SizeSuffix = Array<Type, Size, 1>;

#define EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Size)                     \
                                     \
                                              \
template <typename Type>                                          \
using Array##Size##X = Array<Type, Size, Dynamic>;                \
                                     \
                                              \
template <typename Type>                                          \
using Array##X##Size = Array<Type, Dynamic, Size>;

EIGEN_MAKE_ARRAY_TYPEDEFS(2, 2)
EIGEN_MAKE_ARRAY_TYPEDEFS(3, 3)
EIGEN_MAKE_ARRAY_TYPEDEFS(4, 4)
EIGEN_MAKE_ARRAY_TYPEDEFS(Dynamic, X)
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(2)
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(3)
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(4)

#undef EIGEN_MAKE_ARRAY_TYPEDEFS
#undef EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS

#endif 

#define EIGEN_USING_ARRAY_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, SizeSuffix) \
using Eigen::Matrix##SizeSuffix##TypeSuffix; \
using Eigen::Vector##SizeSuffix##TypeSuffix; \
using Eigen::RowVector##SizeSuffix##TypeSuffix;

#define EIGEN_USING_ARRAY_TYPEDEFS_FOR_TYPE(TypeSuffix) \
EIGEN_USING_ARRAY_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, 2) \
EIGEN_USING_ARRAY_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, 3) \
EIGEN_USING_ARRAY_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, 4) \
EIGEN_USING_ARRAY_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, X) \

#define EIGEN_USING_ARRAY_TYPEDEFS \
EIGEN_USING_ARRAY_TYPEDEFS_FOR_TYPE(i) \
EIGEN_USING_ARRAY_TYPEDEFS_FOR_TYPE(f) \
EIGEN_USING_ARRAY_TYPEDEFS_FOR_TYPE(d) \
EIGEN_USING_ARRAY_TYPEDEFS_FOR_TYPE(cf) \
EIGEN_USING_ARRAY_TYPEDEFS_FOR_TYPE(cd)

} 

#endif 
