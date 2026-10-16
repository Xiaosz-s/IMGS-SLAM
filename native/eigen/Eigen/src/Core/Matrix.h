#ifndef EIGEN_MATRIX_H
#define EIGEN_MATRIX_H

namespace Eigen {

namespace internal {
template<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
struct traits<Matrix<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_> >
{
private:
  enum { size = internal::size_at_compile_time<Rows_,Cols_>::ret };
  typedef typename find_best_packet<Scalar_,size>::type PacketScalar;
  enum {
      row_major_bit = Options_&RowMajor ? RowMajorBit : 0,
      is_dynamic_size_storage = MaxRows_==Dynamic || MaxCols_==Dynamic,
      max_size = is_dynamic_size_storage ? Dynamic : MaxRows_*MaxCols_,
      default_alignment = compute_default_alignment<Scalar_,max_size>::value,
      actual_alignment = ((Options_&DontAlign)==0) ? default_alignment : 0,
      required_alignment = unpacket_traits<PacketScalar>::alignment,
      packet_access_bit = (packet_traits<Scalar_>::Vectorizable && (EIGEN_UNALIGNED_VECTORIZE || (actual_alignment>=required_alignment))) ? PacketAccessBit : 0
    };

public:
  typedef Scalar_ Scalar;
  typedef Dense StorageKind;
  typedef Eigen::Index StorageIndex;
  typedef MatrixXpr XprKind;
  enum {
    RowsAtCompileTime = Rows_,
    ColsAtCompileTime = Cols_,
    MaxRowsAtCompileTime = MaxRows_,
    MaxColsAtCompileTime = MaxCols_,
    Flags = compute_matrix_flags<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>::ret,
    Options = Options_,
    InnerStrideAtCompileTime = 1,
    OuterStrideAtCompileTime = (Options&RowMajor) ? ColsAtCompileTime : RowsAtCompileTime,

    EvaluatorFlags = LinearAccessBit | DirectAccessBit | packet_access_bit | row_major_bit,
    Alignment = actual_alignment
  };
};
}

template<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
class Matrix
  : public PlainObjectBase<Matrix<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_> >
{
  public:

    typedef PlainObjectBase<Matrix> Base;

    enum { Options = Options_ };

    EIGEN_DENSE_PUBLIC_INTERFACE(Matrix)

    typedef typename Base::PlainObject PlainObject;

    using Base::base;
    using Base::coeffRef;

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Matrix& operator=(const Matrix& other)
    {
      return Base::_set(other);
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Matrix& operator=(const DenseBase<OtherDerived>& other)
    {
      return Base::_set(other);
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Matrix& operator=(const EigenBase<OtherDerived> &other)
    {
      return Base::operator=(other);
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Matrix& operator=(const ReturnByValue<OtherDerived>& func)
    {
      return Base::operator=(func);
    }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Matrix() : Base()
    {
      Base::_check_template_params();
      EIGEN_INITIALIZE_COEFFS_IF_THAT_OPTION_IS_ENABLED
    }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    explicit Matrix(internal::constructor_without_unaligned_array_assert)
      : Base(internal::constructor_without_unaligned_array_assert())
    { Base::_check_template_params(); EIGEN_INITIALIZE_COEFFS_IF_THAT_OPTION_IS_ENABLED }

#if EIGEN_HAS_RVALUE_REFERENCES
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Matrix(Matrix&& other) EIGEN_NOEXCEPT_IF(std::is_nothrow_move_constructible<Scalar>::value)
      : Base(std::move(other))
    {
      Base::_check_template_params();
    }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Matrix& operator=(Matrix&& other) EIGEN_NOEXCEPT_IF(std::is_nothrow_move_assignable<Scalar>::value)
    {
      Base::operator=(std::move(other));
      return *this;
    }
#endif

#if EIGEN_HAS_CXX11

    template <typename... ArgTypes>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Matrix(const Scalar& a0, const Scalar& a1, const Scalar& a2,  const Scalar& a3, const ArgTypes&... args)
      : Base(a0, a1, a2, a3, args...) {}

    EIGEN_DEVICE_FUNC
    explicit EIGEN_STRONG_INLINE Matrix(const std::initializer_list<std::initializer_list<Scalar>>& list) : Base(list) {}
#endif 

#ifndef EIGEN_PARSED_BY_DOXYGEN

    template<typename T>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    explicit Matrix(const T& x)
    {
      Base::_check_template_params();
      Base::template _init1<T>(x);
    }

    template<typename T0, typename T1>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Matrix(const T0& x, const T1& y)
    {
      Base::_check_template_params();
      Base::template _init2<T0,T1>(x, y);
    }

#else
    
    EIGEN_DEVICE_FUNC
    explicit Matrix(const Scalar *data);

    EIGEN_STRONG_INLINE explicit Matrix(Index dim);

    Matrix(const Scalar& x);

    EIGEN_DEVICE_FUNC
    Matrix(Index rows, Index cols);

    Matrix(const Scalar& x, const Scalar& y);
    #endif  

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Matrix(const Scalar& x, const Scalar& y, const Scalar& z)
    {
      Base::_check_template_params();
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Matrix, 3)
      m_storage.data()[0] = x;
      m_storage.data()[1] = y;
      m_storage.data()[2] = z;
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Matrix(const Scalar& x, const Scalar& y, const Scalar& z, const Scalar& w)
    {
      Base::_check_template_params();
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Matrix, 4)
      m_storage.data()[0] = x;
      m_storage.data()[1] = y;
      m_storage.data()[2] = z;
      m_storage.data()[3] = w;
    }

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Matrix(const Matrix& other) : Base(other)
    { }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Matrix(const EigenBase<OtherDerived> &other)
      : Base(other.derived())
    { }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index innerStride() const EIGEN_NOEXCEPT { return 1; }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index outerStride() const EIGEN_NOEXCEPT { return this->innerSize(); }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    explicit Matrix(const RotationBase<OtherDerived,ColsAtCompileTime>& r);
    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    Matrix& operator=(const RotationBase<OtherDerived,ColsAtCompileTime>& r);

    #ifdef EIGEN_MATRIX_PLUGIN
    #include EIGEN_MATRIX_PLUGIN
    #endif

  protected:
    template <typename Derived, typename OtherDerived, bool IsVector>
    friend struct internal::conservative_resize_like_impl;

    using Base::m_storage;
};

#define EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, Size, SizeSuffix)   \
                                    \
typedef Matrix<Type, Size, Size> Matrix##SizeSuffix##TypeSuffix;  \
                                    \
typedef Matrix<Type, Size, 1>    Vector##SizeSuffix##TypeSuffix;  \
                                    \
typedef Matrix<Type, 1, Size>    RowVector##SizeSuffix##TypeSuffix;

#define EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, Size)         \
                                    \
typedef Matrix<Type, Size, Dynamic> Matrix##Size##X##TypeSuffix;  \
                                    \
typedef Matrix<Type, Dynamic, Size> Matrix##X##Size##TypeSuffix;

#define EIGEN_MAKE_TYPEDEFS_ALL_SIZES(Type, TypeSuffix) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 2, 2) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 3, 3) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 4, 4) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, Dynamic, X) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 2) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 3) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 4)

EIGEN_MAKE_TYPEDEFS_ALL_SIZES(int,                  i)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(float,                f)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(double,               d)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(std::complex<float>,  cf)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(std::complex<double>, cd)

#undef EIGEN_MAKE_TYPEDEFS_ALL_SIZES
#undef EIGEN_MAKE_TYPEDEFS
#undef EIGEN_MAKE_FIXED_TYPEDEFS

#if EIGEN_HAS_CXX11

#define EIGEN_MAKE_TYPEDEFS(Size, SizeSuffix)                     \
                                    \
                                              \
template <typename Type>                                          \
using Matrix##SizeSuffix = Matrix<Type, Size, Size>;              \
                                    \
                                              \
template <typename Type>                                          \
using Vector##SizeSuffix = Matrix<Type, Size, 1>;                 \
                                    \
                                              \
template <typename Type>                                          \
using RowVector##SizeSuffix = Matrix<Type, 1, Size>;

#define EIGEN_MAKE_FIXED_TYPEDEFS(Size)                           \
                                    \
                                              \
template <typename Type>                                          \
using Matrix##Size##X = Matrix<Type, Size, Dynamic>;              \
                                    \
                                              \
template <typename Type>                                          \
using Matrix##X##Size = Matrix<Type, Dynamic, Size>;

EIGEN_MAKE_TYPEDEFS(2, 2)
EIGEN_MAKE_TYPEDEFS(3, 3)
EIGEN_MAKE_TYPEDEFS(4, 4)
EIGEN_MAKE_TYPEDEFS(Dynamic, X)
EIGEN_MAKE_FIXED_TYPEDEFS(2)
EIGEN_MAKE_FIXED_TYPEDEFS(3)
EIGEN_MAKE_FIXED_TYPEDEFS(4)

template <typename Type, int Size>
using Vector = Matrix<Type, Size, 1>;

template <typename Type, int Size>
using RowVector = Matrix<Type, 1, Size>;

#undef EIGEN_MAKE_TYPEDEFS
#undef EIGEN_MAKE_FIXED_TYPEDEFS

#endif 

} 

#endif 
