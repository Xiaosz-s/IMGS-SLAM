#ifndef EIGEN_XPRHELPER_H
#define EIGEN_XPRHELPER_H

#if EIGEN_COMP_GNUC && !EIGEN_GNUC_AT(4,3)
  #define EIGEN_EMPTY_STRUCT_CTOR(X) \
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE X() {} \
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE X(const X& ) {}
#else
  #define EIGEN_EMPTY_STRUCT_CTOR(X)
#endif

namespace Eigen {

namespace internal {

template<typename IndexDest, typename IndexSrc>
EIGEN_DEVICE_FUNC
inline IndexDest convert_index(const IndexSrc& idx) {
  
  eigen_internal_assert(idx <= NumTraits<IndexDest>::highest() && "Index value to big for target type");
  return IndexDest(idx);
}

template<typename T> struct is_valid_index_type
{
  enum { value =
#if EIGEN_HAS_TYPE_TRAITS
    internal::is_integral<T>::value || std::is_enum<T>::value
#elif EIGEN_COMP_MSVC
    internal::is_integral<T>::value || __is_enum(T)
#else
    
    internal::is_convertible<T,Index>::value && !internal::is_same<T,float>::value && !is_same<T,double>::value
#endif
  };
};

template<typename RowIndices, typename ColIndices>
struct valid_indexed_view_overload {
  enum { value = !(internal::is_valid_index_type<RowIndices>::value && internal::is_valid_index_type<ColIndices>::value) };
};

template<typename ExprScalar,typename T, bool IsSupported>
struct promote_scalar_arg;

template<typename S,typename T>
struct promote_scalar_arg<S,T,true>
{
  typedef T type;
};

template<typename ExprScalar,typename T,typename PromotedType,
  bool ConvertibleToLiteral = internal::is_convertible<T,PromotedType>::value,
  bool IsSafe = NumTraits<T>::IsInteger || !NumTraits<PromotedType>::IsInteger>
struct promote_scalar_arg_unsupported;

template<typename S,typename T>
struct promote_scalar_arg<S,T,false> : promote_scalar_arg_unsupported<S,T,typename NumTraits<S>::Literal> {};

template<typename S,typename T, typename PromotedType>
struct promote_scalar_arg_unsupported<S,T,PromotedType,true,true>
{
  typedef PromotedType type;
};

template<typename ExprScalar,typename T, typename PromotedType>
struct promote_scalar_arg_unsupported<ExprScalar,T,PromotedType,false,true>
   : promote_scalar_arg_unsupported<ExprScalar,T,ExprScalar>
{};

template<typename S,typename T, typename PromotedType, bool ConvertibleToLiteral>
struct promote_scalar_arg_unsupported<S,T,PromotedType,ConvertibleToLiteral,false> {};

template<typename S,typename T>
struct promote_scalar_arg_unsupported<S,T,S,false,true> {};

class no_assignment_operator
{
  private:
    no_assignment_operator& operator=(const no_assignment_operator&);
  protected:
    EIGEN_DEFAULT_COPY_CONSTRUCTOR(no_assignment_operator)
    EIGEN_DEFAULT_EMPTY_CONSTRUCTOR_AND_DESTRUCTOR(no_assignment_operator)
};

template<typename I1, typename I2>
struct promote_index_type
{
  typedef typename conditional<(sizeof(I1)<sizeof(I2)), I2, I1>::type type;
};

template<typename T, int Value> class variable_if_dynamic
{
  public:
    EIGEN_DEFAULT_EMPTY_CONSTRUCTOR_AND_DESTRUCTOR(variable_if_dynamic)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE explicit variable_if_dynamic(T v) { EIGEN_ONLY_USED_FOR_DEBUG(v); eigen_assert(v == T(Value)); }
    EIGEN_DEVICE_FUNC static EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    T value() { return T(Value); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    operator T() const { return T(Value); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    void setValue(T v) const { EIGEN_ONLY_USED_FOR_DEBUG(v); eigen_assert(v == T(Value)); }
};

template<typename T> class variable_if_dynamic<T, Dynamic>
{
    T m_value;
  public:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE explicit variable_if_dynamic(T value = 0) EIGEN_NO_THROW : m_value(value) {}
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T value() const { return m_value; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE operator T() const { return m_value; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void setValue(T value) { m_value = value; }
};

template<typename T, int Value> class variable_if_dynamicindex
{
  public:
    EIGEN_EMPTY_STRUCT_CTOR(variable_if_dynamicindex)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE explicit variable_if_dynamicindex(T v) { EIGEN_ONLY_USED_FOR_DEBUG(v); eigen_assert(v == T(Value)); }
    EIGEN_DEVICE_FUNC static EIGEN_STRONG_INLINE EIGEN_CONSTEXPR
    T value() { return T(Value); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    void setValue(T) {}
};

template<typename T> class variable_if_dynamicindex<T, DynamicIndex>
{
    T m_value;
    EIGEN_DEVICE_FUNC variable_if_dynamicindex() { eigen_assert(false); }
  public:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE explicit variable_if_dynamicindex(T value) : m_value(value) {}
    EIGEN_DEVICE_FUNC T EIGEN_STRONG_INLINE value() const { return m_value; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void setValue(T value) { m_value = value; }
};

template<typename T> struct functor_traits
{
  enum
  {
    Cost = 10,
    PacketAccess = false,
    IsRepeatable = false
  };
};

template<typename T> struct packet_traits;

template<typename T> struct unpacket_traits;

template<int Size, typename PacketType,
         bool Stop = Size==Dynamic || (Size%unpacket_traits<PacketType>::size)==0 || is_same<PacketType,typename unpacket_traits<PacketType>::half>::value>
struct find_best_packet_helper;

template< int Size, typename PacketType>
struct find_best_packet_helper<Size,PacketType,true>
{
  typedef PacketType type;
};

template<int Size, typename PacketType>
struct find_best_packet_helper<Size,PacketType,false>
{
  typedef typename find_best_packet_helper<Size,typename unpacket_traits<PacketType>::half>::type type;
};

template<typename T, int Size>
struct find_best_packet
{
  typedef typename find_best_packet_helper<Size,typename packet_traits<T>::type>::type type;
};

#if EIGEN_MAX_STATIC_ALIGN_BYTES>0
template<int ArrayBytes, int AlignmentBytes,
         bool Match     =  bool((ArrayBytes%AlignmentBytes)==0),
         bool TryHalf   =  bool(EIGEN_MIN_ALIGN_BYTES<AlignmentBytes) >
struct compute_default_alignment_helper
{
  enum { value = 0 };
};

template<int ArrayBytes, int AlignmentBytes, bool TryHalf>
struct compute_default_alignment_helper<ArrayBytes, AlignmentBytes, true, TryHalf> 
{
  enum { value = AlignmentBytes };
};

template<int ArrayBytes, int AlignmentBytes>
struct compute_default_alignment_helper<ArrayBytes, AlignmentBytes, false, true> 
{
  
  enum { value = compute_default_alignment_helper<ArrayBytes, AlignmentBytes/2>::value };
};
#else

template<int ArrayBytes, int AlignmentBytes>
struct compute_default_alignment_helper
{
  enum { value = 0 };
};
#endif

template<typename T, int Size> struct compute_default_alignment {
  enum { value = compute_default_alignment_helper<Size*sizeof(T),EIGEN_MAX_STATIC_ALIGN_BYTES>::value };
};

template<typename T> struct compute_default_alignment<T,Dynamic> {
  enum { value = EIGEN_MAX_ALIGN_BYTES };
};

template<typename Scalar_, int Rows_, int Cols_,
         int Options_ = AutoAlign |
                          ( (Rows_==1 && Cols_!=1) ? RowMajor
                          : (Cols_==1 && Rows_!=1) ? ColMajor
                          : EIGEN_DEFAULT_MATRIX_STORAGE_ORDER_OPTION ),
         int MaxRows_ = Rows_,
         int MaxCols_ = Cols_
> class make_proper_matrix_type
{
    enum {
      IsColVector = Cols_==1 && Rows_!=1,
      IsRowVector = Rows_==1 && Cols_!=1,
      Options = IsColVector ? (Options_ | ColMajor) & ~RowMajor
              : IsRowVector ? (Options_ | RowMajor) & ~ColMajor
              : Options_
    };
  public:
    typedef Matrix<Scalar_, Rows_, Cols_, Options, MaxRows_, MaxCols_> type;
};

template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
class compute_matrix_flags
{
    enum { row_major_bit = Options&RowMajor ? RowMajorBit : 0 };
  public:

    enum { ret = DirectAccessBit | LvalueBit | NestByRefBit | row_major_bit };
};

template<int Rows_, int Cols_> struct size_at_compile_time
{
  enum { ret = (Rows_==Dynamic || Cols_==Dynamic) ? Dynamic : Rows_ * Cols_ };
};

template<typename XprType> struct size_of_xpr_at_compile_time
{
  enum { ret = size_at_compile_time<traits<XprType>::RowsAtCompileTime,traits<XprType>::ColsAtCompileTime>::ret };
};

template<typename T, typename StorageKind = typename traits<T>::StorageKind> struct plain_matrix_type;
template<typename T, typename BaseClassType, int Flags> struct plain_matrix_type_dense;
template<typename T> struct plain_matrix_type<T,Dense>
{
  typedef typename plain_matrix_type_dense<T,typename traits<T>::XprKind, traits<T>::Flags>::type type;
};
template<typename T> struct plain_matrix_type<T,DiagonalShape>
{
  typedef typename T::PlainObject type;
};

template<typename T, int Flags> struct plain_matrix_type_dense<T,MatrixXpr,Flags>
{
  typedef Matrix<typename traits<T>::Scalar,
                traits<T>::RowsAtCompileTime,
                traits<T>::ColsAtCompileTime,
                AutoAlign | (Flags&RowMajorBit ? RowMajor : ColMajor),
                traits<T>::MaxRowsAtCompileTime,
                traits<T>::MaxColsAtCompileTime
          > type;
};

template<typename T, int Flags> struct plain_matrix_type_dense<T,ArrayXpr,Flags>
{
  typedef Array<typename traits<T>::Scalar,
                traits<T>::RowsAtCompileTime,
                traits<T>::ColsAtCompileTime,
                AutoAlign | (Flags&RowMajorBit ? RowMajor : ColMajor),
                traits<T>::MaxRowsAtCompileTime,
                traits<T>::MaxColsAtCompileTime
          > type;
};

template<typename T, typename StorageKind = typename traits<T>::StorageKind> struct eval;

template<typename T> struct eval<T,Dense>
{
  typedef typename plain_matrix_type<T>::type type;

};

template<typename T> struct eval<T,DiagonalShape>
{
  typedef typename plain_matrix_type<T>::type type;
};

template<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
struct eval<Matrix<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>, Dense>
{
  typedef const Matrix<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>& type;
};

template<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
struct eval<Array<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>, Dense>
{
  typedef const Array<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>& type;
};

template<typename T, typename StorageKind = typename traits<T>::StorageKind> struct plain_object_eval;

template<typename T>
struct plain_object_eval<T,Dense>
{
  typedef typename plain_matrix_type_dense<T,typename traits<T>::XprKind, evaluator<T>::Flags>::type type;
};

template<typename T> struct plain_matrix_type_column_major
{
  enum { Rows = traits<T>::RowsAtCompileTime,
         Cols = traits<T>::ColsAtCompileTime,
         MaxRows = traits<T>::MaxRowsAtCompileTime,
         MaxCols = traits<T>::MaxColsAtCompileTime
  };
  typedef Matrix<typename traits<T>::Scalar,
                Rows,
                Cols,
                (MaxRows==1&&MaxCols!=1) ? RowMajor : ColMajor,
                MaxRows,
                MaxCols
          > type;
};

template<typename T> struct plain_matrix_type_row_major
{
  enum { Rows = traits<T>::RowsAtCompileTime,
         Cols = traits<T>::ColsAtCompileTime,
         MaxRows = traits<T>::MaxRowsAtCompileTime,
         MaxCols = traits<T>::MaxColsAtCompileTime
  };
  typedef Matrix<typename traits<T>::Scalar,
                Rows,
                Cols,
                (MaxCols==1&&MaxRows!=1) ? ColMajor : RowMajor,
                MaxRows,
                MaxCols
          > type;
};

template <typename T>
struct ref_selector
{
  typedef typename conditional<
    bool(traits<T>::Flags & NestByRefBit),
    T const&,
    const T
  >::type type;

  typedef typename conditional<
    bool(traits<T>::Flags & NestByRefBit),
    T &,
    T
  >::type non_const_type;
};

template<typename T1, typename T2>
struct transfer_constness
{
  typedef typename conditional<
    bool(internal::is_const<T1>::value),
    typename internal::add_const_on_value_type<T2>::type,
    T2
  >::type type;
};

template<typename T, int n, typename PlainObject = typename plain_object_eval<T>::type> struct nested_eval
{
  enum {
    ScalarReadCost = NumTraits<typename traits<T>::Scalar>::ReadCost,
    CoeffReadCost = evaluator<T>::CoeffReadCost,  

    NAsInteger = n == Dynamic ? HugeCost : n,
    CostEval   = (NAsInteger+1) * ScalarReadCost + CoeffReadCost,
    CostNoEval = NAsInteger * CoeffReadCost,
    Evaluate = (int(evaluator<T>::Flags) & EvalBeforeNestingBit) || (int(CostEval) < int(CostNoEval))
  };

  typedef typename conditional<Evaluate, PlainObject, typename ref_selector<T>::type>::type type;
};

template<typename T>
EIGEN_DEVICE_FUNC
inline T* const_cast_ptr(const T* ptr)
{
  return const_cast<T*>(ptr);
}

template<typename Derived, typename XprKind = typename traits<Derived>::XprKind>
struct dense_xpr_base
{
  
};

template<typename Derived>
struct dense_xpr_base<Derived, MatrixXpr>
{
  typedef MatrixBase<Derived> type;
};

template<typename Derived>
struct dense_xpr_base<Derived, ArrayXpr>
{
  typedef ArrayBase<Derived> type;
};

template<typename Derived, typename XprKind = typename traits<Derived>::XprKind, typename StorageKind = typename traits<Derived>::StorageKind>
struct generic_xpr_base;

template<typename Derived, typename XprKind>
struct generic_xpr_base<Derived, XprKind, Dense>
{
  typedef typename dense_xpr_base<Derived,XprKind>::type type;
};

template<typename XprType, typename CastType> struct cast_return_type
{
  typedef typename XprType::Scalar CurrentScalarType;
  typedef typename remove_all<CastType>::type _CastType;
  typedef typename _CastType::Scalar NewScalarType;
  typedef typename conditional<is_same<CurrentScalarType,NewScalarType>::value,
                              const XprType&,CastType>::type type;
};

template <typename A, typename B> struct promote_storage_type;

template <typename A> struct promote_storage_type<A,A>
{
  typedef A ret;
};
template <typename A> struct promote_storage_type<A, const A>
{
  typedef A ret;
};
template <typename A> struct promote_storage_type<const A, A>
{
  typedef A ret;
};

template <typename A, typename B, typename Functor> struct cwise_promote_storage_type;

template <typename A, typename Functor>                   struct cwise_promote_storage_type<A,A,Functor>                                      { typedef A      ret; };
template <typename Functor>                               struct cwise_promote_storage_type<Dense,Dense,Functor>                              { typedef Dense  ret; };
template <typename A, typename Functor>                   struct cwise_promote_storage_type<A,Dense,Functor>                                  { typedef Dense  ret; };
template <typename B, typename Functor>                   struct cwise_promote_storage_type<Dense,B,Functor>                                  { typedef Dense  ret; };
template <typename Functor>                               struct cwise_promote_storage_type<Sparse,Dense,Functor>                             { typedef Sparse ret; };
template <typename Functor>                               struct cwise_promote_storage_type<Dense,Sparse,Functor>                             { typedef Sparse ret; };

template <typename LhsKind, typename RhsKind, int LhsOrder, int RhsOrder> struct cwise_promote_storage_order {
  enum { value = LhsOrder };
};

template <typename LhsKind, int LhsOrder, int RhsOrder>   struct cwise_promote_storage_order<LhsKind,Sparse,LhsOrder,RhsOrder>                { enum { value = RhsOrder }; };
template <typename RhsKind, int LhsOrder, int RhsOrder>   struct cwise_promote_storage_order<Sparse,RhsKind,LhsOrder,RhsOrder>                { enum { value = LhsOrder }; };
template <int Order>                                      struct cwise_promote_storage_order<Sparse,Sparse,Order,Order>                       { enum { value = Order }; };

template <typename A, typename B, int ProductTag> struct product_promote_storage_type;

template <typename A, int ProductTag> struct product_promote_storage_type<A,                  A,                  ProductTag> { typedef A     ret;};
template <int ProductTag>             struct product_promote_storage_type<Dense,              Dense,              ProductTag> { typedef Dense ret;};
template <typename A, int ProductTag> struct product_promote_storage_type<A,                  Dense,              ProductTag> { typedef Dense ret; };
template <typename B, int ProductTag> struct product_promote_storage_type<Dense,              B,                  ProductTag> { typedef Dense ret; };

template <typename A, int ProductTag> struct product_promote_storage_type<A,                  DiagonalShape,      ProductTag> { typedef A ret; };
template <typename B, int ProductTag> struct product_promote_storage_type<DiagonalShape,      B,                  ProductTag> { typedef B ret; };
template <int ProductTag>             struct product_promote_storage_type<Dense,              DiagonalShape,      ProductTag> { typedef Dense ret; };
template <int ProductTag>             struct product_promote_storage_type<DiagonalShape,      Dense,              ProductTag> { typedef Dense ret; };

template <typename A, int ProductTag> struct product_promote_storage_type<A,                  PermutationStorage, ProductTag> { typedef A ret; };
template <typename B, int ProductTag> struct product_promote_storage_type<PermutationStorage, B,                  ProductTag> { typedef B ret; };
template <int ProductTag>             struct product_promote_storage_type<Dense,              PermutationStorage, ProductTag> { typedef Dense ret; };
template <int ProductTag>             struct product_promote_storage_type<PermutationStorage, Dense,              ProductTag> { typedef Dense ret; };

template<typename ExpressionType, typename Scalar = typename ExpressionType::Scalar>
struct plain_row_type
{
  typedef Matrix<Scalar, 1, ExpressionType::ColsAtCompileTime,
                 int(ExpressionType::PlainObject::Options) | int(RowMajor), 1, ExpressionType::MaxColsAtCompileTime> MatrixRowType;
  typedef Array<Scalar, 1, ExpressionType::ColsAtCompileTime,
                 int(ExpressionType::PlainObject::Options) | int(RowMajor), 1, ExpressionType::MaxColsAtCompileTime> ArrayRowType;

  typedef typename conditional<
    is_same< typename traits<ExpressionType>::XprKind, MatrixXpr >::value,
    MatrixRowType,
    ArrayRowType
  >::type type;
};

template<typename ExpressionType, typename Scalar = typename ExpressionType::Scalar>
struct plain_col_type
{
  typedef Matrix<Scalar, ExpressionType::RowsAtCompileTime, 1,
                 ExpressionType::PlainObject::Options & ~RowMajor, ExpressionType::MaxRowsAtCompileTime, 1> MatrixColType;
  typedef Array<Scalar, ExpressionType::RowsAtCompileTime, 1,
                 ExpressionType::PlainObject::Options & ~RowMajor, ExpressionType::MaxRowsAtCompileTime, 1> ArrayColType;

  typedef typename conditional<
    is_same< typename traits<ExpressionType>::XprKind, MatrixXpr >::value,
    MatrixColType,
    ArrayColType
  >::type type;
};

template<typename ExpressionType, typename Scalar = typename ExpressionType::Scalar>
struct plain_diag_type
{
  enum { diag_size = EIGEN_SIZE_MIN_PREFER_DYNAMIC(ExpressionType::RowsAtCompileTime, ExpressionType::ColsAtCompileTime),
         max_diag_size = EIGEN_SIZE_MIN_PREFER_FIXED(ExpressionType::MaxRowsAtCompileTime, ExpressionType::MaxColsAtCompileTime)
  };
  typedef Matrix<Scalar, diag_size, 1, ExpressionType::PlainObject::Options & ~RowMajor, max_diag_size, 1> MatrixDiagType;
  typedef Array<Scalar, diag_size, 1, ExpressionType::PlainObject::Options & ~RowMajor, max_diag_size, 1> ArrayDiagType;

  typedef typename conditional<
    is_same< typename traits<ExpressionType>::XprKind, MatrixXpr >::value,
    MatrixDiagType,
    ArrayDiagType
  >::type type;
};

template<typename Expr,typename Scalar = typename Expr::Scalar>
struct plain_constant_type
{
  enum { Options = (traits<Expr>::Flags&RowMajorBit)?RowMajor:0 };

  typedef Array<Scalar,  traits<Expr>::RowsAtCompileTime,   traits<Expr>::ColsAtCompileTime,
                Options, traits<Expr>::MaxRowsAtCompileTime,traits<Expr>::MaxColsAtCompileTime> array_type;

  typedef Matrix<Scalar,  traits<Expr>::RowsAtCompileTime,   traits<Expr>::ColsAtCompileTime,
                 Options, traits<Expr>::MaxRowsAtCompileTime,traits<Expr>::MaxColsAtCompileTime> matrix_type;

  typedef CwiseNullaryOp<scalar_constant_op<Scalar>, const typename conditional<is_same< typename traits<Expr>::XprKind, MatrixXpr >::value, matrix_type, array_type>::type > type;
};

template<typename ExpressionType>
struct is_lvalue
{
  enum { value = (!bool(is_const<ExpressionType>::value)) &&
                 bool(traits<ExpressionType>::Flags & LvalueBit) };
};

template<typename T> struct is_diagonal
{ enum { ret = false }; };

template<typename T> struct is_diagonal<DiagonalBase<T> >
{ enum { ret = true }; };

template<typename T> struct is_diagonal<DiagonalWrapper<T> >
{ enum { ret = true }; };

template<typename T, int S> struct is_diagonal<DiagonalMatrix<T,S> >
{ enum { ret = true }; };

template<typename T> struct is_identity
{ enum { value = false }; };

template<typename T> struct is_identity<CwiseNullaryOp<internal::scalar_identity_op<typename T::Scalar>, T> >
{ enum { value = true }; };

template<typename S1, typename S2> struct glue_shapes;
template<> struct glue_shapes<DenseShape,TriangularShape> { typedef TriangularShape type;  };

template<typename T1, typename T2>
struct possibly_same_dense {
  enum { value = has_direct_access<T1>::ret && has_direct_access<T2>::ret && is_same<typename T1::Scalar,typename T2::Scalar>::value };
};

template<typename T1, typename T2>
EIGEN_DEVICE_FUNC
bool is_same_dense(const T1 &mat1, const T2 &mat2, typename enable_if<possibly_same_dense<T1,T2>::value>::type * = 0)
{
  return (mat1.data()==mat2.data()) && (mat1.innerStride()==mat2.innerStride()) && (mat1.outerStride()==mat2.outerStride());
}

template<typename T1, typename T2>
EIGEN_DEVICE_FUNC
bool is_same_dense(const T1 &, const T2 &, typename enable_if<!possibly_same_dense<T1,T2>::value>::type * = 0)
{
  return false;
}

template<typename T,bool Vectorized=false,typename EnableIf = void>
struct scalar_div_cost {
  enum { value = 8*NumTraits<T>::MulCost };
};

template<typename T,bool Vectorized>
struct scalar_div_cost<std::complex<T>, Vectorized> {
  enum { value = 2*scalar_div_cost<T>::value
               + 6*NumTraits<T>::MulCost
               + 3*NumTraits<T>::AddCost
  };
};

template<bool Vectorized>
struct scalar_div_cost<signed long,Vectorized,typename conditional<sizeof(long)==8,void,false_type>::type> { enum { value = 24 }; };
template<bool Vectorized>
struct scalar_div_cost<unsigned long,Vectorized,typename conditional<sizeof(long)==8,void,false_type>::type> { enum { value = 21 }; };

#ifdef EIGEN_DEBUG_ASSIGN
std::string demangle_traversal(int t)
{
  if(t==DefaultTraversal) return "DefaultTraversal";
  if(t==LinearTraversal) return "LinearTraversal";
  if(t==InnerVectorizedTraversal) return "InnerVectorizedTraversal";
  if(t==LinearVectorizedTraversal) return "LinearVectorizedTraversal";
  if(t==SliceVectorizedTraversal) return "SliceVectorizedTraversal";
  return "?";
}
std::string demangle_unrolling(int t)
{
  if(t==NoUnrolling) return "NoUnrolling";
  if(t==InnerUnrolling) return "InnerUnrolling";
  if(t==CompleteUnrolling) return "CompleteUnrolling";
  return "?";
}
std::string demangle_flags(int f)
{
  std::string res;
  if(f&RowMajorBit)                 res += " | RowMajor";
  if(f&PacketAccessBit)             res += " | Packet";
  if(f&LinearAccessBit)             res += " | Linear";
  if(f&LvalueBit)                   res += " | Lvalue";
  if(f&DirectAccessBit)             res += " | Direct";
  if(f&NestByRefBit)                res += " | NestByRef";
  if(f&NoPreferredStorageOrderBit)  res += " | NoPreferredStorageOrderBit";

  return res;
}
#endif

} 

template<typename ScalarA, typename ScalarB, typename BinaryOp=internal::scalar_product_op<ScalarA,ScalarB> >
struct ScalarBinaryOpTraits
#ifndef EIGEN_PARSED_BY_DOXYGEN
  
  : internal::scalar_product_traits<ScalarA,ScalarB>
#endif 
{};

template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<T,T,BinaryOp>
{
  typedef T ReturnType;
};

template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<T, typename NumTraits<typename internal::enable_if<NumTraits<T>::IsComplex,T>::type>::Real, BinaryOp>
{
  typedef T ReturnType;
};
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<typename NumTraits<typename internal::enable_if<NumTraits<T>::IsComplex,T>::type>::Real, T, BinaryOp>
{
  typedef T ReturnType;
};

template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<T,void,BinaryOp>
{
  typedef T ReturnType;
};

template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<void,T,BinaryOp>
{
  typedef T ReturnType;
};

template<typename BinaryOp>
struct ScalarBinaryOpTraits<void,void,BinaryOp>
{
  typedef void ReturnType;
};

#define EIGEN_CHECK_BINARY_COMPATIBILIY(BINOP,LHS,RHS) \
  EIGEN_STATIC_ASSERT((Eigen::internal::has_ReturnType<ScalarBinaryOpTraits<LHS, RHS,BINOP> >::value), \
    YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)

} 

#endif 
