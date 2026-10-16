#ifndef EIGEN_DENSEBASE_H
#define EIGEN_DENSEBASE_H

namespace Eigen {

namespace internal {

static inline void check_DenseIndex_is_signed() {
  EIGEN_STATIC_ASSERT(NumTraits<DenseIndex>::IsSigned,THE_INDEX_TYPE_MUST_BE_A_SIGNED_TYPE)
}

} 

template<typename Derived> class DenseBase
#ifndef EIGEN_PARSED_BY_DOXYGEN
  : public DenseCoeffsBase<Derived, internal::accessors_level<Derived>::value>
#else
  : public DenseCoeffsBase<Derived,DirectWriteAccessors>
#endif 
{
  public:

    typedef Eigen::InnerIterator<Derived> InnerIterator;

    typedef typename internal::traits<Derived>::StorageKind StorageKind;

    typedef typename internal::traits<Derived>::StorageIndex StorageIndex;

    typedef typename internal::traits<Derived>::Scalar Scalar;

    typedef Scalar value_type;

    typedef typename NumTraits<Scalar>::Real RealScalar;
    typedef DenseCoeffsBase<Derived, internal::accessors_level<Derived>::value> Base;

    using Base::derived;
    using Base::const_cast_derived;
    using Base::rows;
    using Base::cols;
    using Base::size;
    using Base::rowIndexByOuterInner;
    using Base::colIndexByOuterInner;
    using Base::coeff;
    using Base::coeffByOuterInner;
    using Base::operator();
    using Base::operator[];
    using Base::x;
    using Base::y;
    using Base::z;
    using Base::w;
    using Base::stride;
    using Base::innerStride;
    using Base::outerStride;
    using Base::rowStride;
    using Base::colStride;
    typedef typename Base::CoeffReturnType CoeffReturnType;

    enum {

      RowsAtCompileTime = internal::traits<Derived>::RowsAtCompileTime,

      ColsAtCompileTime = internal::traits<Derived>::ColsAtCompileTime,

      SizeAtCompileTime = (internal::size_at_compile_time<internal::traits<Derived>::RowsAtCompileTime,
                                                   internal::traits<Derived>::ColsAtCompileTime>::ret),

      MaxRowsAtCompileTime = internal::traits<Derived>::MaxRowsAtCompileTime,

      MaxColsAtCompileTime = internal::traits<Derived>::MaxColsAtCompileTime,

      MaxSizeAtCompileTime = (internal::size_at_compile_time<internal::traits<Derived>::MaxRowsAtCompileTime,
                                                      internal::traits<Derived>::MaxColsAtCompileTime>::ret),

      IsVectorAtCompileTime = internal::traits<Derived>::RowsAtCompileTime == 1
                           || internal::traits<Derived>::ColsAtCompileTime == 1,

      NumDimensions = int(MaxSizeAtCompileTime) == 1 ? 0 : bool(IsVectorAtCompileTime) ? 1 : 2,

      Flags = internal::traits<Derived>::Flags,

      IsRowMajor = int(Flags) & RowMajorBit, 

      InnerSizeAtCompileTime = int(IsVectorAtCompileTime) ? int(SizeAtCompileTime)
                             : int(IsRowMajor) ? int(ColsAtCompileTime) : int(RowsAtCompileTime),

      InnerStrideAtCompileTime = internal::inner_stride_at_compile_time<Derived>::ret,
      OuterStrideAtCompileTime = internal::outer_stride_at_compile_time<Derived>::ret
    };

    typedef typename internal::find_best_packet<Scalar,SizeAtCompileTime>::type PacketScalar;

    enum { IsPlainObjectBase = 0 };

    typedef Matrix<typename internal::traits<Derived>::Scalar,
                internal::traits<Derived>::RowsAtCompileTime,
                internal::traits<Derived>::ColsAtCompileTime,
                AutoAlign | (internal::traits<Derived>::Flags&RowMajorBit ? RowMajor : ColMajor),
                internal::traits<Derived>::MaxRowsAtCompileTime,
                internal::traits<Derived>::MaxColsAtCompileTime
          > PlainMatrix;

    typedef Array<typename internal::traits<Derived>::Scalar,
                internal::traits<Derived>::RowsAtCompileTime,
                internal::traits<Derived>::ColsAtCompileTime,
                AutoAlign | (internal::traits<Derived>::Flags&RowMajorBit ? RowMajor : ColMajor),
                internal::traits<Derived>::MaxRowsAtCompileTime,
                internal::traits<Derived>::MaxColsAtCompileTime
          > PlainArray;

    typedef typename internal::conditional<internal::is_same<typename internal::traits<Derived>::XprKind,MatrixXpr >::value,
                                 PlainMatrix, PlainArray>::type PlainObject;

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index nonZeros() const { return size(); }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    Index outerSize() const
    {
      return IsVectorAtCompileTime ? 1
           : int(IsRowMajor) ? this->rows() : this->cols();
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    Index innerSize() const
    {
      return IsVectorAtCompileTime ? this->size()
           : int(IsRowMajor) ? this->cols() : this->rows();
    }

    EIGEN_DEVICE_FUNC
    void resize(Index newSize)
    {
      EIGEN_ONLY_USED_FOR_DEBUG(newSize);
      eigen_assert(newSize == this->size()
                && "DenseBase::resize() does not actually allow to resize.");
    }

    EIGEN_DEVICE_FUNC
    void resize(Index rows, Index cols)
    {
      EIGEN_ONLY_USED_FOR_DEBUG(rows);
      EIGEN_ONLY_USED_FOR_DEBUG(cols);
      eigen_assert(rows == this->rows() && cols == this->cols()
                && "DenseBase::resize() does not actually allow to resize.");
    }

#ifndef EIGEN_PARSED_BY_DOXYGEN
    
    typedef CwiseNullaryOp<internal::scalar_constant_op<Scalar>,PlainObject> ConstantReturnType;
    
    EIGEN_DEPRECATED typedef CwiseNullaryOp<internal::linspaced_op<Scalar>,PlainObject> SequentialLinSpacedReturnType;
    
    typedef CwiseNullaryOp<internal::linspaced_op<Scalar>,PlainObject> RandomAccessLinSpacedReturnType;
    
    typedef Matrix<typename NumTraits<typename internal::traits<Derived>::Scalar>::Real, internal::traits<Derived>::ColsAtCompileTime, 1> EigenvaluesReturnType;

#endif 

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator=(const DenseBase<OtherDerived>& other);

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator=(const DenseBase& other);

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    Derived& operator=(const EigenBase<OtherDerived> &other);

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    Derived& operator+=(const EigenBase<OtherDerived> &other);

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    Derived& operator-=(const EigenBase<OtherDerived> &other);

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    Derived& operator=(const ReturnByValue<OtherDerived>& func);

    template<typename OtherDerived>
    
    EIGEN_DEPRECATED EIGEN_DEVICE_FUNC
    Derived& lazyAssign(const DenseBase<OtherDerived>& other);

    EIGEN_DEVICE_FUNC
    CommaInitializer<Derived> operator<< (const Scalar& s);

    template<unsigned int Added,unsigned int Removed>
    
    EIGEN_DEPRECATED
    const Derived& flagged() const
    { return derived(); }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    CommaInitializer<Derived> operator<< (const DenseBase<OtherDerived>& other);

    typedef Transpose<Derived> TransposeReturnType;
    EIGEN_DEVICE_FUNC
    TransposeReturnType transpose();
    typedef typename internal::add_const<Transpose<const Derived> >::type ConstTransposeReturnType;
    EIGEN_DEVICE_FUNC
    ConstTransposeReturnType transpose() const;
    EIGEN_DEVICE_FUNC
    void transposeInPlace();

    EIGEN_DEVICE_FUNC static const ConstantReturnType
    Constant(Index rows, Index cols, const Scalar& value);
    EIGEN_DEVICE_FUNC static const ConstantReturnType
    Constant(Index size, const Scalar& value);
    EIGEN_DEVICE_FUNC static const ConstantReturnType
    Constant(const Scalar& value);

    EIGEN_DEPRECATED EIGEN_DEVICE_FUNC static const RandomAccessLinSpacedReturnType
    LinSpaced(Sequential_t, Index size, const Scalar& low, const Scalar& high);
    EIGEN_DEPRECATED EIGEN_DEVICE_FUNC static const RandomAccessLinSpacedReturnType
    LinSpaced(Sequential_t, const Scalar& low, const Scalar& high);

    EIGEN_DEVICE_FUNC static const RandomAccessLinSpacedReturnType
    LinSpaced(Index size, const Scalar& low, const Scalar& high);
    EIGEN_DEVICE_FUNC static const RandomAccessLinSpacedReturnType
    LinSpaced(const Scalar& low, const Scalar& high);

    template<typename CustomNullaryOp> EIGEN_DEVICE_FUNC
    static const CwiseNullaryOp<CustomNullaryOp, PlainObject>
    NullaryExpr(Index rows, Index cols, const CustomNullaryOp& func);
    template<typename CustomNullaryOp> EIGEN_DEVICE_FUNC
    static const CwiseNullaryOp<CustomNullaryOp, PlainObject>
    NullaryExpr(Index size, const CustomNullaryOp& func);
    template<typename CustomNullaryOp> EIGEN_DEVICE_FUNC
    static const CwiseNullaryOp<CustomNullaryOp, PlainObject>
    NullaryExpr(const CustomNullaryOp& func);

    EIGEN_DEVICE_FUNC static const ConstantReturnType Zero(Index rows, Index cols);
    EIGEN_DEVICE_FUNC static const ConstantReturnType Zero(Index size);
    EIGEN_DEVICE_FUNC static const ConstantReturnType Zero();
    EIGEN_DEVICE_FUNC static const ConstantReturnType Ones(Index rows, Index cols);
    EIGEN_DEVICE_FUNC static const ConstantReturnType Ones(Index size);
    EIGEN_DEVICE_FUNC static const ConstantReturnType Ones();

    EIGEN_DEVICE_FUNC void fill(const Scalar& value);
    EIGEN_DEVICE_FUNC Derived& setConstant(const Scalar& value);
    EIGEN_DEVICE_FUNC Derived& setLinSpaced(Index size, const Scalar& low, const Scalar& high);
    EIGEN_DEVICE_FUNC Derived& setLinSpaced(const Scalar& low, const Scalar& high);
    EIGEN_DEVICE_FUNC Derived& setZero();
    EIGEN_DEVICE_FUNC Derived& setOnes();
    EIGEN_DEVICE_FUNC Derived& setRandom();

    template<typename OtherDerived> EIGEN_DEVICE_FUNC
    bool isApprox(const DenseBase<OtherDerived>& other,
                  const RealScalar& prec = NumTraits<Scalar>::dummy_precision()) const;
    EIGEN_DEVICE_FUNC
    bool isMuchSmallerThan(const RealScalar& other,
                           const RealScalar& prec = NumTraits<Scalar>::dummy_precision()) const;
    template<typename OtherDerived> EIGEN_DEVICE_FUNC
    bool isMuchSmallerThan(const DenseBase<OtherDerived>& other,
                           const RealScalar& prec = NumTraits<Scalar>::dummy_precision()) const;

    EIGEN_DEVICE_FUNC bool isApproxToConstant(const Scalar& value, const RealScalar& prec = NumTraits<Scalar>::dummy_precision()) const;
    EIGEN_DEVICE_FUNC bool isConstant(const Scalar& value, const RealScalar& prec = NumTraits<Scalar>::dummy_precision()) const;
    EIGEN_DEVICE_FUNC bool isZero(const RealScalar& prec = NumTraits<Scalar>::dummy_precision()) const;
    EIGEN_DEVICE_FUNC bool isOnes(const RealScalar& prec = NumTraits<Scalar>::dummy_precision()) const;

    inline bool hasNaN() const;
    inline bool allFinite() const;

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator*=(const Scalar& other);
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    Derived& operator/=(const Scalar& other);

    typedef typename internal::add_const_on_value_type<typename internal::eval<Derived>::type>::type EvalReturnType;

    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE EvalReturnType eval() const
    {

      return typename internal::eval<Derived>::type(derived());
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    void swap(const DenseBase<OtherDerived>& other)
    {
      EIGEN_STATIC_ASSERT(!OtherDerived::IsPlainObjectBase,THIS_EXPRESSION_IS_NOT_A_LVALUE__IT_IS_READ_ONLY);
      eigen_assert(rows()==other.rows() && cols()==other.cols());
      call_assignment(derived(), other.const_cast_derived(), internal::swap_assign_op<Scalar>());
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    void swap(PlainObjectBase<OtherDerived>& other)
    {
      eigen_assert(rows()==other.rows() && cols()==other.cols());
      call_assignment(derived(), other.derived(), internal::swap_assign_op<Scalar>());
    }

    EIGEN_DEVICE_FUNC inline const NestByValue<Derived> nestByValue() const;
    EIGEN_DEVICE_FUNC inline const ForceAlignedAccess<Derived> forceAlignedAccess() const;
    EIGEN_DEVICE_FUNC inline ForceAlignedAccess<Derived> forceAlignedAccess();
    template<bool Enable> EIGEN_DEVICE_FUNC
    inline const typename internal::conditional<Enable,ForceAlignedAccess<Derived>,Derived&>::type forceAlignedAccessIf() const;
    template<bool Enable> EIGEN_DEVICE_FUNC
    inline typename internal::conditional<Enable,ForceAlignedAccess<Derived>,Derived&>::type forceAlignedAccessIf();

    EIGEN_DEVICE_FUNC Scalar sum() const;
    EIGEN_DEVICE_FUNC Scalar mean() const;
    EIGEN_DEVICE_FUNC Scalar trace() const;

    EIGEN_DEVICE_FUNC Scalar prod() const;

    template<int NaNPropagation>
    EIGEN_DEVICE_FUNC typename internal::traits<Derived>::Scalar minCoeff() const;
    template<int NaNPropagation>
    EIGEN_DEVICE_FUNC typename internal::traits<Derived>::Scalar maxCoeff() const;

    EIGEN_DEVICE_FUNC inline typename internal::traits<Derived>::Scalar minCoeff() const {
      return minCoeff<PropagateFast>();
    }
    EIGEN_DEVICE_FUNC inline typename internal::traits<Derived>::Scalar maxCoeff() const {
      return maxCoeff<PropagateFast>();
    }

    template<int NaNPropagation, typename IndexType>
    EIGEN_DEVICE_FUNC
    typename internal::traits<Derived>::Scalar minCoeff(IndexType* row, IndexType* col) const;
    template<int NaNPropagation, typename IndexType>
    EIGEN_DEVICE_FUNC
    typename internal::traits<Derived>::Scalar maxCoeff(IndexType* row, IndexType* col) const;
    template<int NaNPropagation, typename IndexType>
    EIGEN_DEVICE_FUNC
    typename internal::traits<Derived>::Scalar minCoeff(IndexType* index) const;
    template<int NaNPropagation, typename IndexType>
    EIGEN_DEVICE_FUNC
    typename internal::traits<Derived>::Scalar maxCoeff(IndexType* index) const;

    template<typename IndexType>
    EIGEN_DEVICE_FUNC inline
    typename internal::traits<Derived>::Scalar minCoeff(IndexType* row, IndexType* col) const {
      return minCoeff<PropagateFast>(row, col);
    }
    template<typename IndexType>
    EIGEN_DEVICE_FUNC inline
    typename internal::traits<Derived>::Scalar maxCoeff(IndexType* row, IndexType* col) const {
      return maxCoeff<PropagateFast>(row, col);
    }
    template<typename IndexType>
     EIGEN_DEVICE_FUNC inline
    typename internal::traits<Derived>::Scalar minCoeff(IndexType* index) const {
      return minCoeff<PropagateFast>(index);
    }
    template<typename IndexType>
    EIGEN_DEVICE_FUNC inline
    typename internal::traits<Derived>::Scalar maxCoeff(IndexType* index) const {
      return maxCoeff<PropagateFast>(index);
    }
  
    template<typename BinaryOp>
    EIGEN_DEVICE_FUNC
    Scalar redux(const BinaryOp& func) const;

    template<typename Visitor>
    EIGEN_DEVICE_FUNC
    void visit(Visitor& func) const;

    inline const WithFormat<Derived> format(const IOFormat& fmt) const
    {
      return WithFormat<Derived>(derived(), fmt);
    }

    EIGEN_DEVICE_FUNC
    CoeffReturnType value() const
    {
      EIGEN_STATIC_ASSERT_SIZE_1x1(Derived)
      eigen_assert(this->rows() == 1 && this->cols() == 1);
      return derived().coeff(0,0);
    }

    EIGEN_DEVICE_FUNC bool all() const;
    EIGEN_DEVICE_FUNC bool any() const;
    EIGEN_DEVICE_FUNC Index count() const;

    typedef VectorwiseOp<Derived, Horizontal> RowwiseReturnType;
    typedef const VectorwiseOp<const Derived, Horizontal> ConstRowwiseReturnType;
    typedef VectorwiseOp<Derived, Vertical> ColwiseReturnType;
    typedef const VectorwiseOp<const Derived, Vertical> ConstColwiseReturnType;

    EIGEN_DEVICE_FUNC inline ConstRowwiseReturnType rowwise() const {
      return ConstRowwiseReturnType(derived());
    }
    EIGEN_DEVICE_FUNC RowwiseReturnType rowwise();

    EIGEN_DEVICE_FUNC inline ConstColwiseReturnType colwise() const {
      return ConstColwiseReturnType(derived());
    }
    EIGEN_DEVICE_FUNC ColwiseReturnType colwise();

    typedef CwiseNullaryOp<internal::scalar_random_op<Scalar>,PlainObject> RandomReturnType;
    static const RandomReturnType Random(Index rows, Index cols);
    static const RandomReturnType Random(Index size);
    static const RandomReturnType Random();

    template<typename ThenDerived,typename ElseDerived>
    inline EIGEN_DEVICE_FUNC const Select<Derived,ThenDerived,ElseDerived>
    select(const DenseBase<ThenDerived>& thenMatrix,
           const DenseBase<ElseDerived>& elseMatrix) const;

    template<typename ThenDerived>
    inline EIGEN_DEVICE_FUNC const Select<Derived,ThenDerived, typename ThenDerived::ConstantReturnType>
    select(const DenseBase<ThenDerived>& thenMatrix, const typename ThenDerived::Scalar& elseScalar) const;

    template<typename ElseDerived>
    inline EIGEN_DEVICE_FUNC const Select<Derived, typename ElseDerived::ConstantReturnType, ElseDerived >
    select(const typename ElseDerived::Scalar& thenScalar, const DenseBase<ElseDerived>& elseMatrix) const;

    template<int p> RealScalar lpNorm() const;

    template<int RowFactor, int ColFactor>
    EIGEN_DEVICE_FUNC
    const Replicate<Derived,RowFactor,ColFactor> replicate() const;

    EIGEN_DEVICE_FUNC
    const Replicate<Derived, Dynamic, Dynamic> replicate(Index rowFactor, Index colFactor) const
    {
      return Replicate<Derived, Dynamic, Dynamic>(derived(), rowFactor, colFactor);
    }

    typedef Reverse<Derived, BothDirections> ReverseReturnType;
    typedef const Reverse<const Derived, BothDirections> ConstReverseReturnType;
    EIGEN_DEVICE_FUNC ReverseReturnType reverse();

    EIGEN_DEVICE_FUNC ConstReverseReturnType reverse() const
    {
      return ConstReverseReturnType(derived());
    }
    EIGEN_DEVICE_FUNC void reverseInPlace();

    #ifdef EIGEN_PARSED_BY_DOXYGEN

    typedef random_access_iterator_type iterator;
    
    typedef random_access_iterator_type const_iterator;
    #else
    typedef typename internal::conditional< (Flags&DirectAccessBit)==DirectAccessBit,
                                            internal::pointer_based_stl_iterator<Derived>,
                                            internal::generic_randaccess_stl_iterator<Derived>
                                          >::type iterator_type;

    typedef typename internal::conditional< (Flags&DirectAccessBit)==DirectAccessBit,
                                            internal::pointer_based_stl_iterator<const Derived>,
                                            internal::generic_randaccess_stl_iterator<const Derived>
                                          >::type const_iterator_type;

    typedef typename internal::conditional< IsVectorAtCompileTime,
                                            iterator_type,
                                            void
                                          >::type iterator;

    typedef typename internal::conditional< IsVectorAtCompileTime,
                                            const_iterator_type,
                                            void
                                          >::type const_iterator;
    #endif

    inline iterator begin();
    inline const_iterator begin() const;
    inline const_iterator cbegin() const;
    inline iterator end();
    inline const_iterator end() const;
    inline const_iterator cend() const;

#define EIGEN_CURRENT_STORAGE_BASE_CLASS Eigen::DenseBase
#define EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL
#define EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(COND)
#define EIGEN_DOC_UNARY_ADDONS(X,Y)
#   include "../plugins/CommonCwiseUnaryOps.h"
#   include "../plugins/BlockMethods.h"
#   include "../plugins/IndexedViewMethods.h"
#   include "../plugins/ReshapedMethods.h"
#   ifdef EIGEN_DENSEBASE_PLUGIN
#     include EIGEN_DENSEBASE_PLUGIN
#   endif
#undef EIGEN_CURRENT_STORAGE_BASE_CLASS
#undef EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL
#undef EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF
#undef EIGEN_DOC_UNARY_ADDONS

    template<typename Dest>
    EIGEN_DEVICE_FUNC
    inline void evalTo(Dest& ) const
    {
      EIGEN_STATIC_ASSERT((internal::is_same<Dest,void>::value),THE_EVAL_EVALTO_FUNCTION_SHOULD_NEVER_BE_CALLED_FOR_DENSE_OBJECTS);
    }

  protected:
    EIGEN_DEFAULT_COPY_CONSTRUCTOR(DenseBase)
    
    EIGEN_DEVICE_FUNC DenseBase()
    {

#ifdef EIGEN_INTERNAL_DEBUGGING
      EIGEN_STATIC_ASSERT((EIGEN_IMPLIES(MaxRowsAtCompileTime==1 && MaxColsAtCompileTime!=1, int(IsRowMajor))
                        && EIGEN_IMPLIES(MaxColsAtCompileTime==1 && MaxRowsAtCompileTime!=1, int(!IsRowMajor))),
                          INVALID_STORAGE_ORDER_FOR_THIS_VECTOR_EXPRESSION)
#endif
    }

  private:
    EIGEN_DEVICE_FUNC explicit DenseBase(int);
    EIGEN_DEVICE_FUNC DenseBase(int,int);
    template<typename OtherDerived> EIGEN_DEVICE_FUNC explicit DenseBase(const DenseBase<OtherDerived>&);
};

} 

#endif 
