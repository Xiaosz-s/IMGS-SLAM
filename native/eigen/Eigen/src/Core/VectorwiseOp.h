#ifndef EIGEN_PARTIAL_REDUX_H
#define EIGEN_PARTIAL_REDUX_H

namespace Eigen {

template< typename MatrixType, typename MemberOp, int Direction>
class PartialReduxExpr;

namespace internal {
template<typename MatrixType, typename MemberOp, int Direction>
struct traits<PartialReduxExpr<MatrixType, MemberOp, Direction> >
 : traits<MatrixType>
{
  typedef typename MemberOp::result_type Scalar;
  typedef typename traits<MatrixType>::StorageKind StorageKind;
  typedef typename traits<MatrixType>::XprKind XprKind;
  typedef typename MatrixType::Scalar InputScalar;
  enum {
    RowsAtCompileTime = Direction==Vertical   ? 1 : MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = Direction==Horizontal ? 1 : MatrixType::ColsAtCompileTime,
    MaxRowsAtCompileTime = Direction==Vertical   ? 1 : MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = Direction==Horizontal ? 1 : MatrixType::MaxColsAtCompileTime,
    Flags = RowsAtCompileTime == 1 ? RowMajorBit : 0,
    TraversalSize = Direction==Vertical ? MatrixType::RowsAtCompileTime :  MatrixType::ColsAtCompileTime
  };
};
}

template< typename MatrixType, typename MemberOp, int Direction>
class PartialReduxExpr : public internal::dense_xpr_base< PartialReduxExpr<MatrixType, MemberOp, Direction> >::type,
                         internal::no_assignment_operator
{
  public:

    typedef typename internal::dense_xpr_base<PartialReduxExpr>::type Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(PartialReduxExpr)

    EIGEN_DEVICE_FUNC
    explicit PartialReduxExpr(const MatrixType& mat, const MemberOp& func = MemberOp())
      : m_matrix(mat), m_functor(func) {}

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    Index rows() const EIGEN_NOEXCEPT { return (Direction==Vertical   ? 1 : m_matrix.rows()); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    Index cols() const EIGEN_NOEXCEPT { return (Direction==Horizontal ? 1 : m_matrix.cols()); }

    EIGEN_DEVICE_FUNC
    typename MatrixType::Nested nestedExpression() const { return m_matrix; }

    EIGEN_DEVICE_FUNC
    const MemberOp& functor() const { return m_functor; }

  protected:
    typename MatrixType::Nested m_matrix;
    const MemberOp m_functor;
};

template<typename A,typename B> struct partial_redux_dummy_func;

#define EIGEN_MAKE_PARTIAL_REDUX_FUNCTOR(MEMBER,COST,VECTORIZABLE,BINARYOP)                \
  template <typename ResultType,typename Scalar>                                                            \
  struct member_##MEMBER {                                                                  \
    EIGEN_EMPTY_STRUCT_CTOR(member_##MEMBER)                                                \
    typedef ResultType result_type;                                                         \
    typedef BINARYOP<Scalar,Scalar> BinaryOp;   \
    template<int Size> struct Cost { enum { value = COST }; };             \
    enum { Vectorizable = VECTORIZABLE };                                                   \
    template<typename XprType>                                                              \
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE                                                   \
    ResultType operator()(const XprType& mat) const                                         \
    { return mat.MEMBER(); }                                                                \
    BinaryOp binaryFunc() const { return BinaryOp(); }                                      \
  }

#define EIGEN_MEMBER_FUNCTOR(MEMBER,COST) \
  EIGEN_MAKE_PARTIAL_REDUX_FUNCTOR(MEMBER,COST,0,partial_redux_dummy_func)

namespace internal {

EIGEN_MEMBER_FUNCTOR(norm, (Size+5) * NumTraits<Scalar>::MulCost + (Size-1)*NumTraits<Scalar>::AddCost);
EIGEN_MEMBER_FUNCTOR(stableNorm, (Size+5) * NumTraits<Scalar>::MulCost + (Size-1)*NumTraits<Scalar>::AddCost);
EIGEN_MEMBER_FUNCTOR(blueNorm, (Size+5) * NumTraits<Scalar>::MulCost + (Size-1)*NumTraits<Scalar>::AddCost);
EIGEN_MEMBER_FUNCTOR(hypotNorm, (Size-1) * functor_traits<scalar_hypot_op<Scalar> >::Cost );
EIGEN_MEMBER_FUNCTOR(all, (Size-1)*NumTraits<Scalar>::AddCost);
EIGEN_MEMBER_FUNCTOR(any, (Size-1)*NumTraits<Scalar>::AddCost);
EIGEN_MEMBER_FUNCTOR(count, (Size-1)*NumTraits<Scalar>::AddCost);

EIGEN_MAKE_PARTIAL_REDUX_FUNCTOR(sum, (Size-1)*NumTraits<Scalar>::AddCost, 1, internal::scalar_sum_op);
EIGEN_MAKE_PARTIAL_REDUX_FUNCTOR(minCoeff, (Size-1)*NumTraits<Scalar>::AddCost, 1, internal::scalar_min_op);
EIGEN_MAKE_PARTIAL_REDUX_FUNCTOR(maxCoeff, (Size-1)*NumTraits<Scalar>::AddCost, 1, internal::scalar_max_op);
EIGEN_MAKE_PARTIAL_REDUX_FUNCTOR(prod, (Size-1)*NumTraits<Scalar>::MulCost, 1, internal::scalar_product_op);

template <int p, typename ResultType,typename Scalar>
struct member_lpnorm {
  typedef ResultType result_type;
  enum { Vectorizable = 0 };
  template<int Size> struct Cost
  { enum { value = (Size+5) * NumTraits<Scalar>::MulCost + (Size-1)*NumTraits<Scalar>::AddCost }; };
  EIGEN_DEVICE_FUNC member_lpnorm() {}
  template<typename XprType>
  EIGEN_DEVICE_FUNC inline ResultType operator()(const XprType& mat) const
  { return mat.template lpNorm<p>(); }
};

template <typename BinaryOpT, typename Scalar>
struct member_redux {
  typedef BinaryOpT BinaryOp;
  typedef typename result_of<
                     BinaryOp(const Scalar&,const Scalar&)
                   >::type  result_type;

  enum { Vectorizable = functor_traits<BinaryOp>::PacketAccess };
  template<int Size> struct Cost { enum { value = (Size-1) * functor_traits<BinaryOp>::Cost }; };
  EIGEN_DEVICE_FUNC explicit member_redux(const BinaryOp func) : m_functor(func) {}
  template<typename Derived>
  EIGEN_DEVICE_FUNC inline result_type operator()(const DenseBase<Derived>& mat) const
  { return mat.redux(m_functor); }
  const BinaryOp& binaryFunc() const { return m_functor; }
  const BinaryOp m_functor;
};
}

template<typename ExpressionType, int Direction> class VectorwiseOp
{
  public:

    typedef typename ExpressionType::Scalar Scalar;
    typedef typename ExpressionType::RealScalar RealScalar;
    typedef Eigen::Index Index; 
    typedef typename internal::ref_selector<ExpressionType>::non_const_type ExpressionTypeNested;
    typedef typename internal::remove_all<ExpressionTypeNested>::type ExpressionTypeNestedCleaned;

    template<template<typename OutScalar,typename InputScalar> class Functor,
                      typename ReturnScalar=Scalar> struct ReturnType
    {
      typedef PartialReduxExpr<ExpressionType,
                               Functor<ReturnScalar,Scalar>,
                               Direction
                              > Type;
    };

    template<typename BinaryOp> struct ReduxReturnType
    {
      typedef PartialReduxExpr<ExpressionType,
                               internal::member_redux<BinaryOp,Scalar>,
                               Direction
                              > Type;
    };

    enum {
      isVertical   = (Direction==Vertical) ? 1 : 0,
      isHorizontal = (Direction==Horizontal) ? 1 : 0
    };

  protected:

    template<typename OtherDerived> struct ExtendedType {
      typedef Replicate<OtherDerived,
                        isVertical   ? 1 : ExpressionType::RowsAtCompileTime,
                        isHorizontal ? 1 : ExpressionType::ColsAtCompileTime> Type;
    };

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    typename ExtendedType<OtherDerived>::Type
    extendedTo(const DenseBase<OtherDerived>& other) const
    {
      EIGEN_STATIC_ASSERT(EIGEN_IMPLIES(isVertical, OtherDerived::MaxColsAtCompileTime==1),
                          YOU_PASSED_A_ROW_VECTOR_BUT_A_COLUMN_VECTOR_WAS_EXPECTED)
      EIGEN_STATIC_ASSERT(EIGEN_IMPLIES(isHorizontal, OtherDerived::MaxRowsAtCompileTime==1),
                          YOU_PASSED_A_COLUMN_VECTOR_BUT_A_ROW_VECTOR_WAS_EXPECTED)
      return typename ExtendedType<OtherDerived>::Type
                      (other.derived(),
                       isVertical   ? 1 : m_matrix.rows(),
                       isHorizontal ? 1 : m_matrix.cols());
    }

    template<typename OtherDerived> struct OppositeExtendedType {
      typedef Replicate<OtherDerived,
                        isHorizontal ? 1 : ExpressionType::RowsAtCompileTime,
                        isVertical   ? 1 : ExpressionType::ColsAtCompileTime> Type;
    };

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    typename OppositeExtendedType<OtherDerived>::Type
    extendedToOpposite(const DenseBase<OtherDerived>& other) const
    {
      EIGEN_STATIC_ASSERT(EIGEN_IMPLIES(isHorizontal, OtherDerived::MaxColsAtCompileTime==1),
                          YOU_PASSED_A_ROW_VECTOR_BUT_A_COLUMN_VECTOR_WAS_EXPECTED)
      EIGEN_STATIC_ASSERT(EIGEN_IMPLIES(isVertical, OtherDerived::MaxRowsAtCompileTime==1),
                          YOU_PASSED_A_COLUMN_VECTOR_BUT_A_ROW_VECTOR_WAS_EXPECTED)
      return typename OppositeExtendedType<OtherDerived>::Type
                      (other.derived(),
                       isHorizontal  ? 1 : m_matrix.rows(),
                       isVertical    ? 1 : m_matrix.cols());
    }

  public:
    EIGEN_DEVICE_FUNC
    explicit inline VectorwiseOp(ExpressionType& matrix) : m_matrix(matrix) {}

    EIGEN_DEVICE_FUNC
    inline const ExpressionType& _expression() const { return m_matrix; }

    #ifdef EIGEN_PARSED_BY_DOXYGEN

    random_access_iterator_type iterator;
    
    random_access_iterator_type const_iterator;
    #else
    typedef internal::subvector_stl_iterator<ExpressionType,               DirectionType(Direction)> iterator;
    typedef internal::subvector_stl_iterator<const ExpressionType,         DirectionType(Direction)> const_iterator;
    typedef internal::subvector_stl_reverse_iterator<ExpressionType,       DirectionType(Direction)> reverse_iterator;
    typedef internal::subvector_stl_reverse_iterator<const ExpressionType, DirectionType(Direction)> const_reverse_iterator;
    #endif

    iterator                 begin()       { return iterator      (m_matrix, 0); }
    
    const_iterator           begin() const { return const_iterator(m_matrix, 0); }
    
    const_iterator          cbegin() const { return const_iterator(m_matrix, 0); }

    reverse_iterator        rbegin()       { return reverse_iterator       (m_matrix, m_matrix.template subVectors<DirectionType(Direction)>()-1); }
	
    const_reverse_iterator  rbegin() const { return const_reverse_iterator (m_matrix, m_matrix.template subVectors<DirectionType(Direction)>()-1); }
	
	const_reverse_iterator crbegin() const { return const_reverse_iterator (m_matrix, m_matrix.template subVectors<DirectionType(Direction)>()-1); }

    iterator                 end()         { return iterator      (m_matrix, m_matrix.template subVectors<DirectionType(Direction)>()); }
    
    const_iterator           end()  const  { return const_iterator(m_matrix, m_matrix.template subVectors<DirectionType(Direction)>()); }
    
    const_iterator          cend()  const  { return const_iterator(m_matrix, m_matrix.template subVectors<DirectionType(Direction)>()); }

    reverse_iterator        rend()         { return reverse_iterator       (m_matrix, -1); }
    
    const_reverse_iterator  rend()  const  { return const_reverse_iterator (m_matrix, -1); }
    
    const_reverse_iterator crend()  const  { return const_reverse_iterator (m_matrix, -1); }

    template<typename BinaryOp>
    EIGEN_DEVICE_FUNC
    const typename ReduxReturnType<BinaryOp>::Type
    redux(const BinaryOp& func = BinaryOp()) const
    {
      eigen_assert(redux_length()>0 && "you are using an empty matrix");
      return typename ReduxReturnType<BinaryOp>::Type(_expression(), internal::member_redux<BinaryOp,Scalar>(func));
    }

    typedef typename ReturnType<internal::member_minCoeff>::Type MinCoeffReturnType;
    typedef typename ReturnType<internal::member_maxCoeff>::Type MaxCoeffReturnType;
    typedef PartialReduxExpr<const CwiseUnaryOp<internal::scalar_abs2_op<Scalar>, const ExpressionTypeNestedCleaned>,internal::member_sum<RealScalar,RealScalar>,Direction> SquaredNormReturnType;
    typedef CwiseUnaryOp<internal::scalar_sqrt_op<RealScalar>, const SquaredNormReturnType> NormReturnType;
    typedef typename ReturnType<internal::member_blueNorm,RealScalar>::Type BlueNormReturnType;
    typedef typename ReturnType<internal::member_stableNorm,RealScalar>::Type StableNormReturnType;
    typedef typename ReturnType<internal::member_hypotNorm,RealScalar>::Type HypotNormReturnType;
    typedef typename ReturnType<internal::member_sum>::Type SumReturnType;
    typedef EIGEN_EXPR_BINARYOP_SCALAR_RETURN_TYPE(SumReturnType,Scalar,quotient) MeanReturnType;
    typedef typename ReturnType<internal::member_all>::Type AllReturnType;
    typedef typename ReturnType<internal::member_any>::Type AnyReturnType;
    typedef PartialReduxExpr<ExpressionType, internal::member_count<Index,Scalar>, Direction> CountReturnType;
    typedef typename ReturnType<internal::member_prod>::Type ProdReturnType;
    typedef Reverse<const ExpressionType, Direction> ConstReverseReturnType;
    typedef Reverse<ExpressionType, Direction> ReverseReturnType;

    template<int p> struct LpNormReturnType {
      typedef PartialReduxExpr<ExpressionType, internal::member_lpnorm<p,RealScalar,Scalar>,Direction> Type;
    };

    EIGEN_DEVICE_FUNC
    const MinCoeffReturnType minCoeff() const
    {
      eigen_assert(redux_length()>0 && "you are using an empty matrix");
      return MinCoeffReturnType(_expression());
    }

    EIGEN_DEVICE_FUNC
    const MaxCoeffReturnType maxCoeff() const
    {
      eigen_assert(redux_length()>0 && "you are using an empty matrix");
      return MaxCoeffReturnType(_expression());
    }

    EIGEN_DEVICE_FUNC
    const SquaredNormReturnType squaredNorm() const
    { return SquaredNormReturnType(m_matrix.cwiseAbs2()); }

    EIGEN_DEVICE_FUNC
    const NormReturnType norm() const
    { return NormReturnType(squaredNorm()); }

    template<int p>
    EIGEN_DEVICE_FUNC
    const typename LpNormReturnType<p>::Type lpNorm() const
    { return typename LpNormReturnType<p>::Type(_expression()); }

    EIGEN_DEVICE_FUNC
    const BlueNormReturnType blueNorm() const
    { return BlueNormReturnType(_expression()); }

    EIGEN_DEVICE_FUNC
    const StableNormReturnType stableNorm() const
    { return StableNormReturnType(_expression()); }

    EIGEN_DEVICE_FUNC
    const HypotNormReturnType hypotNorm() const
    { return HypotNormReturnType(_expression()); }

    EIGEN_DEVICE_FUNC
    const SumReturnType sum() const
    { return SumReturnType(_expression()); }

    EIGEN_DEVICE_FUNC
    const MeanReturnType mean() const
    { return sum() / Scalar(Direction==Vertical?m_matrix.rows():m_matrix.cols()); }

    EIGEN_DEVICE_FUNC
    const AllReturnType all() const
    { return AllReturnType(_expression()); }

    EIGEN_DEVICE_FUNC
    const AnyReturnType any() const
    { return AnyReturnType(_expression()); }

    EIGEN_DEVICE_FUNC
    const CountReturnType count() const
    { return CountReturnType(_expression()); }

    EIGEN_DEVICE_FUNC
    const ProdReturnType prod() const
    { return ProdReturnType(_expression()); }

    EIGEN_DEVICE_FUNC
    const ConstReverseReturnType reverse() const
    { return ConstReverseReturnType( _expression() ); }

    EIGEN_DEVICE_FUNC
    ReverseReturnType reverse()
    { return ReverseReturnType( _expression() ); }

    typedef Replicate<ExpressionType,(isVertical?Dynamic:1),(isHorizontal?Dynamic:1)> ReplicateReturnType;
    EIGEN_DEVICE_FUNC
    const ReplicateReturnType replicate(Index factor) const;

    template<int Factor> const Replicate<ExpressionType,isVertical*Factor+isHorizontal,isHorizontal*Factor+isVertical>
    EIGEN_DEVICE_FUNC
    replicate(Index factor = Factor) const
    {
      return Replicate<ExpressionType,(isVertical?Factor:1),(isHorizontal?Factor:1)>
          (_expression(),isVertical?factor:1,isHorizontal?factor:1);
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    ExpressionType& operator=(const DenseBase<OtherDerived>& other)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(OtherDerived)
      EIGEN_STATIC_ASSERT_SAME_XPR_KIND(ExpressionType, OtherDerived)
      
      return m_matrix = extendedTo(other.derived());
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    ExpressionType& operator+=(const DenseBase<OtherDerived>& other)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(OtherDerived)
      EIGEN_STATIC_ASSERT_SAME_XPR_KIND(ExpressionType, OtherDerived)
      return m_matrix += extendedTo(other.derived());
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    ExpressionType& operator-=(const DenseBase<OtherDerived>& other)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(OtherDerived)
      EIGEN_STATIC_ASSERT_SAME_XPR_KIND(ExpressionType, OtherDerived)
      return m_matrix -= extendedTo(other.derived());
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    ExpressionType& operator*=(const DenseBase<OtherDerived>& other)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(OtherDerived)
      EIGEN_STATIC_ASSERT_ARRAYXPR(ExpressionType)
      EIGEN_STATIC_ASSERT_SAME_XPR_KIND(ExpressionType, OtherDerived)
      m_matrix *= extendedTo(other.derived());
      return m_matrix;
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    ExpressionType& operator/=(const DenseBase<OtherDerived>& other)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(OtherDerived)
      EIGEN_STATIC_ASSERT_ARRAYXPR(ExpressionType)
      EIGEN_STATIC_ASSERT_SAME_XPR_KIND(ExpressionType, OtherDerived)
      m_matrix /= extendedTo(other.derived());
      return m_matrix;
    }

    template<typename OtherDerived> EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC
    CwiseBinaryOp<internal::scalar_sum_op<Scalar,typename OtherDerived::Scalar>,
                  const ExpressionTypeNestedCleaned,
                  const typename ExtendedType<OtherDerived>::Type>
    operator+(const DenseBase<OtherDerived>& other) const
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(OtherDerived)
      EIGEN_STATIC_ASSERT_SAME_XPR_KIND(ExpressionType, OtherDerived)
      return m_matrix + extendedTo(other.derived());
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    CwiseBinaryOp<internal::scalar_difference_op<Scalar,typename OtherDerived::Scalar>,
                  const ExpressionTypeNestedCleaned,
                  const typename ExtendedType<OtherDerived>::Type>
    operator-(const DenseBase<OtherDerived>& other) const
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(OtherDerived)
      EIGEN_STATIC_ASSERT_SAME_XPR_KIND(ExpressionType, OtherDerived)
      return m_matrix - extendedTo(other.derived());
    }

    template<typename OtherDerived> EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC
    CwiseBinaryOp<internal::scalar_product_op<Scalar>,
                  const ExpressionTypeNestedCleaned,
                  const typename ExtendedType<OtherDerived>::Type>
    EIGEN_DEVICE_FUNC
    operator*(const DenseBase<OtherDerived>& other) const
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(OtherDerived)
      EIGEN_STATIC_ASSERT_ARRAYXPR(ExpressionType)
      EIGEN_STATIC_ASSERT_SAME_XPR_KIND(ExpressionType, OtherDerived)
      return m_matrix * extendedTo(other.derived());
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    CwiseBinaryOp<internal::scalar_quotient_op<Scalar>,
                  const ExpressionTypeNestedCleaned,
                  const typename ExtendedType<OtherDerived>::Type>
    operator/(const DenseBase<OtherDerived>& other) const
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(OtherDerived)
      EIGEN_STATIC_ASSERT_ARRAYXPR(ExpressionType)
      EIGEN_STATIC_ASSERT_SAME_XPR_KIND(ExpressionType, OtherDerived)
      return m_matrix / extendedTo(other.derived());
    }

    EIGEN_DEVICE_FUNC
    CwiseBinaryOp<internal::scalar_quotient_op<Scalar>,
                  const ExpressionTypeNestedCleaned,
                  const typename OppositeExtendedType<NormReturnType>::Type>
    normalized() const { return m_matrix.cwiseQuotient(extendedToOpposite(this->norm())); }

    EIGEN_DEVICE_FUNC void normalize() {
      m_matrix = this->normalized();
    }

    EIGEN_DEVICE_FUNC inline void reverseInPlace();

    typedef Homogeneous<ExpressionType,Direction> HomogeneousReturnType;
    EIGEN_DEVICE_FUNC
    HomogeneousReturnType homogeneous() const;

    typedef typename ExpressionType::PlainObject CrossReturnType;
    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    const CrossReturnType cross(const MatrixBase<OtherDerived>& other) const;

    enum {
      HNormalized_Size = Direction==Vertical ? internal::traits<ExpressionType>::RowsAtCompileTime
                                             : internal::traits<ExpressionType>::ColsAtCompileTime,
      HNormalized_SizeMinusOne = HNormalized_Size==Dynamic ? Dynamic : HNormalized_Size-1
    };
    typedef Block<const ExpressionType,
                  Direction==Vertical   ? int(HNormalized_SizeMinusOne)
                                        : int(internal::traits<ExpressionType>::RowsAtCompileTime),
                  Direction==Horizontal ? int(HNormalized_SizeMinusOne)
                                        : int(internal::traits<ExpressionType>::ColsAtCompileTime)>
            HNormalized_Block;
    typedef Block<const ExpressionType,
                  Direction==Vertical   ? 1 : int(internal::traits<ExpressionType>::RowsAtCompileTime),
                  Direction==Horizontal ? 1 : int(internal::traits<ExpressionType>::ColsAtCompileTime)>
            HNormalized_Factors;
    typedef CwiseBinaryOp<internal::scalar_quotient_op<typename internal::traits<ExpressionType>::Scalar>,
                const HNormalized_Block,
                const Replicate<HNormalized_Factors,
                  Direction==Vertical   ? HNormalized_SizeMinusOne : 1,
                  Direction==Horizontal ? HNormalized_SizeMinusOne : 1> >
            HNormalizedReturnType;

    EIGEN_DEVICE_FUNC
    const HNormalizedReturnType hnormalized() const;

#   ifdef EIGEN_VECTORWISEOP_PLUGIN
#     include EIGEN_VECTORWISEOP_PLUGIN
#   endif

  protected:
    Index redux_length() const
    {
      return Direction==Vertical ? m_matrix.rows() : m_matrix.cols();
    }
    ExpressionTypeNested m_matrix;
};

template<typename Derived>
EIGEN_DEVICE_FUNC inline typename DenseBase<Derived>::ColwiseReturnType
DenseBase<Derived>::colwise()
{
  return ColwiseReturnType(derived());
}

template<typename Derived>
EIGEN_DEVICE_FUNC inline typename DenseBase<Derived>::RowwiseReturnType
DenseBase<Derived>::rowwise()
{
  return RowwiseReturnType(derived());
}

} 

#endif 
