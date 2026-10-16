#ifndef EIGEN_SELFADJOINTMATRIX_H
#define EIGEN_SELFADJOINTMATRIX_H

namespace Eigen {

namespace internal {
template<typename MatrixType, unsigned int UpLo>
struct traits<SelfAdjointView<MatrixType, UpLo> > : traits<MatrixType>
{
  typedef typename ref_selector<MatrixType>::non_const_type MatrixTypeNested;
  typedef typename remove_all<MatrixTypeNested>::type MatrixTypeNestedCleaned;
  typedef MatrixType ExpressionType;
  typedef typename MatrixType::PlainObject FullMatrixType;
  enum {
    Mode = UpLo | SelfAdjoint,
    FlagsLvalueBit = is_lvalue<MatrixType>::value ? LvalueBit : 0,
    Flags =  MatrixTypeNestedCleaned::Flags & (HereditaryBits|FlagsLvalueBit)
           & (~(PacketAccessBit | DirectAccessBit | LinearAccessBit)) 
  };
};
}

template<typename MatrixType_, unsigned int UpLo> class SelfAdjointView
  : public TriangularBase<SelfAdjointView<MatrixType_, UpLo> >
{
  public:

    typedef MatrixType_ MatrixType;
    typedef TriangularBase<SelfAdjointView> Base;
    typedef typename internal::traits<SelfAdjointView>::MatrixTypeNested MatrixTypeNested;
    typedef typename internal::traits<SelfAdjointView>::MatrixTypeNestedCleaned MatrixTypeNestedCleaned;
    typedef MatrixTypeNestedCleaned NestedExpression;

    typedef typename internal::traits<SelfAdjointView>::Scalar Scalar;
    typedef typename MatrixType::StorageIndex StorageIndex;
    typedef typename internal::remove_all<typename MatrixType::ConjugateReturnType>::type MatrixConjugateReturnType;
    typedef SelfAdjointView<typename internal::add_const<MatrixType>::type, UpLo> ConstSelfAdjointView;

    enum {
      Mode = internal::traits<SelfAdjointView>::Mode,
      Flags = internal::traits<SelfAdjointView>::Flags,
      TransposeMode = ((int(Mode) & int(Upper)) ? Lower : 0) | ((int(Mode) & int(Lower)) ? Upper : 0)
    };
    typedef typename MatrixType::PlainObject PlainObject;

    EIGEN_DEVICE_FUNC
    explicit inline SelfAdjointView(MatrixType& matrix) : m_matrix(matrix)
    {
      EIGEN_STATIC_ASSERT(UpLo==Lower || UpLo==Upper,SELFADJOINTVIEW_ACCEPTS_UPPER_AND_LOWER_MODE_ONLY);
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index rows() const EIGEN_NOEXCEPT { return m_matrix.rows(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index cols() const EIGEN_NOEXCEPT { return m_matrix.cols(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index outerStride() const EIGEN_NOEXCEPT { return m_matrix.outerStride(); }
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index innerStride() const EIGEN_NOEXCEPT { return m_matrix.innerStride(); }

    EIGEN_DEVICE_FUNC
    inline Scalar coeff(Index row, Index col) const
    {
      Base::check_coordinates_internal(row, col);
      return m_matrix.coeff(row, col);
    }

    EIGEN_DEVICE_FUNC
    inline Scalar& coeffRef(Index row, Index col)
    {
      EIGEN_STATIC_ASSERT_LVALUE(SelfAdjointView);
      Base::check_coordinates_internal(row, col);
      return m_matrix.coeffRef(row, col);
    }

    EIGEN_DEVICE_FUNC
    const MatrixTypeNestedCleaned& _expression() const { return m_matrix; }

    EIGEN_DEVICE_FUNC
    const MatrixTypeNestedCleaned& nestedExpression() const { return m_matrix; }
    EIGEN_DEVICE_FUNC
    MatrixTypeNestedCleaned& nestedExpression() { return m_matrix; }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    const Product<SelfAdjointView,OtherDerived>
    operator*(const MatrixBase<OtherDerived>& rhs) const
    {
      return Product<SelfAdjointView,OtherDerived>(*this, rhs.derived());
    }

    template<typename OtherDerived> friend
    EIGEN_DEVICE_FUNC
    const Product<OtherDerived,SelfAdjointView>
    operator*(const MatrixBase<OtherDerived>& lhs, const SelfAdjointView& rhs)
    {
      return Product<OtherDerived,SelfAdjointView>(lhs.derived(),rhs);
    }

    friend EIGEN_DEVICE_FUNC
    const SelfAdjointView<const EIGEN_SCALAR_BINARYOP_EXPR_RETURN_TYPE(Scalar,MatrixType,product),UpLo>
    operator*(const Scalar& s, const SelfAdjointView& mat)
    {
      return (s*mat.nestedExpression()).template selfadjointView<UpLo>();
    }

    template<typename DerivedU, typename DerivedV>
    EIGEN_DEVICE_FUNC
    SelfAdjointView& rankUpdate(const MatrixBase<DerivedU>& u, const MatrixBase<DerivedV>& v, const Scalar& alpha = Scalar(1));

    template<typename DerivedU>
    EIGEN_DEVICE_FUNC
    SelfAdjointView& rankUpdate(const MatrixBase<DerivedU>& u, const Scalar& alpha = Scalar(1));

    template<unsigned int TriMode>
    EIGEN_DEVICE_FUNC
    typename internal::conditional<(TriMode&(Upper|Lower))==(UpLo&(Upper|Lower)),
                                   TriangularView<MatrixType,TriMode>,
                                   TriangularView<typename MatrixType::AdjointReturnType,TriMode> >::type
    triangularView() const
    {
      typename internal::conditional<(TriMode&(Upper|Lower))==(UpLo&(Upper|Lower)), MatrixType&, typename MatrixType::ConstTransposeReturnType>::type tmp1(m_matrix);
      typename internal::conditional<(TriMode&(Upper|Lower))==(UpLo&(Upper|Lower)), MatrixType&, typename MatrixType::AdjointReturnType>::type tmp2(tmp1);
      return typename internal::conditional<(TriMode&(Upper|Lower))==(UpLo&(Upper|Lower)),
                                   TriangularView<MatrixType,TriMode>,
                                   TriangularView<typename MatrixType::AdjointReturnType,TriMode> >::type(tmp2);
    }

    typedef SelfAdjointView<const MatrixConjugateReturnType,UpLo> ConjugateReturnType;
    
    EIGEN_DEVICE_FUNC
    inline const ConjugateReturnType conjugate() const
    { return ConjugateReturnType(m_matrix.conjugate()); }

    template<bool Cond>
    EIGEN_DEVICE_FUNC
    inline typename internal::conditional<Cond,ConjugateReturnType,ConstSelfAdjointView>::type
    conjugateIf() const
    {
      typedef typename internal::conditional<Cond,ConjugateReturnType,ConstSelfAdjointView>::type ReturnType;
      return ReturnType(m_matrix.template conjugateIf<Cond>());
    }

    typedef SelfAdjointView<const typename MatrixType::AdjointReturnType,TransposeMode> AdjointReturnType;
    
    EIGEN_DEVICE_FUNC
    inline const AdjointReturnType adjoint() const
    { return AdjointReturnType(m_matrix.adjoint()); }

    typedef SelfAdjointView<typename MatrixType::TransposeReturnType,TransposeMode> TransposeReturnType;
     
    EIGEN_DEVICE_FUNC
    inline TransposeReturnType transpose()
    {
      EIGEN_STATIC_ASSERT_LVALUE(MatrixType)
      typename MatrixType::TransposeReturnType tmp(m_matrix);
      return TransposeReturnType(tmp);
    }

    typedef SelfAdjointView<const typename MatrixType::ConstTransposeReturnType,TransposeMode> ConstTransposeReturnType;
    
    EIGEN_DEVICE_FUNC
    inline const ConstTransposeReturnType transpose() const
    {
      return ConstTransposeReturnType(m_matrix.transpose());
    }

    EIGEN_DEVICE_FUNC
    typename MatrixType::ConstDiagonalReturnType diagonal() const
    {
      return typename MatrixType::ConstDiagonalReturnType(m_matrix);
    }

    const LLT<PlainObject, UpLo> llt() const;
    const LDLT<PlainObject, UpLo> ldlt() const;

    typedef typename NumTraits<Scalar>::Real RealScalar;
    
    typedef Matrix<RealScalar, internal::traits<MatrixType>::ColsAtCompileTime, 1> EigenvaluesReturnType;

    EIGEN_DEVICE_FUNC
    EigenvaluesReturnType eigenvalues() const;
    EIGEN_DEVICE_FUNC
    RealScalar operatorNorm() const;

  protected:
    MatrixTypeNested m_matrix;
};

namespace internal {

template<typename MatrixType, unsigned int Mode>
struct evaluator_traits<SelfAdjointView<MatrixType,Mode> >
{
  typedef typename storage_kind_to_evaluator_kind<typename MatrixType::StorageKind>::Kind Kind;
  typedef SelfAdjointShape Shape;
};

template<int UpLo, int SetOpposite, typename DstEvaluatorTypeT, typename SrcEvaluatorTypeT, typename Functor, int Version>
class triangular_dense_assignment_kernel<UpLo,SelfAdjoint,SetOpposite,DstEvaluatorTypeT,SrcEvaluatorTypeT,Functor,Version>
  : public generic_dense_assignment_kernel<DstEvaluatorTypeT, SrcEvaluatorTypeT, Functor, Version>
{
protected:
  typedef generic_dense_assignment_kernel<DstEvaluatorTypeT, SrcEvaluatorTypeT, Functor, Version> Base;
  typedef typename Base::DstXprType DstXprType;
  typedef typename Base::SrcXprType SrcXprType;
  using Base::m_dst;
  using Base::m_src;
  using Base::m_functor;
public:

  typedef typename Base::DstEvaluatorType DstEvaluatorType;
  typedef typename Base::SrcEvaluatorType SrcEvaluatorType;
  typedef typename Base::Scalar Scalar;
  typedef typename Base::AssignmentTraits AssignmentTraits;

  EIGEN_DEVICE_FUNC triangular_dense_assignment_kernel(DstEvaluatorType &dst, const SrcEvaluatorType &src, const Functor &func, DstXprType& dstExpr)
    : Base(dst, src, func, dstExpr)
  {}

  EIGEN_DEVICE_FUNC void assignCoeff(Index row, Index col)
  {
    eigen_internal_assert(row!=col);
    Scalar tmp = m_src.coeff(row,col);
    m_functor.assignCoeff(m_dst.coeffRef(row,col), tmp);
    m_functor.assignCoeff(m_dst.coeffRef(col,row), numext::conj(tmp));
  }

  EIGEN_DEVICE_FUNC void assignDiagonalCoeff(Index id)
  {
    Base::assignCoeff(id,id);
  }

  EIGEN_DEVICE_FUNC void assignOppositeCoeff(Index, Index)
  { eigen_internal_assert(false && "should never be called"); }
};

} 

template<typename Derived>
template<unsigned int UpLo>
EIGEN_DEVICE_FUNC typename MatrixBase<Derived>::template ConstSelfAdjointViewReturnType<UpLo>::Type
MatrixBase<Derived>::selfadjointView() const
{
  return typename ConstSelfAdjointViewReturnType<UpLo>::Type(derived());
}

template<typename Derived>
template<unsigned int UpLo>
EIGEN_DEVICE_FUNC typename MatrixBase<Derived>::template SelfAdjointViewReturnType<UpLo>::Type
MatrixBase<Derived>::selfadjointView()
{
  return typename SelfAdjointViewReturnType<UpLo>::Type(derived());
}

} 

#endif 
