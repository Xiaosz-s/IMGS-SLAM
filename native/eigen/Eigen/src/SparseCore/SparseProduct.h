#ifndef EIGEN_SPARSEPRODUCT_H
#define EIGEN_SPARSEPRODUCT_H

namespace Eigen { 

template<typename Derived>
template<typename OtherDerived>
inline const Product<Derived,OtherDerived,AliasFreeProduct>
SparseMatrixBase<Derived>::operator*(const SparseMatrixBase<OtherDerived> &other) const
{
  return Product<Derived,OtherDerived,AliasFreeProduct>(derived(), other.derived());
}

namespace internal {

template<typename Lhs, typename Rhs, int ProductType>
struct generic_product_impl<Lhs, Rhs, SparseShape, SparseShape, ProductType>
{
  template<typename Dest>
  static void evalTo(Dest& dst, const Lhs& lhs, const Rhs& rhs)
  {
    evalTo(dst, lhs, rhs, typename evaluator_traits<Dest>::Shape());
  }

  template<typename Dest,typename ActualLhs>
  static void addTo(Dest& dst, const ActualLhs& lhs, const Rhs& rhs, typename enable_if<is_same<typename evaluator_traits<Dest>::Shape,DenseShape>::value,int*>::type* = 0)
  {
    typedef typename nested_eval<ActualLhs,Dynamic>::type LhsNested;
    typedef typename nested_eval<Rhs,Dynamic>::type RhsNested;
    LhsNested lhsNested(lhs);
    RhsNested rhsNested(rhs);
    internal::sparse_sparse_to_dense_product_selector<typename remove_all<LhsNested>::type,
                                                      typename remove_all<RhsNested>::type, Dest>::run(lhsNested,rhsNested,dst);
  }

  template<typename Dest>
  static void subTo(Dest& dst, const Lhs& lhs, const Rhs& rhs, typename enable_if<is_same<typename evaluator_traits<Dest>::Shape,DenseShape>::value,int*>::type* = 0)
  {
    addTo(dst, -lhs, rhs);
  }

protected:

  template<typename Dest>
  static void evalTo(Dest& dst, const Lhs& lhs, const Rhs& rhs, SparseShape)
  {
    typedef typename nested_eval<Lhs,Dynamic>::type LhsNested;
    typedef typename nested_eval<Rhs,Dynamic>::type RhsNested;
    LhsNested lhsNested(lhs);
    RhsNested rhsNested(rhs);
    internal::conservative_sparse_sparse_product_selector<typename remove_all<LhsNested>::type,
                                                          typename remove_all<RhsNested>::type, Dest>::run(lhsNested,rhsNested,dst);
  }

  template<typename Dest>
  static void evalTo(Dest& dst, const Lhs& lhs, const Rhs& rhs, DenseShape)
  {
    dst.setZero();
    addTo(dst, lhs, rhs);
  }
};

template<typename Lhs, typename Rhs, int ProductType>
struct generic_product_impl<Lhs, Rhs, SparseShape, SparseTriangularShape, ProductType>
 : public generic_product_impl<Lhs, Rhs, SparseShape, SparseShape, ProductType>
{};

template<typename Lhs, typename Rhs, int ProductType>
struct generic_product_impl<Lhs, Rhs, SparseTriangularShape, SparseShape, ProductType>
 : public generic_product_impl<Lhs, Rhs, SparseShape, SparseShape, ProductType>
{};

template< typename DstXprType, typename Lhs, typename Rhs>
struct Assignment<DstXprType, Product<Lhs,Rhs,AliasFreeProduct>, internal::assign_op<typename DstXprType::Scalar,typename Product<Lhs,Rhs,AliasFreeProduct>::Scalar>, Sparse2Dense>
{
  typedef Product<Lhs,Rhs,AliasFreeProduct> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<typename DstXprType::Scalar,typename SrcXprType::Scalar> &)
  {
    Index dstRows = src.rows();
    Index dstCols = src.cols();
    if((dst.rows()!=dstRows) || (dst.cols()!=dstCols))
      dst.resize(dstRows, dstCols);
    
    generic_product_impl<Lhs, Rhs>::evalTo(dst,src.lhs(),src.rhs());
  }
};

template< typename DstXprType, typename Lhs, typename Rhs>
struct Assignment<DstXprType, Product<Lhs,Rhs,AliasFreeProduct>, internal::add_assign_op<typename DstXprType::Scalar,typename Product<Lhs,Rhs,AliasFreeProduct>::Scalar>, Sparse2Dense>
{
  typedef Product<Lhs,Rhs,AliasFreeProduct> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::add_assign_op<typename DstXprType::Scalar,typename SrcXprType::Scalar> &)
  {
    generic_product_impl<Lhs, Rhs>::addTo(dst,src.lhs(),src.rhs());
  }
};

template< typename DstXprType, typename Lhs, typename Rhs>
struct Assignment<DstXprType, Product<Lhs,Rhs,AliasFreeProduct>, internal::sub_assign_op<typename DstXprType::Scalar,typename Product<Lhs,Rhs,AliasFreeProduct>::Scalar>, Sparse2Dense>
{
  typedef Product<Lhs,Rhs,AliasFreeProduct> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::sub_assign_op<typename DstXprType::Scalar,typename SrcXprType::Scalar> &)
  {
    generic_product_impl<Lhs, Rhs>::subTo(dst,src.lhs(),src.rhs());
  }
};

template<typename Lhs, typename Rhs, int Options>
struct unary_evaluator<SparseView<Product<Lhs, Rhs, Options> >, IteratorBased>
 : public evaluator<typename Product<Lhs, Rhs, DefaultProduct>::PlainObject>
{
  typedef SparseView<Product<Lhs, Rhs, Options> > XprType;
  typedef typename XprType::PlainObject PlainObject;
  typedef evaluator<PlainObject> Base;

  explicit unary_evaluator(const XprType& xpr)
    : m_result(xpr.rows(), xpr.cols())
  {
    using std::abs;
    ::new (static_cast<Base*>(this)) Base(m_result);
    typedef typename nested_eval<Lhs,Dynamic>::type LhsNested;
    typedef typename nested_eval<Rhs,Dynamic>::type RhsNested;
    LhsNested lhsNested(xpr.nestedExpression().lhs());
    RhsNested rhsNested(xpr.nestedExpression().rhs());

    internal::sparse_sparse_product_with_pruning_selector<typename remove_all<LhsNested>::type,
                                                          typename remove_all<RhsNested>::type, PlainObject>::run(lhsNested,rhsNested,m_result,
                                                                                                                  abs(xpr.reference())*xpr.epsilon());
  }

protected:
  PlainObject m_result;
};

} 

template<typename Scalar, int Options_, typename StorageIndex_>
template<typename Lhs, typename Rhs>
SparseMatrix<Scalar,Options_,StorageIndex_>& SparseMatrix<Scalar,Options_,StorageIndex_>::operator=(const Product<Lhs,Rhs,AliasFreeProduct>& src)
{
  
  SparseMatrix dst(src.rows(),src.cols());
  internal::generic_product_impl<Lhs, Rhs>::evalTo(dst,src.lhs(),src.rhs());
  this->swap(dst);
  return *this;
}

} 

#endif 
