#ifndef EIGEN_SPARSESPARSEPRODUCTWITHPRUNING_H
#define EIGEN_SPARSESPARSEPRODUCTWITHPRUNING_H

namespace Eigen { 

namespace internal {

template<typename Lhs, typename Rhs, typename ResultType>
static void sparse_sparse_product_with_pruning_impl(const Lhs& lhs, const Rhs& rhs, ResultType& res, const typename ResultType::RealScalar& tolerance)
{

  typedef typename remove_all<Rhs>::type::Scalar RhsScalar;
  typedef typename remove_all<ResultType>::type::Scalar ResScalar;
  typedef typename remove_all<Lhs>::type::StorageIndex StorageIndex;

  Index rows = lhs.innerSize();
  Index cols = rhs.outerSize();
  
  eigen_assert(lhs.outerSize() == rhs.innerSize());

  AmbiVector<ResScalar,StorageIndex> tempVector(rows);

  if(ResultType::IsRowMajor)
    res.resize(cols, rows);
  else
    res.resize(rows, cols);
  
  evaluator<Lhs> lhsEval(lhs);
  evaluator<Rhs> rhsEval(rhs);

  Index estimated_nnz_prod = lhsEval.nonZerosEstimate() + rhsEval.nonZerosEstimate();

  res.reserve(estimated_nnz_prod);
  double ratioColRes = double(estimated_nnz_prod)/(double(lhs.rows())*double(rhs.cols()));
  for (Index j=0; j<cols; ++j)
  {

    tempVector.init(ratioColRes);
    tempVector.setZero();
    for (typename evaluator<Rhs>::InnerIterator rhsIt(rhsEval, j); rhsIt; ++rhsIt)
    {
      
      tempVector.restart();
      RhsScalar x = rhsIt.value();
      for (typename evaluator<Lhs>::InnerIterator lhsIt(lhsEval, rhsIt.index()); lhsIt; ++lhsIt)
      {
        tempVector.coeffRef(lhsIt.index()) += lhsIt.value() * x;
      }
    }
    res.startVec(j);
    for (typename AmbiVector<ResScalar,StorageIndex>::Iterator it(tempVector,tolerance); it; ++it)
      res.insertBackByOuterInner(j,it.index()) = it.value();
  }
  res.finalize();
}

template<typename Lhs, typename Rhs, typename ResultType,
  int LhsStorageOrder = traits<Lhs>::Flags&RowMajorBit,
  int RhsStorageOrder = traits<Rhs>::Flags&RowMajorBit,
  int ResStorageOrder = traits<ResultType>::Flags&RowMajorBit>
struct sparse_sparse_product_with_pruning_selector;

template<typename Lhs, typename Rhs, typename ResultType>
struct sparse_sparse_product_with_pruning_selector<Lhs,Rhs,ResultType,ColMajor,ColMajor,ColMajor>
{
  typedef typename ResultType::RealScalar RealScalar;

  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res, const RealScalar& tolerance)
  {
    typename remove_all<ResultType>::type _res(res.rows(), res.cols());
    internal::sparse_sparse_product_with_pruning_impl<Lhs,Rhs,ResultType>(lhs, rhs, _res, tolerance);
    res.swap(_res);
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct sparse_sparse_product_with_pruning_selector<Lhs,Rhs,ResultType,ColMajor,ColMajor,RowMajor>
{
  typedef typename ResultType::RealScalar RealScalar;
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res, const RealScalar& tolerance)
  {
    
    typedef SparseMatrix<typename ResultType::Scalar,ColMajor,typename ResultType::StorageIndex> SparseTemporaryType;
    SparseTemporaryType _res(res.rows(), res.cols());
    internal::sparse_sparse_product_with_pruning_impl<Lhs,Rhs,SparseTemporaryType>(lhs, rhs, _res, tolerance);
    res = _res;
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct sparse_sparse_product_with_pruning_selector<Lhs,Rhs,ResultType,RowMajor,RowMajor,RowMajor>
{
  typedef typename ResultType::RealScalar RealScalar;
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res, const RealScalar& tolerance)
  {
    
    typename remove_all<ResultType>::type _res(res.rows(), res.cols());
    internal::sparse_sparse_product_with_pruning_impl<Rhs,Lhs,ResultType>(rhs, lhs, _res, tolerance);
    res.swap(_res);
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct sparse_sparse_product_with_pruning_selector<Lhs,Rhs,ResultType,RowMajor,RowMajor,ColMajor>
{
  typedef typename ResultType::RealScalar RealScalar;
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res, const RealScalar& tolerance)
  {
    typedef SparseMatrix<typename Lhs::Scalar,ColMajor,typename Lhs::StorageIndex> ColMajorMatrixLhs;
    typedef SparseMatrix<typename Rhs::Scalar,ColMajor,typename Lhs::StorageIndex> ColMajorMatrixRhs;
    ColMajorMatrixLhs colLhs(lhs);
    ColMajorMatrixRhs colRhs(rhs);
    internal::sparse_sparse_product_with_pruning_impl<ColMajorMatrixLhs,ColMajorMatrixRhs,ResultType>(colLhs, colRhs, res, tolerance);

  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct sparse_sparse_product_with_pruning_selector<Lhs,Rhs,ResultType,ColMajor,RowMajor,RowMajor>
{
  typedef typename ResultType::RealScalar RealScalar;
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res, const RealScalar& tolerance)
  {
    typedef SparseMatrix<typename Lhs::Scalar,RowMajor,typename Lhs::StorageIndex> RowMajorMatrixLhs;
    RowMajorMatrixLhs rowLhs(lhs);
    sparse_sparse_product_with_pruning_selector<RowMajorMatrixLhs,Rhs,ResultType,RowMajor,RowMajor>(rowLhs,rhs,res,tolerance);
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct sparse_sparse_product_with_pruning_selector<Lhs,Rhs,ResultType,RowMajor,ColMajor,RowMajor>
{
  typedef typename ResultType::RealScalar RealScalar;
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res, const RealScalar& tolerance)
  {
    typedef SparseMatrix<typename Rhs::Scalar,RowMajor,typename Lhs::StorageIndex> RowMajorMatrixRhs;
    RowMajorMatrixRhs rowRhs(rhs);
    sparse_sparse_product_with_pruning_selector<Lhs,RowMajorMatrixRhs,ResultType,RowMajor,RowMajor,RowMajor>(lhs,rowRhs,res,tolerance);
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct sparse_sparse_product_with_pruning_selector<Lhs,Rhs,ResultType,ColMajor,RowMajor,ColMajor>
{
  typedef typename ResultType::RealScalar RealScalar;
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res, const RealScalar& tolerance)
  {
    typedef SparseMatrix<typename Rhs::Scalar,ColMajor,typename Lhs::StorageIndex> ColMajorMatrixRhs;
    ColMajorMatrixRhs colRhs(rhs);
    internal::sparse_sparse_product_with_pruning_impl<Lhs,ColMajorMatrixRhs,ResultType>(lhs, colRhs, res, tolerance);
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct sparse_sparse_product_with_pruning_selector<Lhs,Rhs,ResultType,RowMajor,ColMajor,ColMajor>
{
  typedef typename ResultType::RealScalar RealScalar;
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res, const RealScalar& tolerance)
  {
    typedef SparseMatrix<typename Lhs::Scalar,ColMajor,typename Lhs::StorageIndex> ColMajorMatrixLhs;
    ColMajorMatrixLhs colLhs(lhs);
    internal::sparse_sparse_product_with_pruning_impl<ColMajorMatrixLhs,Rhs,ResultType>(colLhs, rhs, res, tolerance);
  }
};

} 

} 

#endif 
