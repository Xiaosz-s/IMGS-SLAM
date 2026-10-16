#ifndef EIGEN_SOLVE_H
#define EIGEN_SOLVE_H

namespace Eigen {

template<typename Decomposition, typename RhsType, typename StorageKind> class SolveImpl;

namespace internal {

template<typename Decomposition, typename RhsType,typename StorageKind> struct solve_traits;

template<typename Decomposition, typename RhsType>
struct solve_traits<Decomposition,RhsType,Dense>
{
  typedef typename make_proper_matrix_type<typename RhsType::Scalar,
                 Decomposition::ColsAtCompileTime,
                 RhsType::ColsAtCompileTime,
                 RhsType::PlainObject::Options,
                 Decomposition::MaxColsAtCompileTime,
                 RhsType::MaxColsAtCompileTime>::type PlainObject;
};

template<typename Decomposition, typename RhsType>
struct traits<Solve<Decomposition, RhsType> >
  : traits<typename solve_traits<Decomposition,RhsType,typename internal::traits<RhsType>::StorageKind>::PlainObject>
{
  typedef typename solve_traits<Decomposition,RhsType,typename internal::traits<RhsType>::StorageKind>::PlainObject PlainObject;
  typedef typename promote_index_type<typename Decomposition::StorageIndex, typename RhsType::StorageIndex>::type StorageIndex;
  typedef traits<PlainObject> BaseTraits;
  enum {
    Flags = BaseTraits::Flags & RowMajorBit,
    CoeffReadCost = HugeCost
  };
};

}

template<typename Decomposition, typename RhsType>
class Solve : public SolveImpl<Decomposition,RhsType,typename internal::traits<RhsType>::StorageKind>
{
public:
  typedef typename internal::traits<Solve>::PlainObject PlainObject;
  typedef typename internal::traits<Solve>::StorageIndex StorageIndex;

  Solve(const Decomposition &dec, const RhsType &rhs)
    : m_dec(dec), m_rhs(rhs)
  {}

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return m_dec.cols(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return m_rhs.cols(); }

  EIGEN_DEVICE_FUNC const Decomposition& dec() const { return m_dec; }
  EIGEN_DEVICE_FUNC const RhsType&       rhs() const { return m_rhs; }

protected:
  const Decomposition &m_dec;
  const typename internal::ref_selector<RhsType>::type m_rhs;
};

template<typename Decomposition, typename RhsType>
class SolveImpl<Decomposition,RhsType,Dense>
  : public MatrixBase<Solve<Decomposition,RhsType> >
{
  typedef Solve<Decomposition,RhsType> Derived;

public:

  typedef MatrixBase<Solve<Decomposition,RhsType> > Base;
  EIGEN_DENSE_PUBLIC_INTERFACE(Derived)

private:

  Scalar coeff(Index row, Index col) const;
  Scalar coeff(Index i) const;
};

template<typename Decomposition, typename RhsType, typename StorageKind>
class SolveImpl : public internal::generic_xpr_base<Solve<Decomposition,RhsType>, MatrixXpr, StorageKind>::type
{
  public:
    typedef typename internal::generic_xpr_base<Solve<Decomposition,RhsType>, MatrixXpr, StorageKind>::type Base;
};

namespace internal {

template<typename Decomposition, typename RhsType>
struct evaluator<Solve<Decomposition,RhsType> >
  : public evaluator<typename Solve<Decomposition,RhsType>::PlainObject>
{
  typedef Solve<Decomposition,RhsType> SolveType;
  typedef typename SolveType::PlainObject PlainObject;
  typedef evaluator<PlainObject> Base;

  enum { Flags = Base::Flags | EvalBeforeNestingBit };

  EIGEN_DEVICE_FUNC explicit evaluator(const SolveType& solve)
    : m_result(solve.rows(), solve.cols())
  {
    ::new (static_cast<Base*>(this)) Base(m_result);
    solve.dec()._solve_impl(solve.rhs(), m_result);
  }

protected:
  PlainObject m_result;
};

template<typename DstXprType, typename DecType, typename RhsType, typename Scalar>
struct Assignment<DstXprType, Solve<DecType,RhsType>, internal::assign_op<Scalar,Scalar>, Dense2Dense>
{
  typedef Solve<DecType,RhsType> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<Scalar,Scalar> &)
  {
    Index dstRows = src.rows();
    Index dstCols = src.cols();
    if((dst.rows()!=dstRows) || (dst.cols()!=dstCols))
      dst.resize(dstRows, dstCols);

    src.dec()._solve_impl(src.rhs(), dst);
  }
};

template<typename DstXprType, typename DecType, typename RhsType, typename Scalar>
struct Assignment<DstXprType, Solve<Transpose<const DecType>,RhsType>, internal::assign_op<Scalar,Scalar>, Dense2Dense>
{
  typedef Solve<Transpose<const DecType>,RhsType> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<Scalar,Scalar> &)
  {
    Index dstRows = src.rows();
    Index dstCols = src.cols();
    if((dst.rows()!=dstRows) || (dst.cols()!=dstCols))
      dst.resize(dstRows, dstCols);

    src.dec().nestedExpression().template _solve_impl_transposed<false>(src.rhs(), dst);
  }
};

template<typename DstXprType, typename DecType, typename RhsType, typename Scalar>
struct Assignment<DstXprType, Solve<CwiseUnaryOp<internal::scalar_conjugate_op<typename DecType::Scalar>, const Transpose<const DecType> >,RhsType>,
                  internal::assign_op<Scalar,Scalar>, Dense2Dense>
{
  typedef Solve<CwiseUnaryOp<internal::scalar_conjugate_op<typename DecType::Scalar>, const Transpose<const DecType> >,RhsType> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<Scalar,Scalar> &)
  {
    Index dstRows = src.rows();
    Index dstCols = src.cols();
    if((dst.rows()!=dstRows) || (dst.cols()!=dstCols))
      dst.resize(dstRows, dstCols);

    src.dec().nestedExpression().nestedExpression().template _solve_impl_transposed<true>(src.rhs(), dst);
  }
};

} 

} 

#endif 
