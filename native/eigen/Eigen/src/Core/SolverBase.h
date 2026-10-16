#ifndef EIGEN_SOLVERBASE_H
#define EIGEN_SOLVERBASE_H

namespace Eigen {

namespace internal {

template<typename Derived>
struct solve_assertion {
    template<bool Transpose_, typename Rhs>
    static void run(const Derived& solver, const Rhs& b) { solver.template _check_solve_assertion<Transpose_>(b); }
};

template<typename Derived>
struct solve_assertion<Transpose<Derived> >
{
    typedef Transpose<Derived> type;

    template<bool Transpose_, typename Rhs>
    static void run(const type& transpose, const Rhs& b)
    {
        internal::solve_assertion<typename internal::remove_all<Derived>::type>::template run<true>(transpose.nestedExpression(), b);
    }
};

template<typename Scalar, typename Derived>
struct solve_assertion<CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<Scalar>, const Transpose<Derived> > >
{
    typedef CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<Scalar>, const Transpose<Derived> > type;

    template<bool Transpose_, typename Rhs>
    static void run(const type& adjoint, const Rhs& b)
    {
        internal::solve_assertion<typename internal::remove_all<Transpose<Derived> >::type>::template run<true>(adjoint.nestedExpression(), b);
    }
};
} 

template<typename Derived>
class SolverBase : public EigenBase<Derived>
{
  public:

    typedef EigenBase<Derived> Base;
    typedef typename internal::traits<Derived>::Scalar Scalar;
    typedef Scalar CoeffReturnType;

    template<typename Derived_>
    friend struct internal::solve_assertion;

    enum {
      RowsAtCompileTime = internal::traits<Derived>::RowsAtCompileTime,
      ColsAtCompileTime = internal::traits<Derived>::ColsAtCompileTime,
      SizeAtCompileTime = (internal::size_at_compile_time<internal::traits<Derived>::RowsAtCompileTime,
                                                          internal::traits<Derived>::ColsAtCompileTime>::ret),
      MaxRowsAtCompileTime = internal::traits<Derived>::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = internal::traits<Derived>::MaxColsAtCompileTime,
      MaxSizeAtCompileTime = (internal::size_at_compile_time<internal::traits<Derived>::MaxRowsAtCompileTime,
                                                             internal::traits<Derived>::MaxColsAtCompileTime>::ret),
      IsVectorAtCompileTime = internal::traits<Derived>::MaxRowsAtCompileTime == 1
                           || internal::traits<Derived>::MaxColsAtCompileTime == 1,
      NumDimensions = int(MaxSizeAtCompileTime) == 1 ? 0 : bool(IsVectorAtCompileTime) ? 1 : 2
    };

    SolverBase()
    {}

    ~SolverBase()
    {}

    using Base::derived;

    template<typename Rhs>
    inline const Solve<Derived, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      internal::solve_assertion<typename internal::remove_all<Derived>::type>::template run<false>(derived(), b);
      return Solve<Derived, Rhs>(derived(), b.derived());
    }

    typedef typename internal::add_const<Transpose<const Derived> >::type ConstTransposeReturnType;

    inline ConstTransposeReturnType transpose() const
    {
      return ConstTransposeReturnType(derived());
    }

    typedef typename internal::conditional<NumTraits<Scalar>::IsComplex,
                        CwiseUnaryOp<internal::scalar_conjugate_op<Scalar>, ConstTransposeReturnType>,
                        ConstTransposeReturnType
                     >::type AdjointReturnType;

    inline AdjointReturnType adjoint() const
    {
      return AdjointReturnType(derived().transpose());
    }

  protected:

    template<bool Transpose_, typename Rhs>
    void _check_solve_assertion(const Rhs& b) const {
        EIGEN_ONLY_USED_FOR_DEBUG(b);
        eigen_assert(derived().m_isInitialized && "Solver is not initialized.");
        eigen_assert((Transpose_?derived().cols():derived().rows())==b.rows() && "SolverBase::solve(): invalid number of rows of the right hand side matrix b");
    }
};

namespace internal {

template<typename Derived>
struct generic_xpr_base<Derived, MatrixXpr, SolverStorage>
{
  typedef SolverBase<Derived> type;

};

} 

} 

#endif 
