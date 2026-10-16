#ifndef EIGEN_SPARSEMATRIXBASE_H
#define EIGEN_SPARSEMATRIXBASE_H

namespace Eigen { 

template<typename Derived> class SparseMatrixBase
  : public EigenBase<Derived>
{
  public:

    typedef typename internal::traits<Derived>::Scalar Scalar;

    typedef Scalar value_type;
    
    typedef typename internal::packet_traits<Scalar>::type PacketScalar;
    typedef typename internal::traits<Derived>::StorageKind StorageKind;

    typedef typename internal::traits<Derived>::StorageIndex StorageIndex;

    typedef typename internal::add_const_on_value_type_if_arithmetic<
                         typename internal::packet_traits<Scalar>::type
                     >::type PacketReturnType;

    typedef SparseMatrixBase StorageBaseType;

    typedef Matrix<StorageIndex,Dynamic,1> IndexVector;
    typedef Matrix<Scalar,Dynamic,1> ScalarVector;
    
    template<typename OtherDerived>
    Derived& operator=(const EigenBase<OtherDerived> &other);

    enum {

      RowsAtCompileTime = internal::traits<Derived>::RowsAtCompileTime,

      ColsAtCompileTime = internal::traits<Derived>::ColsAtCompileTime,

      SizeAtCompileTime = (internal::size_at_compile_time<internal::traits<Derived>::RowsAtCompileTime,
                                                   internal::traits<Derived>::ColsAtCompileTime>::ret),

      MaxRowsAtCompileTime = RowsAtCompileTime,
      MaxColsAtCompileTime = ColsAtCompileTime,

      MaxSizeAtCompileTime = (internal::size_at_compile_time<MaxRowsAtCompileTime,
                                                      MaxColsAtCompileTime>::ret),

      IsVectorAtCompileTime = RowsAtCompileTime == 1 || ColsAtCompileTime == 1,

      NumDimensions = int(MaxSizeAtCompileTime) == 1 ? 0 : bool(IsVectorAtCompileTime) ? 1 : 2,

      Flags = internal::traits<Derived>::Flags,

      IsRowMajor = Flags&RowMajorBit ? 1 : 0,
      
      InnerSizeAtCompileTime = int(IsVectorAtCompileTime) ? int(SizeAtCompileTime)
                             : int(IsRowMajor) ? int(ColsAtCompileTime) : int(RowsAtCompileTime),

      #ifndef EIGEN_PARSED_BY_DOXYGEN
      _HasDirectAccess = (int(Flags)&DirectAccessBit) ? 1 : 0 
      #endif
    };

    typedef typename internal::conditional<NumTraits<Scalar>::IsComplex,
                        CwiseUnaryOp<internal::scalar_conjugate_op<Scalar>, Eigen::Transpose<const Derived> >,
                        Transpose<const Derived>
                     >::type AdjointReturnType;
    typedef Transpose<Derived> TransposeReturnType;
    typedef typename internal::add_const<Transpose<const Derived> >::type ConstTransposeReturnType;

    typedef SparseMatrix<Scalar, Flags&RowMajorBit ? RowMajor : ColMajor, StorageIndex> PlainObject;

#ifndef EIGEN_PARSED_BY_DOXYGEN

    typedef typename NumTraits<Scalar>::Real RealScalar;

    typedef typename internal::conditional<_HasDirectAccess, const Scalar&, Scalar>::type CoeffReturnType;

    typedef CwiseNullaryOp<internal::scalar_constant_op<Scalar>,Matrix<Scalar,Dynamic,Dynamic> > ConstantReturnType;

    typedef Matrix<Scalar,RowsAtCompileTime,ColsAtCompileTime> DenseMatrixType;
    
    typedef Matrix<Scalar,EIGEN_SIZE_MAX(RowsAtCompileTime,ColsAtCompileTime),
                          EIGEN_SIZE_MAX(RowsAtCompileTime,ColsAtCompileTime)> SquareMatrixType;

    inline const Derived& derived() const { return *static_cast<const Derived*>(this); }
    inline Derived& derived() { return *static_cast<Derived*>(this); }
    inline Derived& const_cast_derived() const
    { return *static_cast<Derived*>(const_cast<SparseMatrixBase*>(this)); }

    typedef EigenBase<Derived> Base;

#endif 

#define EIGEN_CURRENT_STORAGE_BASE_CLASS Eigen::SparseMatrixBase
#ifdef EIGEN_PARSED_BY_DOXYGEN
#define EIGEN_DOC_UNARY_ADDONS(METHOD,OP)           
#define EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL      
#define EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(COND) 
#else
#define EIGEN_DOC_UNARY_ADDONS(X,Y)
#define EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL
#define EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF(COND)
#endif
#   include "../plugins/CommonCwiseUnaryOps.h"
#   include "../plugins/CommonCwiseBinaryOps.h"
#   include "../plugins/MatrixCwiseUnaryOps.h"
#   include "../plugins/MatrixCwiseBinaryOps.h"
#   include "../plugins/BlockMethods.h"
#   ifdef EIGEN_SPARSEMATRIXBASE_PLUGIN
#     include EIGEN_SPARSEMATRIXBASE_PLUGIN
#   endif
#undef EIGEN_CURRENT_STORAGE_BASE_CLASS
#undef EIGEN_DOC_UNARY_ADDONS
#undef EIGEN_DOC_BLOCK_ADDONS_NOT_INNER_PANEL
#undef EIGEN_DOC_BLOCK_ADDONS_INNER_PANEL_IF

    inline Index rows() const { return derived().rows(); }
    
    inline Index cols() const { return derived().cols(); }

    inline Index size() const { return rows() * cols(); }

    inline bool isVector() const { return rows()==1 || cols()==1; }

    Index outerSize() const { return (int(Flags)&RowMajorBit) ? this->rows() : this->cols(); }

    Index innerSize() const { return (int(Flags)&RowMajorBit) ? this->cols() : this->rows(); }

    bool isRValue() const { return m_isRValue; }
    Derived& markAsRValue() { m_isRValue = true; return derived(); }

    SparseMatrixBase() : m_isRValue(false) {  }

    template<typename OtherDerived>
    Derived& operator=(const ReturnByValue<OtherDerived>& other);

    template<typename OtherDerived>
    inline Derived& operator=(const SparseMatrixBase<OtherDerived>& other);

    inline Derived& operator=(const Derived& other);

  protected:

    template<typename OtherDerived>
    inline Derived& assign(const OtherDerived& other);

    template<typename OtherDerived>
    inline void assignGeneric(const OtherDerived& other);

  public:

    friend std::ostream & operator << (std::ostream & s, const SparseMatrixBase& m)
    {
      typedef typename Derived::Nested Nested;
      typedef typename internal::remove_all<Nested>::type NestedCleaned;

      if (Flags&RowMajorBit)
      {
        Nested nm(m.derived());
        internal::evaluator<NestedCleaned> thisEval(nm);
        for (Index row=0; row<nm.outerSize(); ++row)
        {
          Index col = 0;
          for (typename internal::evaluator<NestedCleaned>::InnerIterator it(thisEval, row); it; ++it)
          {
            for ( ; col<it.index(); ++col)
              s << "0 ";
            s << it.value() << " ";
            ++col;
          }
          for ( ; col<m.cols(); ++col)
            s << "0 ";
          s << std::endl;
        }
      }
      else
      {
        Nested nm(m.derived());
        internal::evaluator<NestedCleaned> thisEval(nm);
        if (m.cols() == 1) {
          Index row = 0;
          for (typename internal::evaluator<NestedCleaned>::InnerIterator it(thisEval, 0); it; ++it)
          {
            for ( ; row<it.index(); ++row)
              s << "0" << std::endl;
            s << it.value() << std::endl;
            ++row;
          }
          for ( ; row<m.rows(); ++row)
            s << "0" << std::endl;
        }
        else
        {
          SparseMatrix<Scalar, RowMajorBit, StorageIndex> trans = m;
          s << static_cast<const SparseMatrixBase<SparseMatrix<Scalar, RowMajorBit, StorageIndex> >&>(trans);
        }
      }
      return s;
    }

    template<typename OtherDerived>
    Derived& operator+=(const SparseMatrixBase<OtherDerived>& other);
    template<typename OtherDerived>
    Derived& operator-=(const SparseMatrixBase<OtherDerived>& other);
    
    template<typename OtherDerived>
    Derived& operator+=(const DiagonalBase<OtherDerived>& other);
    template<typename OtherDerived>
    Derived& operator-=(const DiagonalBase<OtherDerived>& other);

    template<typename OtherDerived>
    Derived& operator+=(const EigenBase<OtherDerived> &other);
    template<typename OtherDerived>
    Derived& operator-=(const EigenBase<OtherDerived> &other);

    Derived& operator*=(const Scalar& other);
    Derived& operator/=(const Scalar& other);

    template<typename OtherDerived> struct CwiseProductDenseReturnType {
      typedef CwiseBinaryOp<internal::scalar_product_op<typename ScalarBinaryOpTraits<
                                                          typename internal::traits<Derived>::Scalar,
                                                          typename internal::traits<OtherDerived>::Scalar
                                                        >::ReturnType>,
                            const Derived,
                            const OtherDerived
                          > Type;
    };

    template<typename OtherDerived>
    EIGEN_STRONG_INLINE const typename CwiseProductDenseReturnType<OtherDerived>::Type
    cwiseProduct(const MatrixBase<OtherDerived> &other) const;

    template<typename OtherDerived>
    const Product<Derived,OtherDerived>
    operator*(const DiagonalBase<OtherDerived> &other) const
    { return Product<Derived,OtherDerived>(derived(), other.derived()); }

    template<typename OtherDerived> friend
    const Product<OtherDerived,Derived>
    operator*(const DiagonalBase<OtherDerived> &lhs, const SparseMatrixBase& rhs)
    { return Product<OtherDerived,Derived>(lhs.derived(), rhs.derived()); }

    template<typename OtherDerived>
    const Product<Derived,OtherDerived,AliasFreeProduct>
    operator*(const SparseMatrixBase<OtherDerived> &other) const;

    template<typename OtherDerived>
    const Product<Derived,OtherDerived>
    operator*(const MatrixBase<OtherDerived> &other) const
    { return Product<Derived,OtherDerived>(derived(), other.derived()); }

    template<typename OtherDerived> friend
    const Product<OtherDerived,Derived>
    operator*(const MatrixBase<OtherDerived> &lhs, const SparseMatrixBase& rhs)
    { return Product<OtherDerived,Derived>(lhs.derived(), rhs.derived()); }

    SparseSymmetricPermutationProduct<Derived,Upper|Lower> twistedBy(const PermutationMatrix<Dynamic,Dynamic,StorageIndex>& perm) const
    {
      return SparseSymmetricPermutationProduct<Derived,Upper|Lower>(derived(), perm);
    }

    template<typename OtherDerived>
    Derived& operator*=(const SparseMatrixBase<OtherDerived>& other);

    template<int Mode>
    inline const TriangularView<const Derived, Mode> triangularView() const;
    
    template<unsigned int UpLo> struct SelfAdjointViewReturnType { typedef SparseSelfAdjointView<Derived, UpLo> Type; };
    template<unsigned int UpLo> struct ConstSelfAdjointViewReturnType { typedef const SparseSelfAdjointView<const Derived, UpLo> Type; };

    template<unsigned int UpLo> inline 
    typename ConstSelfAdjointViewReturnType<UpLo>::Type selfadjointView() const;
    template<unsigned int UpLo> inline
    typename SelfAdjointViewReturnType<UpLo>::Type selfadjointView();

    template<typename OtherDerived> Scalar dot(const MatrixBase<OtherDerived>& other) const;
    template<typename OtherDerived> Scalar dot(const SparseMatrixBase<OtherDerived>& other) const;
    RealScalar squaredNorm() const;
    RealScalar norm()  const;
    RealScalar blueNorm() const;

    TransposeReturnType transpose() { return TransposeReturnType(derived()); }
    const ConstTransposeReturnType transpose() const { return ConstTransposeReturnType(derived()); }
    const AdjointReturnType adjoint() const { return AdjointReturnType(transpose()); }

    DenseMatrixType toDense() const
    {
      return DenseMatrixType(derived());
    }

    template<typename OtherDerived>
    bool isApprox(const SparseMatrixBase<OtherDerived>& other,
                  const RealScalar& prec = NumTraits<Scalar>::dummy_precision()) const;

    template<typename OtherDerived>
    bool isApprox(const MatrixBase<OtherDerived>& other,
                  const RealScalar& prec = NumTraits<Scalar>::dummy_precision()) const
    { return toDense().isApprox(other,prec); }

    inline const typename internal::eval<Derived>::type eval() const
    { return typename internal::eval<Derived>::type(derived()); }

    Scalar sum() const;
    
    inline const SparseView<Derived>
    pruned(const Scalar& reference = Scalar(0), const RealScalar& epsilon = NumTraits<Scalar>::dummy_precision()) const;

  protected:

    bool m_isRValue;

    static inline StorageIndex convert_index(const Index idx) {
      return internal::convert_index<StorageIndex>(idx);
    }
  private:
    template<typename Dest> void evalTo(Dest &) const;
};

} 

#endif 
