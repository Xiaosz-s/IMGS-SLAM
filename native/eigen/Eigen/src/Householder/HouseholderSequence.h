#ifndef EIGEN_HOUSEHOLDER_SEQUENCE_H
#define EIGEN_HOUSEHOLDER_SEQUENCE_H

namespace Eigen {

namespace internal {

template<typename VectorsType, typename CoeffsType, int Side>
struct traits<HouseholderSequence<VectorsType,CoeffsType,Side> >
{
  typedef typename VectorsType::Scalar Scalar;
  typedef typename VectorsType::StorageIndex StorageIndex;
  typedef typename VectorsType::StorageKind StorageKind;
  enum {
    RowsAtCompileTime = Side==OnTheLeft ? traits<VectorsType>::RowsAtCompileTime
                                        : traits<VectorsType>::ColsAtCompileTime,
    ColsAtCompileTime = RowsAtCompileTime,
    MaxRowsAtCompileTime = Side==OnTheLeft ? traits<VectorsType>::MaxRowsAtCompileTime
                                           : traits<VectorsType>::MaxColsAtCompileTime,
    MaxColsAtCompileTime = MaxRowsAtCompileTime,
    Flags = 0
  };
};

struct HouseholderSequenceShape {};

template<typename VectorsType, typename CoeffsType, int Side>
struct evaluator_traits<HouseholderSequence<VectorsType,CoeffsType,Side> >
  : public evaluator_traits_base<HouseholderSequence<VectorsType,CoeffsType,Side> >
{
  typedef HouseholderSequenceShape Shape;
};

template<typename VectorsType, typename CoeffsType, int Side>
struct hseq_side_dependent_impl
{
  typedef Block<const VectorsType, Dynamic, 1> EssentialVectorType;
  typedef HouseholderSequence<VectorsType, CoeffsType, OnTheLeft> HouseholderSequenceType;
  static EIGEN_DEVICE_FUNC inline const EssentialVectorType essentialVector(const HouseholderSequenceType& h, Index k)
  {
    Index start = k+1+h.m_shift;
    return Block<const VectorsType,Dynamic,1>(h.m_vectors, start, k, h.rows()-start, 1);
  }
};

template<typename VectorsType, typename CoeffsType>
struct hseq_side_dependent_impl<VectorsType, CoeffsType, OnTheRight>
{
  typedef Transpose<Block<const VectorsType, 1, Dynamic> > EssentialVectorType;
  typedef HouseholderSequence<VectorsType, CoeffsType, OnTheRight> HouseholderSequenceType;
  static inline const EssentialVectorType essentialVector(const HouseholderSequenceType& h, Index k)
  {
    Index start = k+1+h.m_shift;
    return Block<const VectorsType,1,Dynamic>(h.m_vectors, k, start, 1, h.rows()-start).transpose();
  }
};

template<typename OtherScalarType, typename MatrixType> struct matrix_type_times_scalar_type
{
  typedef typename ScalarBinaryOpTraits<OtherScalarType, typename MatrixType::Scalar>::ReturnType
    ResultScalar;
  typedef Matrix<ResultScalar, MatrixType::RowsAtCompileTime, MatrixType::ColsAtCompileTime,
                 0, MatrixType::MaxRowsAtCompileTime, MatrixType::MaxColsAtCompileTime> Type;
};

} 

template<typename VectorsType, typename CoeffsType, int Side> class HouseholderSequence
  : public EigenBase<HouseholderSequence<VectorsType,CoeffsType,Side> >
{
    typedef typename internal::hseq_side_dependent_impl<VectorsType,CoeffsType,Side>::EssentialVectorType EssentialVectorType;

  public:
    enum {
      RowsAtCompileTime = internal::traits<HouseholderSequence>::RowsAtCompileTime,
      ColsAtCompileTime = internal::traits<HouseholderSequence>::ColsAtCompileTime,
      MaxRowsAtCompileTime = internal::traits<HouseholderSequence>::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = internal::traits<HouseholderSequence>::MaxColsAtCompileTime
    };
    typedef typename internal::traits<HouseholderSequence>::Scalar Scalar;

    typedef HouseholderSequence<
      typename internal::conditional<NumTraits<Scalar>::IsComplex,
        typename internal::remove_all<typename VectorsType::ConjugateReturnType>::type,
        VectorsType>::type,
      typename internal::conditional<NumTraits<Scalar>::IsComplex,
        typename internal::remove_all<typename CoeffsType::ConjugateReturnType>::type,
        CoeffsType>::type,
      Side
    > ConjugateReturnType;

    typedef HouseholderSequence<
      VectorsType,
      typename internal::conditional<NumTraits<Scalar>::IsComplex,
        typename internal::remove_all<typename CoeffsType::ConjugateReturnType>::type,
        CoeffsType>::type,
      Side
    > AdjointReturnType;

    typedef HouseholderSequence<
      typename internal::conditional<NumTraits<Scalar>::IsComplex,
        typename internal::remove_all<typename VectorsType::ConjugateReturnType>::type,
        VectorsType>::type,
      CoeffsType,
      Side
    > TransposeReturnType;

    typedef HouseholderSequence<
      typename internal::add_const<VectorsType>::type,
      typename internal::add_const<CoeffsType>::type,
      Side
    > ConstHouseholderSequence;

    EIGEN_DEVICE_FUNC
    HouseholderSequence(const VectorsType& v, const CoeffsType& h)
      : m_vectors(v), m_coeffs(h), m_reverse(false), m_length(v.diagonalSize()),
        m_shift(0)
    {
    }

    EIGEN_DEVICE_FUNC
    HouseholderSequence(const HouseholderSequence& other)
      : m_vectors(other.m_vectors),
        m_coeffs(other.m_coeffs),
        m_reverse(other.m_reverse),
        m_length(other.m_length),
        m_shift(other.m_shift)
    {
    }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    Index rows() const EIGEN_NOEXCEPT { return Side==OnTheLeft ? m_vectors.rows() : m_vectors.cols(); }

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    Index cols() const EIGEN_NOEXCEPT { return rows(); }

    EIGEN_DEVICE_FUNC
    const EssentialVectorType essentialVector(Index k) const
    {
      eigen_assert(k >= 0 && k < m_length);
      return internal::hseq_side_dependent_impl<VectorsType,CoeffsType,Side>::essentialVector(*this, k);
    }

    TransposeReturnType transpose() const
    {
      return TransposeReturnType(m_vectors.conjugate(), m_coeffs)
              .setReverseFlag(!m_reverse)
              .setLength(m_length)
              .setShift(m_shift);
    }

    ConjugateReturnType conjugate() const
    {
      return ConjugateReturnType(m_vectors.conjugate(), m_coeffs.conjugate())
             .setReverseFlag(m_reverse)
             .setLength(m_length)
             .setShift(m_shift);
    }

    template<bool Cond>
    EIGEN_DEVICE_FUNC
    inline typename internal::conditional<Cond,ConjugateReturnType,ConstHouseholderSequence>::type
    conjugateIf() const
    {
      typedef typename internal::conditional<Cond,ConjugateReturnType,ConstHouseholderSequence>::type ReturnType;
      return ReturnType(m_vectors.template conjugateIf<Cond>(), m_coeffs.template conjugateIf<Cond>());
    }

    AdjointReturnType adjoint() const
    {
      return AdjointReturnType(m_vectors, m_coeffs.conjugate())
              .setReverseFlag(!m_reverse)
              .setLength(m_length)
              .setShift(m_shift);
    }

    AdjointReturnType inverse() const { return adjoint(); }

    template<typename DestType>
    inline EIGEN_DEVICE_FUNC
    void evalTo(DestType& dst) const
    {
      Matrix<Scalar, DestType::RowsAtCompileTime, 1,
             AutoAlign|ColMajor, DestType::MaxRowsAtCompileTime, 1> workspace(rows());
      evalTo(dst, workspace);
    }

    template<typename Dest, typename Workspace>
    EIGEN_DEVICE_FUNC
    void evalTo(Dest& dst, Workspace& workspace) const
    {
      workspace.resize(rows());
      Index vecs = m_length;
      if(internal::is_same_dense(dst,m_vectors))
      {
        
        dst.diagonal().setOnes();
        dst.template triangularView<StrictlyUpper>().setZero();
        for(Index k = vecs-1; k >= 0; --k)
        {
          Index cornerSize = rows() - k - m_shift;
          if(m_reverse)
            dst.bottomRightCorner(cornerSize, cornerSize)
               .applyHouseholderOnTheRight(essentialVector(k), m_coeffs.coeff(k), workspace.data());
          else
            dst.bottomRightCorner(cornerSize, cornerSize)
               .applyHouseholderOnTheLeft(essentialVector(k), m_coeffs.coeff(k), workspace.data());

          dst.col(k).tail(rows()-k-1).setZero();
        }
        
        for(Index k = 0; k<cols()-vecs ; ++k)
          dst.col(k).tail(rows()-k-1).setZero();
      }
      else if(m_length>BlockSize)
      {
        dst.setIdentity(rows(), rows());
        if(m_reverse)
          applyThisOnTheLeft(dst,workspace,true);
        else
          applyThisOnTheLeft(dst,workspace,true);
      }
      else
      {
        dst.setIdentity(rows(), rows());
        for(Index k = vecs-1; k >= 0; --k)
        {
          Index cornerSize = rows() - k - m_shift;
          if(m_reverse)
            dst.bottomRightCorner(cornerSize, cornerSize)
               .applyHouseholderOnTheRight(essentialVector(k), m_coeffs.coeff(k), workspace.data());
          else
            dst.bottomRightCorner(cornerSize, cornerSize)
               .applyHouseholderOnTheLeft(essentialVector(k), m_coeffs.coeff(k), workspace.data());
        }
      }
    }

    template<typename Dest> inline void applyThisOnTheRight(Dest& dst) const
    {
      Matrix<Scalar,1,Dest::RowsAtCompileTime,RowMajor,1,Dest::MaxRowsAtCompileTime> workspace(dst.rows());
      applyThisOnTheRight(dst, workspace);
    }

    template<typename Dest, typename Workspace>
    inline void applyThisOnTheRight(Dest& dst, Workspace& workspace) const
    {
      workspace.resize(dst.rows());
      for(Index k = 0; k < m_length; ++k)
      {
        Index actual_k = m_reverse ? m_length-k-1 : k;
        dst.rightCols(rows()-m_shift-actual_k)
           .applyHouseholderOnTheRight(essentialVector(actual_k), m_coeffs.coeff(actual_k), workspace.data());
      }
    }

    template<typename Dest> inline void applyThisOnTheLeft(Dest& dst, bool inputIsIdentity = false) const
    {
      Matrix<Scalar,1,Dest::ColsAtCompileTime,RowMajor,1,Dest::MaxColsAtCompileTime> workspace;
      applyThisOnTheLeft(dst, workspace, inputIsIdentity);
    }

    template<typename Dest, typename Workspace>
    inline void applyThisOnTheLeft(Dest& dst, Workspace& workspace, bool inputIsIdentity = false) const
    {
      if(inputIsIdentity && m_reverse)
        inputIsIdentity = false;
      
      if(m_length>=BlockSize && dst.cols()>1)
      {
        
        Index blockSize = m_length<Index(2*BlockSize) ? (m_length+1)/2 : Index(BlockSize);
        for(Index i = 0; i < m_length; i+=blockSize)
        {
          Index end = m_reverse ? (std::min)(m_length,i+blockSize) : m_length-i;
          Index k = m_reverse ? i : (std::max)(Index(0),end-blockSize);
          Index bs = end-k;
          Index start = k + m_shift;

          typedef Block<typename internal::remove_all<VectorsType>::type,Dynamic,Dynamic> SubVectorsType;
          SubVectorsType sub_vecs1(m_vectors.const_cast_derived(), Side==OnTheRight ? k : start,
                                                                   Side==OnTheRight ? start : k,
                                                                   Side==OnTheRight ? bs : m_vectors.rows()-start,
                                                                   Side==OnTheRight ? m_vectors.cols()-start : bs);
          typename internal::conditional<Side==OnTheRight, Transpose<SubVectorsType>, SubVectorsType&>::type sub_vecs(sub_vecs1);

          Index dstStart = dst.rows()-rows()+m_shift+k;
          Index dstRows  = rows()-m_shift-k;
          Block<Dest,Dynamic,Dynamic> sub_dst(dst,
                                              dstStart,
                                              inputIsIdentity ? dstStart : 0,
                                              dstRows,
                                              inputIsIdentity ? dstRows : dst.cols());
          apply_block_householder_on_the_left(sub_dst, sub_vecs, m_coeffs.segment(k, bs), !m_reverse);
        }
      }
      else
      {
        workspace.resize(dst.cols());
        for(Index k = 0; k < m_length; ++k)
        {
          Index actual_k = m_reverse ? k : m_length-k-1;
          Index dstStart = rows()-m_shift-actual_k;
          dst.bottomRightCorner(dstStart, inputIsIdentity ? dstStart : dst.cols())
            .applyHouseholderOnTheLeft(essentialVector(actual_k), m_coeffs.coeff(actual_k), workspace.data());
        }
      }
    }

    template<typename OtherDerived>
    typename internal::matrix_type_times_scalar_type<Scalar, OtherDerived>::Type operator*(const MatrixBase<OtherDerived>& other) const
    {
      typename internal::matrix_type_times_scalar_type<Scalar, OtherDerived>::Type
        res(other.template cast<typename internal::matrix_type_times_scalar_type<Scalar,OtherDerived>::ResultScalar>());
      applyThisOnTheLeft(res, internal::is_identity<OtherDerived>::value && res.rows()==res.cols());
      return res;
    }

    template<typename VectorsType_, typename CoeffsType_, int Side_> friend struct internal::hseq_side_dependent_impl;

    EIGEN_DEVICE_FUNC
    HouseholderSequence& setLength(Index length)
    {
      m_length = length;
      return *this;
    }

    EIGEN_DEVICE_FUNC
    HouseholderSequence& setShift(Index shift)
    {
      m_shift = shift;
      return *this;
    }

    EIGEN_DEVICE_FUNC
    Index length() const { return m_length; }  

    EIGEN_DEVICE_FUNC
    Index shift() const { return m_shift; }    

    template <typename VectorsType2, typename CoeffsType2, int Side2> friend class HouseholderSequence;

  protected:

    HouseholderSequence& setReverseFlag(bool reverse)
    {
      m_reverse = reverse;
      return *this;
    }

    bool reverseFlag() const { return m_reverse; }     

    typename VectorsType::Nested m_vectors;
    typename CoeffsType::Nested m_coeffs;
    bool m_reverse;
    Index m_length;
    Index m_shift;
    enum { BlockSize = 48 };
};

template<typename OtherDerived, typename VectorsType, typename CoeffsType, int Side>
typename internal::matrix_type_times_scalar_type<typename VectorsType::Scalar,OtherDerived>::Type operator*(const MatrixBase<OtherDerived>& other, const HouseholderSequence<VectorsType,CoeffsType,Side>& h)
{
  typename internal::matrix_type_times_scalar_type<typename VectorsType::Scalar,OtherDerived>::Type
    res(other.template cast<typename internal::matrix_type_times_scalar_type<typename VectorsType::Scalar,OtherDerived>::ResultScalar>());
  h.applyThisOnTheRight(res);
  return res;
}

template<typename VectorsType, typename CoeffsType>
HouseholderSequence<VectorsType,CoeffsType> householderSequence(const VectorsType& v, const CoeffsType& h)
{
  return HouseholderSequence<VectorsType,CoeffsType,OnTheLeft>(v, h);
}

template<typename VectorsType, typename CoeffsType>
HouseholderSequence<VectorsType,CoeffsType,OnTheRight> rightHouseholderSequence(const VectorsType& v, const CoeffsType& h)
{
  return HouseholderSequence<VectorsType,CoeffsType,OnTheRight>(v, h);
}

} 

#endif 
