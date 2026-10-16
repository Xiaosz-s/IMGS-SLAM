#ifndef EIGEN_ALIGNEDBOX_H
#define EIGEN_ALIGNEDBOX_H

namespace Eigen {

template <typename Scalar_, int _AmbientDim>
class AlignedBox
{
public:
EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF_VECTORIZABLE_FIXED_SIZE(Scalar_,_AmbientDim)
  enum { AmbientDimAtCompileTime = _AmbientDim };
  typedef Scalar_                                   Scalar;
  typedef NumTraits<Scalar>                         ScalarTraits;
  typedef Eigen::Index                              Index; 
  typedef typename ScalarTraits::Real               RealScalar;
  typedef typename ScalarTraits::NonInteger         NonInteger;
  typedef Matrix<Scalar,AmbientDimAtCompileTime,1>  VectorType;
  typedef CwiseBinaryOp<internal::scalar_sum_op<Scalar>, const VectorType, const VectorType> VectorTypeSum;

  enum CornerType
  {
    
    Min=0, Max=1,

    BottomLeft=0, BottomRight=1,
    TopLeft=2, TopRight=3,

    BottomLeftFloor=0, BottomRightFloor=1,
    TopLeftFloor=2, TopRightFloor=3,
    BottomLeftCeil=4, BottomRightCeil=5,
    TopLeftCeil=6, TopRightCeil=7
    
  };

  EIGEN_DEVICE_FUNC inline AlignedBox()
  { if (EIGEN_CONST_CONDITIONAL(AmbientDimAtCompileTime!=Dynamic)) setEmpty(); }

  EIGEN_DEVICE_FUNC inline explicit AlignedBox(Index _dim) : m_min(_dim), m_max(_dim)
  { setEmpty(); }

  template<typename OtherVectorType1, typename OtherVectorType2>
  EIGEN_DEVICE_FUNC inline AlignedBox(const OtherVectorType1& _min, const OtherVectorType2& _max) : m_min(_min), m_max(_max) {}

  template<typename Derived>
  EIGEN_DEVICE_FUNC inline explicit AlignedBox(const MatrixBase<Derived>& p) : m_min(p), m_max(m_min)
  { }

  EIGEN_DEVICE_FUNC ~AlignedBox() {}

  EIGEN_DEVICE_FUNC inline Index dim() const { return AmbientDimAtCompileTime==Dynamic ? m_min.size() : Index(AmbientDimAtCompileTime); }

  EIGEN_DEVICE_FUNC inline bool isNull() const { return isEmpty(); }

  EIGEN_DEVICE_FUNC inline void setNull() { setEmpty(); }

  EIGEN_DEVICE_FUNC inline bool isEmpty() const { return (m_min.array() > m_max.array()).any(); }

  EIGEN_DEVICE_FUNC inline void setEmpty()
  {
    m_min.setConstant( ScalarTraits::highest() );
    m_max.setConstant( ScalarTraits::lowest() );
  }

  EIGEN_DEVICE_FUNC inline const VectorType& (min)() const { return m_min; }
  
  EIGEN_DEVICE_FUNC inline VectorType& (min)() { return m_min; }
  
  EIGEN_DEVICE_FUNC inline const VectorType& (max)() const { return m_max; }
  
  EIGEN_DEVICE_FUNC inline VectorType& (max)() { return m_max; }

  EIGEN_DEVICE_FUNC inline const EIGEN_EXPR_BINARYOP_SCALAR_RETURN_TYPE(VectorTypeSum, RealScalar, quotient)
  center() const
  { return (m_min+m_max)/RealScalar(2); }

  EIGEN_DEVICE_FUNC inline const CwiseBinaryOp< internal::scalar_difference_op<Scalar,Scalar>, const VectorType, const VectorType> sizes() const
  { return m_max - m_min; }

  EIGEN_DEVICE_FUNC inline Scalar volume() const
  { return sizes().prod(); }

  EIGEN_DEVICE_FUNC inline CwiseBinaryOp< internal::scalar_difference_op<Scalar,Scalar>, const VectorType, const VectorType> diagonal() const
  { return sizes(); }

  EIGEN_DEVICE_FUNC inline VectorType corner(CornerType corner) const
  {
    EIGEN_STATIC_ASSERT(_AmbientDim <= 3, THIS_METHOD_IS_ONLY_FOR_VECTORS_OF_A_SPECIFIC_SIZE);

    VectorType res;

    Index mult = 1;
    for(Index d=0; d<dim(); ++d)
    {
      if( mult & corner ) res[d] = m_max[d];
      else                res[d] = m_min[d];
      mult *= 2;
    }
    return res;
  }

  EIGEN_DEVICE_FUNC inline VectorType sample() const
  {
    VectorType r(dim());
    for(Index d=0; d<dim(); ++d)
    {
      if(!ScalarTraits::IsInteger)
      {
        r[d] = m_min[d] + (m_max[d]-m_min[d])
             * internal::random<Scalar>(Scalar(0), Scalar(1));
      }
      else
        r[d] = internal::random(m_min[d], m_max[d]);
    }
    return r;
  }

  template<typename Derived>
  EIGEN_DEVICE_FUNC inline bool contains(const MatrixBase<Derived>& p) const
  {
    typename internal::nested_eval<Derived,2>::type p_n(p.derived());
    return (m_min.array()<=p_n.array()).all() && (p_n.array()<=m_max.array()).all();
  }

  EIGEN_DEVICE_FUNC inline bool contains(const AlignedBox& b) const
  { return (m_min.array()<=(b.min)().array()).all() && ((b.max)().array()<=m_max.array()).all(); }

  EIGEN_DEVICE_FUNC inline bool intersects(const AlignedBox& b) const
  { return (m_min.array()<=(b.max)().array()).all() && ((b.min)().array()<=m_max.array()).all(); }

  template<typename Derived>
  EIGEN_DEVICE_FUNC inline AlignedBox& extend(const MatrixBase<Derived>& p)
  {
    typename internal::nested_eval<Derived,2>::type p_n(p.derived());
    m_min = m_min.cwiseMin(p_n);
    m_max = m_max.cwiseMax(p_n);
    return *this;
  }

  EIGEN_DEVICE_FUNC inline AlignedBox& extend(const AlignedBox& b)
  {
    m_min = m_min.cwiseMin(b.m_min);
    m_max = m_max.cwiseMax(b.m_max);
    return *this;
  }

  EIGEN_DEVICE_FUNC inline AlignedBox& clamp(const AlignedBox& b)
  {
    m_min = m_min.cwiseMax(b.m_min);
    m_max = m_max.cwiseMin(b.m_max);
    return *this;
  }

  EIGEN_DEVICE_FUNC inline AlignedBox intersection(const AlignedBox& b) const
  {return AlignedBox(m_min.cwiseMax(b.m_min), m_max.cwiseMin(b.m_max)); }

  EIGEN_DEVICE_FUNC inline AlignedBox merged(const AlignedBox& b) const
  { return AlignedBox(m_min.cwiseMin(b.m_min), m_max.cwiseMax(b.m_max)); }

  template<typename Derived>
  EIGEN_DEVICE_FUNC inline AlignedBox& translate(const MatrixBase<Derived>& a_t)
  {
    const typename internal::nested_eval<Derived,2>::type t(a_t.derived());
    m_min += t;
    m_max += t;
    return *this;
  }

  template<typename Derived>
  EIGEN_DEVICE_FUNC inline AlignedBox translated(const MatrixBase<Derived>& a_t) const
  {
    AlignedBox result(m_min, m_max);
    result.translate(a_t);
    return result;
  }

  template<typename Derived>
  EIGEN_DEVICE_FUNC inline Scalar squaredExteriorDistance(const MatrixBase<Derived>& p) const;

  EIGEN_DEVICE_FUNC inline Scalar squaredExteriorDistance(const AlignedBox& b) const;

  template<typename Derived>
  EIGEN_DEVICE_FUNC inline NonInteger exteriorDistance(const MatrixBase<Derived>& p) const
  { EIGEN_USING_STD(sqrt) return sqrt(NonInteger(squaredExteriorDistance(p))); }

  EIGEN_DEVICE_FUNC inline NonInteger exteriorDistance(const AlignedBox& b) const
  { EIGEN_USING_STD(sqrt) return sqrt(NonInteger(squaredExteriorDistance(b))); }

  template<int Mode, int Options>
  EIGEN_DEVICE_FUNC inline void transform(
      const typename Transform<Scalar, AmbientDimAtCompileTime, Mode, Options>::TranslationType& translation)
  {
    this->translate(translation);
  }

  template<int Mode, int Options>
  EIGEN_DEVICE_FUNC inline void transform(const Transform<Scalar, AmbientDimAtCompileTime, Mode, Options>& transform)
  {
    
    EIGEN_STATIC_ASSERT(Mode == Affine || Mode == AffineCompact || Mode == Isometry, THIS_METHOD_IS_ONLY_FOR_SPECIFIC_TRANSFORMATIONS);

    const VectorType rotated_extent_2 = transform.linear().cwiseAbs() * sizes();
    
    const VectorType rotated_center_2 = transform.linear() * (this->m_max + this->m_min) +
        Scalar(2) * transform.translation();

    this->m_max = (rotated_center_2 + rotated_extent_2) / Scalar(2);
    this->m_min = (rotated_center_2 - rotated_extent_2) / Scalar(2);
  }

  template<int Mode, int Options>
  EIGEN_DEVICE_FUNC AlignedBox transformed(const Transform<Scalar, AmbientDimAtCompileTime, Mode, Options>& transform) const
  {
    AlignedBox result(m_min, m_max);
    result.transform(transform);
    return result;
  }

  template<typename NewScalarType>
  EIGEN_DEVICE_FUNC inline typename internal::cast_return_type<AlignedBox,
           AlignedBox<NewScalarType,AmbientDimAtCompileTime> >::type cast() const
  {
    return typename internal::cast_return_type<AlignedBox,
                    AlignedBox<NewScalarType,AmbientDimAtCompileTime> >::type(*this);
  }

  template<typename OtherScalarType>
  EIGEN_DEVICE_FUNC inline explicit AlignedBox(const AlignedBox<OtherScalarType,AmbientDimAtCompileTime>& other)
  {
    m_min = (other.min)().template cast<Scalar>();
    m_max = (other.max)().template cast<Scalar>();
  }

  EIGEN_DEVICE_FUNC bool isApprox(const AlignedBox& other, const RealScalar& prec = ScalarTraits::dummy_precision()) const
  { return m_min.isApprox(other.m_min, prec) && m_max.isApprox(other.m_max, prec); }

protected:

  VectorType m_min, m_max;
};

template<typename Scalar,int AmbientDim>
template<typename Derived>
EIGEN_DEVICE_FUNC inline Scalar AlignedBox<Scalar,AmbientDim>::squaredExteriorDistance(const MatrixBase<Derived>& a_p) const
{
  typename internal::nested_eval<Derived,2*AmbientDim>::type p(a_p.derived());
  Scalar dist2(0);
  Scalar aux;
  for (Index k=0; k<dim(); ++k)
  {
    if( m_min[k] > p[k] )
    {
      aux = m_min[k] - p[k];
      dist2 += aux*aux;
    }
    else if( p[k] > m_max[k] )
    {
      aux = p[k] - m_max[k];
      dist2 += aux*aux;
    }
  }
  return dist2;
}

template<typename Scalar,int AmbientDim>
EIGEN_DEVICE_FUNC inline Scalar AlignedBox<Scalar,AmbientDim>::squaredExteriorDistance(const AlignedBox& b) const
{
  Scalar dist2(0);
  Scalar aux;
  for (Index k=0; k<dim(); ++k)
  {
    if( m_min[k] > b.m_max[k] )
    {
      aux = m_min[k] - b.m_max[k];
      dist2 += aux*aux;
    }
    else if( b.m_min[k] > m_max[k] )
    {
      aux = b.m_min[k] - m_max[k];
      dist2 += aux*aux;
    }
  }
  return dist2;
}

#define EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, Size, SizeSuffix)    \
                                 \
typedef AlignedBox<Type, Size>   AlignedBox##SizeSuffix##TypeSuffix;

#define EIGEN_MAKE_TYPEDEFS_ALL_SIZES(Type, TypeSuffix) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 1, 1) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 2, 2) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 3, 3) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 4, 4) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, Dynamic, X)

EIGEN_MAKE_TYPEDEFS_ALL_SIZES(int,                  i)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(float,                f)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(double,               d)

#undef EIGEN_MAKE_TYPEDEFS_ALL_SIZES
#undef EIGEN_MAKE_TYPEDEFS

} 

#endif 
