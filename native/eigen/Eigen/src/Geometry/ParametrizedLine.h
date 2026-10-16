#ifndef EIGEN_PARAMETRIZEDLINE_H
#define EIGEN_PARAMETRIZEDLINE_H

namespace Eigen { 

template <typename Scalar_, int _AmbientDim, int Options_>
class ParametrizedLine
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF_VECTORIZABLE_FIXED_SIZE(Scalar_,_AmbientDim)
  enum {
    AmbientDimAtCompileTime = _AmbientDim,
    Options = Options_
  };
  typedef Scalar_ Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Eigen::Index Index; 
  typedef Matrix<Scalar,AmbientDimAtCompileTime,1,Options> VectorType;

  EIGEN_DEVICE_FUNC inline ParametrizedLine() {}
  
  template<int OtherOptions>
  EIGEN_DEVICE_FUNC ParametrizedLine(const ParametrizedLine<Scalar,AmbientDimAtCompileTime,OtherOptions>& other)
   : m_origin(other.origin()), m_direction(other.direction())
  {}

  EIGEN_DEVICE_FUNC inline explicit ParametrizedLine(Index _dim) : m_origin(_dim), m_direction(_dim) {}

  EIGEN_DEVICE_FUNC ParametrizedLine(const VectorType& origin, const VectorType& direction)
    : m_origin(origin), m_direction(direction) {}

  template <int OtherOptions>
  EIGEN_DEVICE_FUNC explicit ParametrizedLine(const Hyperplane<Scalar_, _AmbientDim, OtherOptions>& hyperplane);

  EIGEN_DEVICE_FUNC static inline ParametrizedLine Through(const VectorType& p0, const VectorType& p1)
  { return ParametrizedLine(p0, (p1-p0).normalized()); }

  EIGEN_DEVICE_FUNC ~ParametrizedLine() {}

  EIGEN_DEVICE_FUNC inline Index dim() const { return m_direction.size(); }

  EIGEN_DEVICE_FUNC const VectorType& origin() const { return m_origin; }
  EIGEN_DEVICE_FUNC VectorType& origin() { return m_origin; }

  EIGEN_DEVICE_FUNC const VectorType& direction() const { return m_direction; }
  EIGEN_DEVICE_FUNC VectorType& direction() { return m_direction; }

  EIGEN_DEVICE_FUNC RealScalar squaredDistance(const VectorType& p) const
  {
    VectorType diff = p - origin();
    return (diff - direction().dot(diff) * direction()).squaredNorm();
  }

  EIGEN_DEVICE_FUNC RealScalar distance(const VectorType& p) const { EIGEN_USING_STD(sqrt) return sqrt(squaredDistance(p)); }

  EIGEN_DEVICE_FUNC VectorType projection(const VectorType& p) const
  { return origin() + direction().dot(p-origin()) * direction(); }

  EIGEN_DEVICE_FUNC VectorType pointAt(const Scalar& t) const;
  
  template <int OtherOptions>
  EIGEN_DEVICE_FUNC Scalar intersectionParameter(const Hyperplane<Scalar_, _AmbientDim, OtherOptions>& hyperplane) const;
 
  template <int OtherOptions>
  EIGEN_DEVICE_FUNC Scalar intersection(const Hyperplane<Scalar_, _AmbientDim, OtherOptions>& hyperplane) const;
  
  template <int OtherOptions>
  EIGEN_DEVICE_FUNC VectorType intersectionPoint(const Hyperplane<Scalar_, _AmbientDim, OtherOptions>& hyperplane) const;

  template<typename XprType>
  EIGEN_DEVICE_FUNC inline ParametrizedLine& transform(const MatrixBase<XprType>& mat, TransformTraits traits = Affine)
  {
    if (traits==Affine)
      direction() = (mat * direction()).normalized();
    else if (traits==Isometry)
      direction() = mat * direction();
    else
    {
      eigen_assert(0 && "invalid traits value in ParametrizedLine::transform()");
    }
    origin() = mat * origin();
    return *this;
  }

  template<int TrOptions>
  EIGEN_DEVICE_FUNC inline ParametrizedLine& transform(const Transform<Scalar,AmbientDimAtCompileTime,Affine,TrOptions>& t,
                                                       TransformTraits traits = Affine)
  {
    transform(t.linear(), traits);
    origin() += t.translation();
    return *this;
  }

  template<typename NewScalarType>
  EIGEN_DEVICE_FUNC inline typename internal::cast_return_type<ParametrizedLine,
           ParametrizedLine<NewScalarType,AmbientDimAtCompileTime,Options> >::type cast() const
  {
    return typename internal::cast_return_type<ParametrizedLine,
                    ParametrizedLine<NewScalarType,AmbientDimAtCompileTime,Options> >::type(*this);
  }

  template<typename OtherScalarType,int OtherOptions>
  EIGEN_DEVICE_FUNC inline explicit ParametrizedLine(const ParametrizedLine<OtherScalarType,AmbientDimAtCompileTime,OtherOptions>& other)
  {
    m_origin = other.origin().template cast<Scalar>();
    m_direction = other.direction().template cast<Scalar>();
  }

  EIGEN_DEVICE_FUNC bool isApprox(const ParametrizedLine& other, const typename NumTraits<Scalar>::Real& prec = NumTraits<Scalar>::dummy_precision()) const
  { return m_origin.isApprox(other.m_origin, prec) && m_direction.isApprox(other.m_direction, prec); }

protected:

  VectorType m_origin, m_direction;
};

template <typename Scalar_, int _AmbientDim, int Options_>
template <int OtherOptions>
EIGEN_DEVICE_FUNC inline ParametrizedLine<Scalar_, _AmbientDim,Options_>::ParametrizedLine(const Hyperplane<Scalar_, _AmbientDim,OtherOptions>& hyperplane)
{
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VectorType, 2)
  direction() = hyperplane.normal().unitOrthogonal();
  origin() = -hyperplane.normal()*hyperplane.offset();
}

template <typename Scalar_, int _AmbientDim, int Options_>
EIGEN_DEVICE_FUNC inline typename ParametrizedLine<Scalar_, _AmbientDim,Options_>::VectorType
ParametrizedLine<Scalar_, _AmbientDim,Options_>::pointAt(const Scalar_& t) const
{
  return origin() + (direction()*t); 
}

template <typename Scalar_, int _AmbientDim, int Options_>
template <int OtherOptions>
EIGEN_DEVICE_FUNC inline Scalar_ ParametrizedLine<Scalar_, _AmbientDim,Options_>::intersectionParameter(const Hyperplane<Scalar_, _AmbientDim, OtherOptions>& hyperplane) const
{
  return -(hyperplane.offset()+hyperplane.normal().dot(origin()))
          / hyperplane.normal().dot(direction());
}

template <typename Scalar_, int _AmbientDim, int Options_>
template <int OtherOptions>
EIGEN_DEVICE_FUNC inline Scalar_ ParametrizedLine<Scalar_, _AmbientDim,Options_>::intersection(const Hyperplane<Scalar_, _AmbientDim, OtherOptions>& hyperplane) const
{
  return intersectionParameter(hyperplane);
}

template <typename Scalar_, int _AmbientDim, int Options_>
template <int OtherOptions>
EIGEN_DEVICE_FUNC inline typename ParametrizedLine<Scalar_, _AmbientDim,Options_>::VectorType
ParametrizedLine<Scalar_, _AmbientDim,Options_>::intersectionPoint(const Hyperplane<Scalar_, _AmbientDim, OtherOptions>& hyperplane) const
{
  return pointAt(intersectionParameter(hyperplane));
}

} 

#endif 
