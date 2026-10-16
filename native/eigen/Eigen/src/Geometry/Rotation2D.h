#ifndef EIGEN_ROTATION2D_H
#define EIGEN_ROTATION2D_H

namespace Eigen { 

namespace internal {

template<typename Scalar_> struct traits<Rotation2D<Scalar_> >
{
  typedef Scalar_ Scalar;
};
} 

template<typename Scalar_>
class Rotation2D : public RotationBase<Rotation2D<Scalar_>,2>
{
  typedef RotationBase<Rotation2D<Scalar_>,2> Base;

public:

  using Base::operator*;

  enum { Dim = 2 };
  
  typedef Scalar_ Scalar;
  typedef Matrix<Scalar,2,1> Vector2;
  typedef Matrix<Scalar,2,2> Matrix2;

protected:

  Scalar m_angle;

public:

  EIGEN_DEVICE_FUNC explicit inline Rotation2D(const Scalar& a) : m_angle(a) {}

  EIGEN_DEVICE_FUNC Rotation2D() {}

  template<typename Derived>
  EIGEN_DEVICE_FUNC explicit Rotation2D(const MatrixBase<Derived>& m)
  {
    fromRotationMatrix(m.derived());
  }

  EIGEN_DEVICE_FUNC inline Scalar angle() const { return m_angle; }

  EIGEN_DEVICE_FUNC inline Scalar& angle() { return m_angle; }

  EIGEN_DEVICE_FUNC inline Scalar smallestPositiveAngle() const {
    Scalar tmp = numext::fmod(m_angle,Scalar(2*EIGEN_PI));
    return tmp<Scalar(0) ? tmp + Scalar(2*EIGEN_PI) : tmp;
  }

  EIGEN_DEVICE_FUNC inline Scalar smallestAngle() const {
    Scalar tmp = numext::fmod(m_angle,Scalar(2*EIGEN_PI));
    if(tmp>Scalar(EIGEN_PI))       tmp -= Scalar(2*EIGEN_PI);
    else if(tmp<-Scalar(EIGEN_PI)) tmp += Scalar(2*EIGEN_PI);
    return tmp;
  }

  EIGEN_DEVICE_FUNC inline Rotation2D inverse() const { return Rotation2D(-m_angle); }

  EIGEN_DEVICE_FUNC inline Rotation2D operator*(const Rotation2D& other) const
  { return Rotation2D(m_angle + other.m_angle); }

  EIGEN_DEVICE_FUNC inline Rotation2D& operator*=(const Rotation2D& other)
  { m_angle += other.m_angle; return *this; }

  EIGEN_DEVICE_FUNC Vector2 operator* (const Vector2& vec) const
  { return toRotationMatrix() * vec; }
  
  template<typename Derived>
  EIGEN_DEVICE_FUNC Rotation2D& fromRotationMatrix(const MatrixBase<Derived>& m);
  EIGEN_DEVICE_FUNC Matrix2 toRotationMatrix() const;

  template<typename Derived>
  EIGEN_DEVICE_FUNC Rotation2D& operator=(const  MatrixBase<Derived>& m)
  { return fromRotationMatrix(m.derived()); }

  EIGEN_DEVICE_FUNC inline Rotation2D slerp(const Scalar& t, const Rotation2D& other) const
  {
    Scalar dist = Rotation2D(other.m_angle-m_angle).smallestAngle();
    return Rotation2D(m_angle + dist*t);
  }

  template<typename NewScalarType>
  EIGEN_DEVICE_FUNC inline typename internal::cast_return_type<Rotation2D,Rotation2D<NewScalarType> >::type cast() const
  { return typename internal::cast_return_type<Rotation2D,Rotation2D<NewScalarType> >::type(*this); }

  template<typename OtherScalarType>
  EIGEN_DEVICE_FUNC inline explicit Rotation2D(const Rotation2D<OtherScalarType>& other)
  {
    m_angle = Scalar(other.angle());
  }

  EIGEN_DEVICE_FUNC static inline Rotation2D Identity() { return Rotation2D(0); }

  EIGEN_DEVICE_FUNC bool isApprox(const Rotation2D& other, const typename NumTraits<Scalar>::Real& prec = NumTraits<Scalar>::dummy_precision()) const
  { return internal::isApprox(m_angle,other.m_angle, prec); }
  
};

typedef Rotation2D<float> Rotation2Df;

typedef Rotation2D<double> Rotation2Dd;

template<typename Scalar>
template<typename Derived>
EIGEN_DEVICE_FUNC Rotation2D<Scalar>& Rotation2D<Scalar>::fromRotationMatrix(const MatrixBase<Derived>& mat)
{
  EIGEN_USING_STD(atan2)
  EIGEN_STATIC_ASSERT(Derived::RowsAtCompileTime==2 && Derived::ColsAtCompileTime==2,YOU_MADE_A_PROGRAMMING_MISTAKE)
  m_angle = atan2(mat.coeff(1,0), mat.coeff(0,0));
  return *this;
}

template<typename Scalar>
typename Rotation2D<Scalar>::Matrix2
EIGEN_DEVICE_FUNC Rotation2D<Scalar>::toRotationMatrix(void) const
{
  EIGEN_USING_STD(sin)
  EIGEN_USING_STD(cos)
  Scalar sinA = sin(m_angle);
  Scalar cosA = cos(m_angle);
  return (Matrix2() << cosA, -sinA, sinA, cosA).finished();
}

} 

#endif 
