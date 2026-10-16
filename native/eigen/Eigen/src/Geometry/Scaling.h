#ifndef EIGEN_SCALING_H
#define EIGEN_SCALING_H

namespace Eigen { 

namespace internal
{

  template <typename Scalar, int Dim, int Mode>
  struct uniformscaling_times_affine_returntype
  {
    enum
    {
      NewMode = int(Mode) == int(Isometry) ? Affine : Mode
    };
    typedef Transform <Scalar, Dim, NewMode> type;
  };
}

template<typename Scalar_>
class UniformScaling
{
public:
  
  typedef Scalar_ Scalar;

protected:

  Scalar m_factor;

public:

  UniformScaling() {}
  
  explicit inline UniformScaling(const Scalar& s) : m_factor(s) {}

  inline const Scalar& factor() const { return m_factor; }
  inline Scalar& factor() { return m_factor; }

  inline UniformScaling operator* (const UniformScaling& other) const
  { return UniformScaling(m_factor * other.factor()); }

  template<int Dim>
  inline Transform<Scalar,Dim,Affine> operator* (const Translation<Scalar,Dim>& t) const;

  template<int Dim, int Mode, int Options>
  inline typename
	internal::uniformscaling_times_affine_returntype<Scalar,Dim,Mode>::type
	operator* (const Transform<Scalar, Dim, Mode, Options>& t) const
  {
    typename internal::uniformscaling_times_affine_returntype<Scalar,Dim,Mode>::type res = t;
    res.prescale(factor());
    return res;
  }

  template<typename Derived>
  inline typename Eigen::internal::plain_matrix_type<Derived>::type operator* (const MatrixBase<Derived>& other) const
  { return other * m_factor; }

  template<typename Derived,int Dim>
  inline Matrix<Scalar,Dim,Dim> operator*(const RotationBase<Derived,Dim>& r) const
  { return r.toRotationMatrix() * m_factor; }

  inline UniformScaling inverse() const
  { return UniformScaling(Scalar(1)/m_factor); }

  template<typename NewScalarType>
  inline UniformScaling<NewScalarType> cast() const
  { return UniformScaling<NewScalarType>(NewScalarType(m_factor)); }

  template<typename OtherScalarType>
  inline explicit UniformScaling(const UniformScaling<OtherScalarType>& other)
  { m_factor = Scalar(other.factor()); }

  bool isApprox(const UniformScaling& other, const typename NumTraits<Scalar>::Real& prec = NumTraits<Scalar>::dummy_precision()) const
  { return internal::isApprox(m_factor, other.factor(), prec); }

};

template<typename Derived,typename Scalar>
EIGEN_EXPR_BINARYOP_SCALAR_RETURN_TYPE(Derived,Scalar,product)
operator*(const MatrixBase<Derived>& matrix, const UniformScaling<Scalar>& s)
{ return matrix.derived() * s.factor(); }

inline UniformScaling<float> Scaling(float s) { return UniformScaling<float>(s); }

inline UniformScaling<double> Scaling(double s) { return UniformScaling<double>(s); }

template<typename RealScalar>
inline UniformScaling<std::complex<RealScalar> > Scaling(const std::complex<RealScalar>& s)
{ return UniformScaling<std::complex<RealScalar> >(s); }

template<typename Scalar>
inline DiagonalMatrix<Scalar,2> Scaling(const Scalar& sx, const Scalar& sy)
{ return DiagonalMatrix<Scalar,2>(sx, sy); }

template<typename Scalar>
inline DiagonalMatrix<Scalar,3> Scaling(const Scalar& sx, const Scalar& sy, const Scalar& sz)
{ return DiagonalMatrix<Scalar,3>(sx, sy, sz); }

template<typename Derived>
inline const DiagonalWrapper<const Derived> Scaling(const MatrixBase<Derived>& coeffs)
{ return coeffs.asDiagonal(); }

typedef DiagonalMatrix<float, 2> AlignedScaling2f;

typedef DiagonalMatrix<double,2> AlignedScaling2d;

typedef DiagonalMatrix<float, 3> AlignedScaling3f;

typedef DiagonalMatrix<double,3> AlignedScaling3d;

template<typename Scalar>
template<int Dim>
inline Transform<Scalar,Dim,Affine>
UniformScaling<Scalar>::operator* (const Translation<Scalar,Dim>& t) const
{
  Transform<Scalar,Dim,Affine> res;
  res.matrix().setZero();
  res.linear().diagonal().fill(factor());
  res.translation() = factor() * t.vector();
  res(Dim,Dim) = Scalar(1);
  return res;
}

} 

#endif 
