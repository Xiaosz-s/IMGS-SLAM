#ifndef EIGEN_ROTATIONBASE_H
#define EIGEN_ROTATIONBASE_H

namespace Eigen { 

namespace internal {
template<typename RotationDerived, typename MatrixType, bool IsVector=MatrixType::IsVectorAtCompileTime>
struct rotation_base_generic_product_selector;
}

template<typename Derived, int Dim_>
class RotationBase
{
  public:
    enum { Dim = Dim_ };
    
    typedef typename internal::traits<Derived>::Scalar Scalar;

    typedef Matrix<Scalar,Dim,Dim> RotationMatrixType;
    typedef Matrix<Scalar,Dim,1> VectorType;

  public:
    EIGEN_DEVICE_FUNC inline const Derived& derived() const { return *static_cast<const Derived*>(this); }
    EIGEN_DEVICE_FUNC inline Derived& derived() { return *static_cast<Derived*>(this); }

    EIGEN_DEVICE_FUNC inline RotationMatrixType toRotationMatrix() const { return derived().toRotationMatrix(); }

    EIGEN_DEVICE_FUNC inline RotationMatrixType matrix() const { return derived().toRotationMatrix(); }

    EIGEN_DEVICE_FUNC inline Derived inverse() const { return derived().inverse(); }

    EIGEN_DEVICE_FUNC inline Transform<Scalar,Dim,Isometry> operator*(const Translation<Scalar,Dim>& t) const
    { return Transform<Scalar,Dim,Isometry>(*this) * t; }

    EIGEN_DEVICE_FUNC inline RotationMatrixType operator*(const UniformScaling<Scalar>& s) const
    { return toRotationMatrix() * s.factor(); }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename internal::rotation_base_generic_product_selector<Derived,OtherDerived,OtherDerived::IsVectorAtCompileTime>::ReturnType
    operator*(const EigenBase<OtherDerived>& e) const
    { return internal::rotation_base_generic_product_selector<Derived,OtherDerived>::run(derived(), e.derived()); }

    template<typename OtherDerived> friend
    EIGEN_DEVICE_FUNC inline RotationMatrixType operator*(const EigenBase<OtherDerived>& l, const Derived& r)
    { return l.derived() * r.toRotationMatrix(); }

    EIGEN_DEVICE_FUNC friend inline Transform<Scalar,Dim,Affine> operator*(const DiagonalMatrix<Scalar,Dim>& l, const Derived& r)
    { 
      Transform<Scalar,Dim,Affine> res(r);
      res.linear().applyOnTheLeft(l);
      return res;
    }

    template<int Mode, int Options>
    EIGEN_DEVICE_FUNC inline Transform<Scalar,Dim,Mode> operator*(const Transform<Scalar,Dim,Mode,Options>& t) const
    { return toRotationMatrix() * t; }

    template<typename OtherVectorType>
    EIGEN_DEVICE_FUNC inline VectorType _transformVector(const OtherVectorType& v) const
    { return toRotationMatrix() * v; }
};

namespace internal {

template<typename RotationDerived, typename MatrixType>
struct rotation_base_generic_product_selector<RotationDerived,MatrixType,false>
{
  enum { Dim = RotationDerived::Dim };
  typedef Matrix<typename RotationDerived::Scalar,Dim,Dim> ReturnType;
  EIGEN_DEVICE_FUNC static inline ReturnType run(const RotationDerived& r, const MatrixType& m)
  { return r.toRotationMatrix() * m; }
};

template<typename RotationDerived, typename Scalar, int Dim, int MaxDim>
struct rotation_base_generic_product_selector< RotationDerived, DiagonalMatrix<Scalar,Dim,MaxDim>, false >
{
  typedef Transform<Scalar,Dim,Affine> ReturnType;
  EIGEN_DEVICE_FUNC static inline ReturnType run(const RotationDerived& r, const DiagonalMatrix<Scalar,Dim,MaxDim>& m)
  {
    ReturnType res(r);
    res.linear() *= m;
    return res;
  }
};

template<typename RotationDerived,typename OtherVectorType>
struct rotation_base_generic_product_selector<RotationDerived,OtherVectorType,true>
{
  enum { Dim = RotationDerived::Dim };
  typedef Matrix<typename RotationDerived::Scalar,Dim,1> ReturnType;
  EIGEN_DEVICE_FUNC static EIGEN_STRONG_INLINE ReturnType run(const RotationDerived& r, const OtherVectorType& v)
  {
    return r._transformVector(v);
  }
};

} 

template<typename Scalar_, int Rows_, int Cols_, int _Storage, int MaxRows_, int MaxCols_>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC Matrix<Scalar_, Rows_, Cols_, _Storage, MaxRows_, MaxCols_>
::Matrix(const RotationBase<OtherDerived,ColsAtCompileTime>& r)
{
  EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix,int(OtherDerived::Dim),int(OtherDerived::Dim))
  *this = r.toRotationMatrix();
}

template<typename Scalar_, int Rows_, int Cols_, int _Storage, int MaxRows_, int MaxCols_>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC Matrix<Scalar_, Rows_, Cols_, _Storage, MaxRows_, MaxCols_>&
Matrix<Scalar_, Rows_, Cols_, _Storage, MaxRows_, MaxCols_>
::operator=(const RotationBase<OtherDerived,ColsAtCompileTime>& r)
{
  EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix,int(OtherDerived::Dim),int(OtherDerived::Dim))
  return *this = r.toRotationMatrix();
}

namespace internal {

template<typename Scalar, int Dim>
EIGEN_DEVICE_FUNC static inline Matrix<Scalar,2,2> toRotationMatrix(const Scalar& s)
{
  EIGEN_STATIC_ASSERT(Dim==2,YOU_MADE_A_PROGRAMMING_MISTAKE)
  return Rotation2D<Scalar>(s).toRotationMatrix();
}

template<typename Scalar, int Dim, typename OtherDerived>
EIGEN_DEVICE_FUNC static inline Matrix<Scalar,Dim,Dim> toRotationMatrix(const RotationBase<OtherDerived,Dim>& r)
{
  return r.toRotationMatrix();
}

template<typename Scalar, int Dim, typename OtherDerived>
EIGEN_DEVICE_FUNC static inline const MatrixBase<OtherDerived>& toRotationMatrix(const MatrixBase<OtherDerived>& mat)
{
  EIGEN_STATIC_ASSERT(OtherDerived::RowsAtCompileTime==Dim && OtherDerived::ColsAtCompileTime==Dim,
    YOU_MADE_A_PROGRAMMING_MISTAKE)
  return mat;
}

} 

} 

#endif 
