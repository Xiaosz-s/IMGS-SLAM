#ifndef EIGEN_TRANSFORM_H
#define EIGEN_TRANSFORM_H

namespace Eigen {

namespace internal {

template<typename Transform>
struct transform_traits
{
  enum
  {
    Dim = Transform::Dim,
    HDim = Transform::HDim,
    Mode = Transform::Mode,
    IsProjective = (int(Mode)==int(Projective))
  };
};

template< typename TransformType,
          typename MatrixType,
          int Case = transform_traits<TransformType>::IsProjective ? 0
                   : int(MatrixType::RowsAtCompileTime) == int(transform_traits<TransformType>::HDim) ? 1
                   : 2,
          int RhsCols = MatrixType::ColsAtCompileTime>
struct transform_right_product_impl;

template< typename Other,
          int Mode,
          int Options,
          int Dim,
          int HDim,
          int OtherRows=Other::RowsAtCompileTime,
          int OtherCols=Other::ColsAtCompileTime>
struct transform_left_product_impl;

template< typename Lhs,
          typename Rhs,
          bool AnyProjective =
            transform_traits<Lhs>::IsProjective ||
            transform_traits<Rhs>::IsProjective>
struct transform_transform_product_impl;

template< typename Other,
          int Mode,
          int Options,
          int Dim,
          int HDim,
          int OtherRows=Other::RowsAtCompileTime,
          int OtherCols=Other::ColsAtCompileTime>
struct transform_construct_from_matrix;

template<typename TransformType> struct transform_take_affine_part;

template<typename Scalar_, int Dim_, int _Mode, int Options_>
struct traits<Transform<Scalar_,Dim_,_Mode,Options_> >
{
  typedef Scalar_ Scalar;
  typedef Eigen::Index StorageIndex;
  typedef Dense StorageKind;
  enum {
    Dim1 = Dim_==Dynamic ? Dim_ : Dim_ + 1,
    RowsAtCompileTime = _Mode==Projective ? Dim1 : Dim_,
    ColsAtCompileTime = Dim1,
    MaxRowsAtCompileTime = RowsAtCompileTime,
    MaxColsAtCompileTime = ColsAtCompileTime,
    Flags = 0
  };
};

template<int Mode> struct transform_make_affine;

} 

template<typename Scalar_, int Dim_, int _Mode, int Options_>
class Transform
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF_VECTORIZABLE_FIXED_SIZE(Scalar_,Dim_==Dynamic ? Dynamic : (Dim_+1)*(Dim_+1))
  enum {
    Mode = _Mode,
    Options = Options_,
    Dim = Dim_,     
    HDim = Dim_+1,  
    Rows = int(Mode)==(AffineCompact) ? Dim : HDim
  };
  
  typedef Scalar_ Scalar;
  typedef Eigen::Index StorageIndex;
  typedef Eigen::Index Index; 
  
  typedef typename internal::make_proper_matrix_type<Scalar,Rows,HDim,Options>::type MatrixType;
  
  typedef const MatrixType ConstMatrixType;
  
  typedef Matrix<Scalar,Dim,Dim,Options> LinearMatrixType;
  
  typedef Block<MatrixType,Dim,Dim,int(Mode)==(AffineCompact) && (int(Options)&RowMajor)==0> LinearPart;
  
  typedef const Block<ConstMatrixType,Dim,Dim,int(Mode)==(AffineCompact) && (int(Options)&RowMajor)==0> ConstLinearPart;
  
  typedef typename internal::conditional<int(Mode)==int(AffineCompact),
                              MatrixType&,
                              Block<MatrixType,Dim,HDim> >::type AffinePart;
  
  typedef typename internal::conditional<int(Mode)==int(AffineCompact),
                              const MatrixType&,
                              const Block<const MatrixType,Dim,HDim> >::type ConstAffinePart;
  
  typedef Matrix<Scalar,Dim,1> VectorType;
  
  typedef Block<MatrixType,Dim,1,!(internal::traits<MatrixType>::Flags & RowMajorBit)> TranslationPart;
  
  typedef const Block<ConstMatrixType,Dim,1,!(internal::traits<MatrixType>::Flags & RowMajorBit)> ConstTranslationPart;
  
  typedef Translation<Scalar,Dim> TranslationType;

  enum { TransformTimeDiagonalMode = ((Mode==int(Isometry))?Affine:int(Mode)) };
  
  typedef Transform<Scalar,Dim,TransformTimeDiagonalMode> TransformTimeDiagonalReturnType;

protected:

  MatrixType m_matrix;

public:

  EIGEN_DEVICE_FUNC inline Transform()
  {
    check_template_params();
    internal::transform_make_affine<(int(Mode)==Affine || int(Mode)==Isometry) ? Affine : AffineCompact>::run(m_matrix);
  }

  EIGEN_DEVICE_FUNC inline explicit Transform(const TranslationType& t)
  {
    check_template_params();
    *this = t;
  }
  EIGEN_DEVICE_FUNC inline explicit Transform(const UniformScaling<Scalar>& s)
  {
    check_template_params();
    *this = s;
  }
  template<typename Derived>
  EIGEN_DEVICE_FUNC inline explicit Transform(const RotationBase<Derived, Dim>& r)
  {
    check_template_params();
    *this = r;
  }

  typedef internal::transform_take_affine_part<Transform> take_affine_part;

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC inline explicit Transform(const EigenBase<OtherDerived>& other)
  {
    EIGEN_STATIC_ASSERT((internal::is_same<Scalar,typename OtherDerived::Scalar>::value),
      YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY);

    check_template_params();
    internal::transform_construct_from_matrix<OtherDerived,Mode,Options,Dim,HDim>::run(this, other.derived());
  }

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC inline Transform& operator=(const EigenBase<OtherDerived>& other)
  {
    EIGEN_STATIC_ASSERT((internal::is_same<Scalar,typename OtherDerived::Scalar>::value),
      YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY);

    internal::transform_construct_from_matrix<OtherDerived,Mode,Options,Dim,HDim>::run(this, other.derived());
    return *this;
  }

  template<int OtherOptions>
  EIGEN_DEVICE_FUNC inline Transform(const Transform<Scalar,Dim,Mode,OtherOptions>& other)
  {
    check_template_params();
    
    m_matrix = other.matrix();
  }

  template<int OtherMode,int OtherOptions>
  EIGEN_DEVICE_FUNC inline Transform(const Transform<Scalar,Dim,OtherMode,OtherOptions>& other)
  {
    check_template_params();

    EIGEN_STATIC_ASSERT(EIGEN_IMPLIES(OtherMode==int(Projective), Mode==int(Projective)),
                        YOU_PERFORMED_AN_INVALID_TRANSFORMATION_CONVERSION)

    EIGEN_STATIC_ASSERT(EIGEN_IMPLIES(OtherMode==int(Affine)||OtherMode==int(AffineCompact), Mode!=int(Isometry)),
                        YOU_PERFORMED_AN_INVALID_TRANSFORMATION_CONVERSION)

    enum { ModeIsAffineCompact = Mode == int(AffineCompact),
           OtherModeIsAffineCompact = OtherMode == int(AffineCompact)
    };

    if(EIGEN_CONST_CONDITIONAL(ModeIsAffineCompact == OtherModeIsAffineCompact))
    {

      m_matrix.template block<Dim,Dim+1>(0,0) = other.matrix().template block<Dim,Dim+1>(0,0);
      makeAffine();
    }
    else if(EIGEN_CONST_CONDITIONAL(OtherModeIsAffineCompact))
    {
      typedef typename Transform<Scalar,Dim,OtherMode,OtherOptions>::MatrixType OtherMatrixType;
      internal::transform_construct_from_matrix<OtherMatrixType,Mode,Options,Dim,HDim>::run(this, other.matrix());
    }
    else
    {

      linear() = other.linear();
      translation() = other.translation();
    }
  }

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC Transform(const ReturnByValue<OtherDerived>& other)
  {
    check_template_params();
    other.evalTo(*this);
  }

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC Transform& operator=(const ReturnByValue<OtherDerived>& other)
  {
    other.evalTo(*this);
    return *this;
  }

  #ifdef EIGEN_QT_SUPPORT
  inline Transform(const QMatrix& other);
  inline Transform& operator=(const QMatrix& other);
  inline QMatrix toQMatrix(void) const;
  inline Transform(const QTransform& other);
  inline Transform& operator=(const QTransform& other);
  inline QTransform toQTransform(void) const;
  #endif

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return int(Mode)==int(Projective) ? m_matrix.cols() : (m_matrix.cols()-1); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return m_matrix.cols(); }

  EIGEN_DEVICE_FUNC inline Scalar operator() (Index row, Index col) const { return m_matrix(row,col); }

  EIGEN_DEVICE_FUNC inline Scalar& operator() (Index row, Index col) { return m_matrix(row,col); }

  EIGEN_DEVICE_FUNC inline const MatrixType& matrix() const { return m_matrix; }
  
  EIGEN_DEVICE_FUNC inline MatrixType& matrix() { return m_matrix; }

  EIGEN_DEVICE_FUNC inline ConstLinearPart linear() const { return ConstLinearPart(m_matrix,0,0); }
  
  EIGEN_DEVICE_FUNC inline LinearPart linear() { return LinearPart(m_matrix,0,0); }

  EIGEN_DEVICE_FUNC inline ConstAffinePart affine() const { return take_affine_part::run(m_matrix); }
  
  EIGEN_DEVICE_FUNC inline AffinePart affine() { return take_affine_part::run(m_matrix); }

  EIGEN_DEVICE_FUNC inline ConstTranslationPart translation() const { return ConstTranslationPart(m_matrix,0,Dim); }
  
  EIGEN_DEVICE_FUNC inline TranslationPart translation() { return TranslationPart(m_matrix,0,Dim); }

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const typename internal::transform_right_product_impl<Transform, OtherDerived>::ResultType
  operator * (const EigenBase<OtherDerived> &other) const
  { return internal::transform_right_product_impl<Transform, OtherDerived>::run(*this,other.derived()); }

  template<typename OtherDerived> friend
  EIGEN_DEVICE_FUNC inline const typename internal::transform_left_product_impl<OtherDerived,Mode,Options,Dim_,Dim_+1>::ResultType
    operator * (const EigenBase<OtherDerived> &a, const Transform &b)
  { return internal::transform_left_product_impl<OtherDerived,Mode,Options,Dim,HDim>::run(a.derived(),b); }

  template<typename DiagonalDerived>
  EIGEN_DEVICE_FUNC inline const TransformTimeDiagonalReturnType
    operator * (const DiagonalBase<DiagonalDerived> &b) const
  {
    TransformTimeDiagonalReturnType res(*this);
    res.linearExt() *= b;
    return res;
  }

  template<typename DiagonalDerived>
  EIGEN_DEVICE_FUNC friend inline TransformTimeDiagonalReturnType
    operator * (const DiagonalBase<DiagonalDerived> &a, const Transform &b)
  {
    TransformTimeDiagonalReturnType res;
    res.linear().noalias() = a*b.linear();
    res.translation().noalias() = a*b.translation();
    if (EIGEN_CONST_CONDITIONAL(Mode!=int(AffineCompact)))
      res.matrix().row(Dim) = b.matrix().row(Dim);
    return res;
  }

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC inline Transform& operator*=(const EigenBase<OtherDerived>& other) { return *this = *this * other; }

  EIGEN_DEVICE_FUNC inline const Transform operator * (const Transform& other) const
  {
    return internal::transform_transform_product_impl<Transform,Transform>::run(*this,other);
  }

  #if EIGEN_COMP_ICC
private:

  template<int OtherMode,int OtherOptions> struct icc_11_workaround
  {
    typedef internal::transform_transform_product_impl<Transform,Transform<Scalar,Dim,OtherMode,OtherOptions> > ProductType;
    typedef typename ProductType::ResultType ResultType;
  };

public:
  
  template<int OtherMode,int OtherOptions>
  inline typename icc_11_workaround<OtherMode,OtherOptions>::ResultType
    operator * (const Transform<Scalar,Dim,OtherMode,OtherOptions>& other) const
  {
    typedef typename icc_11_workaround<OtherMode,OtherOptions>::ProductType ProductType;
    return ProductType::run(*this,other);
  }
  #else
  
  template<int OtherMode,int OtherOptions>
  EIGEN_DEVICE_FUNC inline typename internal::transform_transform_product_impl<Transform,Transform<Scalar,Dim,OtherMode,OtherOptions> >::ResultType
    operator * (const Transform<Scalar,Dim,OtherMode,OtherOptions>& other) const
  {
    return internal::transform_transform_product_impl<Transform,Transform<Scalar,Dim,OtherMode,OtherOptions> >::run(*this,other);
  }
  #endif

  EIGEN_DEVICE_FUNC void setIdentity() { m_matrix.setIdentity(); }

  EIGEN_DEVICE_FUNC static const Transform Identity()
  {
    return Transform(MatrixType::Identity());
  }

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC
  inline Transform& scale(const MatrixBase<OtherDerived> &other);

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC
  inline Transform& prescale(const MatrixBase<OtherDerived> &other);

  EIGEN_DEVICE_FUNC inline Transform& scale(const Scalar& s);
  EIGEN_DEVICE_FUNC inline Transform& prescale(const Scalar& s);

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC
  inline Transform& translate(const MatrixBase<OtherDerived> &other);

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC
  inline Transform& pretranslate(const MatrixBase<OtherDerived> &other);

  template<typename RotationType>
  EIGEN_DEVICE_FUNC
  inline Transform& rotate(const RotationType& rotation);

  template<typename RotationType>
  EIGEN_DEVICE_FUNC
  inline Transform& prerotate(const RotationType& rotation);

  EIGEN_DEVICE_FUNC Transform& shear(const Scalar& sx, const Scalar& sy);
  EIGEN_DEVICE_FUNC Transform& preshear(const Scalar& sx, const Scalar& sy);

  EIGEN_DEVICE_FUNC inline Transform& operator=(const TranslationType& t);

  EIGEN_DEVICE_FUNC
  inline Transform& operator*=(const TranslationType& t) { return translate(t.vector()); }

  EIGEN_DEVICE_FUNC inline Transform operator*(const TranslationType& t) const;

  EIGEN_DEVICE_FUNC
  inline Transform& operator=(const UniformScaling<Scalar>& t);

  EIGEN_DEVICE_FUNC
  inline Transform& operator*=(const UniformScaling<Scalar>& s) { return scale(s.factor()); }

  EIGEN_DEVICE_FUNC
  inline TransformTimeDiagonalReturnType operator*(const UniformScaling<Scalar>& s) const
  {
    TransformTimeDiagonalReturnType res = *this;
    res.scale(s.factor());
    return res;
  }

  EIGEN_DEVICE_FUNC
  inline Transform& operator*=(const DiagonalMatrix<Scalar,Dim>& s) { linearExt() *= s; return *this; }

  template<typename Derived>
  EIGEN_DEVICE_FUNC inline Transform& operator=(const RotationBase<Derived,Dim>& r);
  template<typename Derived>
  EIGEN_DEVICE_FUNC inline Transform& operator*=(const RotationBase<Derived,Dim>& r) { return rotate(r.toRotationMatrix()); }
  template<typename Derived>
  EIGEN_DEVICE_FUNC inline Transform operator*(const RotationBase<Derived,Dim>& r) const;

  typedef typename internal::conditional<int(Mode)==Isometry,ConstLinearPart,const LinearMatrixType>::type RotationReturnType;
  EIGEN_DEVICE_FUNC RotationReturnType rotation() const;

  template<typename RotationMatrixType, typename ScalingMatrixType>
  EIGEN_DEVICE_FUNC
  void computeRotationScaling(RotationMatrixType *rotation, ScalingMatrixType *scaling) const;
  template<typename ScalingMatrixType, typename RotationMatrixType>
  EIGEN_DEVICE_FUNC
  void computeScalingRotation(ScalingMatrixType *scaling, RotationMatrixType *rotation) const;

  template<typename PositionDerived, typename OrientationType, typename ScaleDerived>
  EIGEN_DEVICE_FUNC
  Transform& fromPositionOrientationScale(const MatrixBase<PositionDerived> &position,
    const OrientationType& orientation, const MatrixBase<ScaleDerived> &scale);

  EIGEN_DEVICE_FUNC
  inline Transform inverse(TransformTraits traits = (TransformTraits)Mode) const;

  EIGEN_DEVICE_FUNC const Scalar* data() const { return m_matrix.data(); }
  
  EIGEN_DEVICE_FUNC Scalar* data() { return m_matrix.data(); }

  template<typename NewScalarType>
  EIGEN_DEVICE_FUNC inline typename internal::cast_return_type<Transform,Transform<NewScalarType,Dim,Mode,Options> >::type cast() const
  { return typename internal::cast_return_type<Transform,Transform<NewScalarType,Dim,Mode,Options> >::type(*this); }

  template<typename OtherScalarType>
  EIGEN_DEVICE_FUNC inline explicit Transform(const Transform<OtherScalarType,Dim,Mode,Options>& other)
  {
    check_template_params();
    m_matrix = other.matrix().template cast<Scalar>();
  }

  EIGEN_DEVICE_FUNC bool isApprox(const Transform& other, const typename NumTraits<Scalar>::Real& prec = NumTraits<Scalar>::dummy_precision()) const
  { return m_matrix.isApprox(other.m_matrix, prec); }

  EIGEN_DEVICE_FUNC void makeAffine()
  {
    internal::transform_make_affine<int(Mode)>::run(m_matrix);
  }

  EIGEN_DEVICE_FUNC inline Block<MatrixType,int(Mode)==int(Projective)?HDim:Dim,Dim> linearExt()
  { return m_matrix.template block<int(Mode)==int(Projective)?HDim:Dim,Dim>(0,0); }

  EIGEN_DEVICE_FUNC inline const Block<MatrixType,int(Mode)==int(Projective)?HDim:Dim,Dim> linearExt() const
  { return m_matrix.template block<int(Mode)==int(Projective)?HDim:Dim,Dim>(0,0); }

  EIGEN_DEVICE_FUNC inline Block<MatrixType,int(Mode)==int(Projective)?HDim:Dim,1> translationExt()
  { return m_matrix.template block<int(Mode)==int(Projective)?HDim:Dim,1>(0,Dim); }

  EIGEN_DEVICE_FUNC inline const Block<MatrixType,int(Mode)==int(Projective)?HDim:Dim,1> translationExt() const
  { return m_matrix.template block<int(Mode)==int(Projective)?HDim:Dim,1>(0,Dim); }

  #ifdef EIGEN_TRANSFORM_PLUGIN
  #include EIGEN_TRANSFORM_PLUGIN
  #endif

protected:
  #ifndef EIGEN_PARSED_BY_DOXYGEN
    EIGEN_DEVICE_FUNC static EIGEN_STRONG_INLINE void check_template_params()
    {
      EIGEN_STATIC_ASSERT((Options & (DontAlign|RowMajor)) == Options, INVALID_MATRIX_TEMPLATE_PARAMETERS)
    }
  #endif

};

typedef Transform<float,2,Isometry> Isometry2f;

typedef Transform<float,3,Isometry> Isometry3f;

typedef Transform<double,2,Isometry> Isometry2d;

typedef Transform<double,3,Isometry> Isometry3d;

typedef Transform<float,2,Affine> Affine2f;

typedef Transform<float,3,Affine> Affine3f;

typedef Transform<double,2,Affine> Affine2d;

typedef Transform<double,3,Affine> Affine3d;

typedef Transform<float,2,AffineCompact> AffineCompact2f;

typedef Transform<float,3,AffineCompact> AffineCompact3f;

typedef Transform<double,2,AffineCompact> AffineCompact2d;

typedef Transform<double,3,AffineCompact> AffineCompact3d;

typedef Transform<float,2,Projective> Projective2f;

typedef Transform<float,3,Projective> Projective3f;

typedef Transform<double,2,Projective> Projective2d;

typedef Transform<double,3,Projective> Projective3d;

#ifdef EIGEN_QT_SUPPORT

template<typename Scalar, int Dim, int Mode,int Options>
Transform<Scalar,Dim,Mode,Options>::Transform(const QMatrix& other)
{
  check_template_params();
  *this = other;
}

template<typename Scalar, int Dim, int Mode,int Options>
Transform<Scalar,Dim,Mode,Options>& Transform<Scalar,Dim,Mode,Options>::operator=(const QMatrix& other)
{
  EIGEN_STATIC_ASSERT(Dim==2, YOU_MADE_A_PROGRAMMING_MISTAKE)
  if (EIGEN_CONST_CONDITIONAL(Mode == int(AffineCompact)))
    m_matrix << other.m11(), other.m21(), other.dx(),
                other.m12(), other.m22(), other.dy();
  else
    m_matrix << other.m11(), other.m21(), other.dx(),
                other.m12(), other.m22(), other.dy(),
                0, 0, 1;
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
QMatrix Transform<Scalar,Dim,Mode,Options>::toQMatrix(void) const
{
  check_template_params();
  EIGEN_STATIC_ASSERT(Dim==2, YOU_MADE_A_PROGRAMMING_MISTAKE)
  return QMatrix(m_matrix.coeff(0,0), m_matrix.coeff(1,0),
                 m_matrix.coeff(0,1), m_matrix.coeff(1,1),
                 m_matrix.coeff(0,2), m_matrix.coeff(1,2));
}

template<typename Scalar, int Dim, int Mode,int Options>
Transform<Scalar,Dim,Mode,Options>::Transform(const QTransform& other)
{
  check_template_params();
  *this = other;
}

template<typename Scalar, int Dim, int Mode, int Options>
Transform<Scalar,Dim,Mode,Options>& Transform<Scalar,Dim,Mode,Options>::operator=(const QTransform& other)
{
  check_template_params();
  EIGEN_STATIC_ASSERT(Dim==2, YOU_MADE_A_PROGRAMMING_MISTAKE)
  if (EIGEN_CONST_CONDITIONAL(Mode == int(AffineCompact)))
    m_matrix << other.m11(), other.m21(), other.dx(),
                other.m12(), other.m22(), other.dy();
  else
    m_matrix << other.m11(), other.m21(), other.dx(),
                other.m12(), other.m22(), other.dy(),
                other.m13(), other.m23(), other.m33();
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
QTransform Transform<Scalar,Dim,Mode,Options>::toQTransform(void) const
{
  EIGEN_STATIC_ASSERT(Dim==2, YOU_MADE_A_PROGRAMMING_MISTAKE)
  if (EIGEN_CONST_CONDITIONAL(Mode == int(AffineCompact)))
    return QTransform(m_matrix.coeff(0,0), m_matrix.coeff(1,0),
                      m_matrix.coeff(0,1), m_matrix.coeff(1,1),
                      m_matrix.coeff(0,2), m_matrix.coeff(1,2));
  else
    return QTransform(m_matrix.coeff(0,0), m_matrix.coeff(1,0), m_matrix.coeff(2,0),
                      m_matrix.coeff(0,1), m_matrix.coeff(1,1), m_matrix.coeff(2,1),
                      m_matrix.coeff(0,2), m_matrix.coeff(1,2), m_matrix.coeff(2,2));
}
#endif

template<typename Scalar, int Dim, int Mode, int Options>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC Transform<Scalar,Dim,Mode,Options>&
Transform<Scalar,Dim,Mode,Options>::scale(const MatrixBase<OtherDerived> &other)
{
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(OtherDerived,int(Dim))
  EIGEN_STATIC_ASSERT(Mode!=int(Isometry), THIS_METHOD_IS_ONLY_FOR_SPECIFIC_TRANSFORMATIONS)
  linearExt().noalias() = (linearExt() * other.asDiagonal());
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
EIGEN_DEVICE_FUNC inline Transform<Scalar,Dim,Mode,Options>& Transform<Scalar,Dim,Mode,Options>::scale(const Scalar& s)
{
  EIGEN_STATIC_ASSERT(Mode!=int(Isometry), THIS_METHOD_IS_ONLY_FOR_SPECIFIC_TRANSFORMATIONS)
  linearExt() *= s;
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC Transform<Scalar,Dim,Mode,Options>&
Transform<Scalar,Dim,Mode,Options>::prescale(const MatrixBase<OtherDerived> &other)
{
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(OtherDerived,int(Dim))
  EIGEN_STATIC_ASSERT(Mode!=int(Isometry), THIS_METHOD_IS_ONLY_FOR_SPECIFIC_TRANSFORMATIONS)
  affine().noalias() = (other.asDiagonal() * affine());
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
EIGEN_DEVICE_FUNC inline Transform<Scalar,Dim,Mode,Options>& Transform<Scalar,Dim,Mode,Options>::prescale(const Scalar& s)
{
  EIGEN_STATIC_ASSERT(Mode!=int(Isometry), THIS_METHOD_IS_ONLY_FOR_SPECIFIC_TRANSFORMATIONS)
  m_matrix.template topRows<Dim>() *= s;
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC Transform<Scalar,Dim,Mode,Options>&
Transform<Scalar,Dim,Mode,Options>::translate(const MatrixBase<OtherDerived> &other)
{
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(OtherDerived,int(Dim))
  translationExt() += linearExt() * other;
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC Transform<Scalar,Dim,Mode,Options>&
Transform<Scalar,Dim,Mode,Options>::pretranslate(const MatrixBase<OtherDerived> &other)
{
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(OtherDerived,int(Dim))
  if(EIGEN_CONST_CONDITIONAL(int(Mode)==int(Projective)))
    affine() += other * m_matrix.row(Dim);
  else
    translation() += other;
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
template<typename RotationType>
EIGEN_DEVICE_FUNC Transform<Scalar,Dim,Mode,Options>&
Transform<Scalar,Dim,Mode,Options>::rotate(const RotationType& rotation)
{
  linearExt() *= internal::toRotationMatrix<Scalar,Dim>(rotation);
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
template<typename RotationType>
EIGEN_DEVICE_FUNC Transform<Scalar,Dim,Mode,Options>&
Transform<Scalar,Dim,Mode,Options>::prerotate(const RotationType& rotation)
{
  m_matrix.template block<Dim,HDim>(0,0) = internal::toRotationMatrix<Scalar,Dim>(rotation)
                                         * m_matrix.template block<Dim,HDim>(0,0);
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
EIGEN_DEVICE_FUNC Transform<Scalar,Dim,Mode,Options>&
Transform<Scalar,Dim,Mode,Options>::shear(const Scalar& sx, const Scalar& sy)
{
  EIGEN_STATIC_ASSERT(int(Dim)==2, YOU_MADE_A_PROGRAMMING_MISTAKE)
  EIGEN_STATIC_ASSERT(Mode!=int(Isometry), THIS_METHOD_IS_ONLY_FOR_SPECIFIC_TRANSFORMATIONS)
  VectorType tmp = linear().col(0)*sy + linear().col(1);
  linear() << linear().col(0) + linear().col(1)*sx, tmp;
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
EIGEN_DEVICE_FUNC Transform<Scalar,Dim,Mode,Options>&
Transform<Scalar,Dim,Mode,Options>::preshear(const Scalar& sx, const Scalar& sy)
{
  EIGEN_STATIC_ASSERT(int(Dim)==2, YOU_MADE_A_PROGRAMMING_MISTAKE)
  EIGEN_STATIC_ASSERT(Mode!=int(Isometry), THIS_METHOD_IS_ONLY_FOR_SPECIFIC_TRANSFORMATIONS)
  m_matrix.template block<Dim,HDim>(0,0) = LinearMatrixType(1, sx, sy, 1) * m_matrix.template block<Dim,HDim>(0,0);
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
EIGEN_DEVICE_FUNC inline Transform<Scalar,Dim,Mode,Options>& Transform<Scalar,Dim,Mode,Options>::operator=(const TranslationType& t)
{
  linear().setIdentity();
  translation() = t.vector();
  makeAffine();
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
EIGEN_DEVICE_FUNC inline Transform<Scalar,Dim,Mode,Options> Transform<Scalar,Dim,Mode,Options>::operator*(const TranslationType& t) const
{
  Transform res = *this;
  res.translate(t.vector());
  return res;
}

template<typename Scalar, int Dim, int Mode, int Options>
EIGEN_DEVICE_FUNC inline Transform<Scalar,Dim,Mode,Options>& Transform<Scalar,Dim,Mode,Options>::operator=(const UniformScaling<Scalar>& s)
{
  m_matrix.setZero();
  linear().diagonal().fill(s.factor());
  makeAffine();
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
template<typename Derived>
EIGEN_DEVICE_FUNC inline Transform<Scalar,Dim,Mode,Options>& Transform<Scalar,Dim,Mode,Options>::operator=(const RotationBase<Derived,Dim>& r)
{
  linear() = internal::toRotationMatrix<Scalar,Dim>(r);
  translation().setZero();
  makeAffine();
  return *this;
}

template<typename Scalar, int Dim, int Mode, int Options>
template<typename Derived>
EIGEN_DEVICE_FUNC inline Transform<Scalar,Dim,Mode,Options> Transform<Scalar,Dim,Mode,Options>::operator*(const RotationBase<Derived,Dim>& r) const
{
  Transform res = *this;
  res.rotate(r.derived());
  return res;
}

namespace internal {
template<int Mode> struct transform_rotation_impl {
  template<typename TransformType>
  EIGEN_DEVICE_FUNC static inline
  const typename TransformType::LinearMatrixType run(const TransformType& t)
  {
    typedef typename TransformType::LinearMatrixType LinearMatrixType;
    LinearMatrixType result;
    t.computeRotationScaling(&result, (LinearMatrixType*)0);
    return result;
  }
};
template<> struct transform_rotation_impl<Isometry> {
  template<typename TransformType>
  EIGEN_DEVICE_FUNC static inline
  typename TransformType::ConstLinearPart run(const TransformType& t)
  {
    return t.linear();
  }
};
}

template<typename Scalar, int Dim, int Mode, int Options>
EIGEN_DEVICE_FUNC
typename Transform<Scalar,Dim,Mode,Options>::RotationReturnType
Transform<Scalar,Dim,Mode,Options>::rotation() const
{
  return internal::transform_rotation_impl<Mode>::run(*this);
}

template<typename Scalar, int Dim, int Mode, int Options>
template<typename RotationMatrixType, typename ScalingMatrixType>
EIGEN_DEVICE_FUNC void Transform<Scalar,Dim,Mode,Options>::computeRotationScaling(RotationMatrixType *rotation, ScalingMatrixType *scaling) const
{
  
  JacobiSVD<LinearMatrixType> svd(linear(), ComputeFullU | ComputeFullV);

  Scalar x = (svd.matrixU() * svd.matrixV().adjoint()).determinant() < Scalar(0) ? Scalar(-1) : Scalar(1); 
  VectorType sv(svd.singularValues());
  sv.coeffRef(Dim-1) *= x;
  if(scaling) *scaling = svd.matrixV() * sv.asDiagonal() * svd.matrixV().adjoint();
  if(rotation)
  {
    LinearMatrixType m(svd.matrixU());
    m.col(Dim-1) *= x;
    *rotation = m * svd.matrixV().adjoint();
  }
}

template<typename Scalar, int Dim, int Mode, int Options>
template<typename ScalingMatrixType, typename RotationMatrixType>
EIGEN_DEVICE_FUNC void Transform<Scalar,Dim,Mode,Options>::computeScalingRotation(ScalingMatrixType *scaling, RotationMatrixType *rotation) const
{
  
  JacobiSVD<LinearMatrixType> svd(linear(), ComputeFullU | ComputeFullV);

  Scalar x = (svd.matrixU() * svd.matrixV().adjoint()).determinant() < Scalar(0) ? Scalar(-1) : Scalar(1); 
  VectorType sv(svd.singularValues());
  sv.coeffRef(Dim-1) *= x;
  if(scaling) *scaling = svd.matrixU() * sv.asDiagonal() * svd.matrixU().adjoint();
  if(rotation)
  {
    LinearMatrixType m(svd.matrixU());
    m.col(Dim-1) *= x;
    *rotation = m * svd.matrixV().adjoint();
  }
}

template<typename Scalar, int Dim, int Mode, int Options>
template<typename PositionDerived, typename OrientationType, typename ScaleDerived>
EIGEN_DEVICE_FUNC Transform<Scalar,Dim,Mode,Options>&
Transform<Scalar,Dim,Mode,Options>::fromPositionOrientationScale(const MatrixBase<PositionDerived> &position,
  const OrientationType& orientation, const MatrixBase<ScaleDerived> &scale)
{
  linear() = internal::toRotationMatrix<Scalar,Dim>(orientation);
  linear() *= scale.asDiagonal();
  translation() = position;
  makeAffine();
  return *this;
}

namespace internal {

template<int Mode>
struct transform_make_affine
{
  template<typename MatrixType>
  EIGEN_DEVICE_FUNC static void run(MatrixType &mat)
  {
    static const int Dim = MatrixType::ColsAtCompileTime-1;
    mat.template block<1,Dim>(Dim,0).setZero();
    mat.coeffRef(Dim,Dim) = typename MatrixType::Scalar(1);
  }
};

template<>
struct transform_make_affine<AffineCompact>
{
  template<typename MatrixType> EIGEN_DEVICE_FUNC static void run(MatrixType &) { }
};

template<typename TransformType, int Mode=TransformType::Mode>
struct projective_transform_inverse
{
  EIGEN_DEVICE_FUNC static inline void run(const TransformType&, TransformType&)
  {}
};

template<typename TransformType>
struct projective_transform_inverse<TransformType, Projective>
{
  EIGEN_DEVICE_FUNC static inline void run(const TransformType& m, TransformType& res)
  {
    res.matrix() = m.matrix().inverse();
  }
};

} 

template<typename Scalar, int Dim, int Mode, int Options>
EIGEN_DEVICE_FUNC Transform<Scalar,Dim,Mode,Options>
Transform<Scalar,Dim,Mode,Options>::inverse(TransformTraits hint) const
{
  Transform res;
  if (hint == Projective)
  {
    internal::projective_transform_inverse<Transform>::run(*this, res);
  }
  else
  {
    if (hint == Isometry)
    {
      res.matrix().template topLeftCorner<Dim,Dim>() = linear().transpose();
    }
    else if(hint&Affine)
    {
      res.matrix().template topLeftCorner<Dim,Dim>() = linear().inverse();
    }
    else
    {
      eigen_assert(false && "Invalid transform traits in Transform::Inverse");
    }
    
    res.matrix().template topRightCorner<Dim,1>()
      = - res.matrix().template topLeftCorner<Dim,Dim>() * translation();
    res.makeAffine(); 
  }
  return res;
}

namespace internal {

template<typename TransformType> struct transform_take_affine_part {
  typedef typename TransformType::MatrixType MatrixType;
  typedef typename TransformType::AffinePart AffinePart;
  typedef typename TransformType::ConstAffinePart ConstAffinePart;
  static inline AffinePart run(MatrixType& m)
  { return m.template block<TransformType::Dim,TransformType::HDim>(0,0); }
  static inline ConstAffinePart run(const MatrixType& m)
  { return m.template block<TransformType::Dim,TransformType::HDim>(0,0); }
};

template<typename Scalar, int Dim, int Options>
struct transform_take_affine_part<Transform<Scalar,Dim,AffineCompact, Options> > {
  typedef typename Transform<Scalar,Dim,AffineCompact,Options>::MatrixType MatrixType;
  static inline MatrixType& run(MatrixType& m) { return m; }
  static inline const MatrixType& run(const MatrixType& m) { return m; }
};

template<typename Other, int Mode, int Options, int Dim, int HDim>
struct transform_construct_from_matrix<Other, Mode,Options,Dim,HDim, Dim,Dim>
{
  static inline void run(Transform<typename Other::Scalar,Dim,Mode,Options> *transform, const Other& other)
  {
    transform->linear() = other;
    transform->translation().setZero();
    transform->makeAffine();
  }
};

template<typename Other, int Mode, int Options, int Dim, int HDim>
struct transform_construct_from_matrix<Other, Mode,Options,Dim,HDim, Dim,HDim>
{
  static inline void run(Transform<typename Other::Scalar,Dim,Mode,Options> *transform, const Other& other)
  {
    transform->affine() = other;
    transform->makeAffine();
  }
};

template<typename Other, int Mode, int Options, int Dim, int HDim>
struct transform_construct_from_matrix<Other, Mode,Options,Dim,HDim, HDim,HDim>
{
  static inline void run(Transform<typename Other::Scalar,Dim,Mode,Options> *transform, const Other& other)
  { transform->matrix() = other; }
};

template<typename Other, int Options, int Dim, int HDim>
struct transform_construct_from_matrix<Other, AffineCompact,Options,Dim,HDim, HDim,HDim>
{
  static inline void run(Transform<typename Other::Scalar,Dim,AffineCompact,Options> *transform, const Other& other)
  { transform->matrix() = other.template block<Dim,HDim>(0,0); }
};

template<int LhsMode,int RhsMode>
struct transform_product_result
{
  enum
  {
    Mode =
      (LhsMode == (int)Projective    || RhsMode == (int)Projective    ) ? Projective :
      (LhsMode == (int)Affine        || RhsMode == (int)Affine        ) ? Affine :
      (LhsMode == (int)AffineCompact || RhsMode == (int)AffineCompact ) ? AffineCompact :
      (LhsMode == (int)Isometry      || RhsMode == (int)Isometry      ) ? Isometry : Projective
  };
};

template< typename TransformType, typename MatrixType, int RhsCols>
struct transform_right_product_impl< TransformType, MatrixType, 0, RhsCols>
{
  typedef typename MatrixType::PlainObject ResultType;

  static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ResultType run(const TransformType& T, const MatrixType& other)
  {
    return T.matrix() * other;
  }
};

template< typename TransformType, typename MatrixType, int RhsCols>
struct transform_right_product_impl< TransformType, MatrixType, 1, RhsCols>
{
  enum {
    Dim = TransformType::Dim,
    HDim = TransformType::HDim,
    OtherRows = MatrixType::RowsAtCompileTime,
    OtherCols = MatrixType::ColsAtCompileTime
  };

  typedef typename MatrixType::PlainObject ResultType;

  static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ResultType run(const TransformType& T, const MatrixType& other)
  {
    EIGEN_STATIC_ASSERT(OtherRows==HDim, YOU_MIXED_MATRICES_OF_DIFFERENT_SIZES);

    typedef Block<ResultType, Dim, OtherCols, int(MatrixType::RowsAtCompileTime)==Dim> TopLeftLhs;

    ResultType res(other.rows(),other.cols());
    TopLeftLhs(res, 0, 0, Dim, other.cols()).noalias() = T.affine() * other;
    res.row(OtherRows-1) = other.row(OtherRows-1);

    return res;
  }
};

template< typename TransformType, typename MatrixType, int RhsCols>
struct transform_right_product_impl< TransformType, MatrixType, 2, RhsCols>
{
  enum {
    Dim = TransformType::Dim,
    HDim = TransformType::HDim,
    OtherRows = MatrixType::RowsAtCompileTime,
    OtherCols = MatrixType::ColsAtCompileTime
  };

  typedef typename MatrixType::PlainObject ResultType;

  static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ResultType run(const TransformType& T, const MatrixType& other)
  {
    EIGEN_STATIC_ASSERT(OtherRows==Dim, YOU_MIXED_MATRICES_OF_DIFFERENT_SIZES);

    typedef Block<ResultType, Dim, OtherCols, true> TopLeftLhs;
    ResultType res(Replicate<typename TransformType::ConstTranslationPart, 1, OtherCols>(T.translation(),1,other.cols()));
    TopLeftLhs(res, 0, 0, Dim, other.cols()).noalias() += T.linear() * other;

    return res;
  }
};

template< typename TransformType, typename MatrixType >
struct transform_right_product_impl< TransformType, MatrixType, 2, 1> 
{
  typedef typename TransformType::MatrixType TransformMatrix;
  enum {
    Dim = TransformType::Dim,
    HDim = TransformType::HDim,
    OtherRows = MatrixType::RowsAtCompileTime,
    WorkingRows = EIGEN_PLAIN_ENUM_MIN(TransformMatrix::RowsAtCompileTime,HDim)
  };

  typedef typename MatrixType::PlainObject ResultType;

  static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ResultType run(const TransformType& T, const MatrixType& other)
  {
    EIGEN_STATIC_ASSERT(OtherRows==Dim, YOU_MIXED_MATRICES_OF_DIFFERENT_SIZES);

    Matrix<typename ResultType::Scalar, Dim+1, 1> rhs;
    rhs.template head<Dim>() = other; rhs[Dim] = typename ResultType::Scalar(1);
    Matrix<typename ResultType::Scalar, WorkingRows, 1> res(T.matrix() * rhs);
    return res.template head<Dim>();
  }
};

template<typename Other,int Mode, int Options, int Dim, int HDim>
struct transform_left_product_impl<Other,Mode,Options,Dim,HDim, HDim,HDim>
{
  typedef Transform<typename Other::Scalar,Dim,Mode,Options> TransformType;
  typedef typename TransformType::MatrixType MatrixType;
  typedef Transform<typename Other::Scalar,Dim,Projective,Options> ResultType;
  static ResultType run(const Other& other,const TransformType& tr)
  { return ResultType(other * tr.matrix()); }
};

template<typename Other, int Options, int Dim, int HDim>
struct transform_left_product_impl<Other,AffineCompact,Options,Dim,HDim, HDim,HDim>
{
  typedef Transform<typename Other::Scalar,Dim,AffineCompact,Options> TransformType;
  typedef typename TransformType::MatrixType MatrixType;
  typedef Transform<typename Other::Scalar,Dim,Projective,Options> ResultType;
  static ResultType run(const Other& other,const TransformType& tr)
  {
    ResultType res;
    res.matrix().noalias() = other.template block<HDim,Dim>(0,0) * tr.matrix();
    res.matrix().col(Dim) += other.col(Dim);
    return res;
  }
};

template<typename Other,int Mode, int Options, int Dim, int HDim>
struct transform_left_product_impl<Other,Mode,Options,Dim,HDim, Dim,HDim>
{
  typedef Transform<typename Other::Scalar,Dim,Mode,Options> TransformType;
  typedef typename TransformType::MatrixType MatrixType;
  typedef TransformType ResultType;
  static ResultType run(const Other& other,const TransformType& tr)
  {
    ResultType res;
    res.affine().noalias() = other * tr.matrix();
    res.matrix().row(Dim) = tr.matrix().row(Dim);
    return res;
  }
};

template<typename Other, int Options, int Dim, int HDim>
struct transform_left_product_impl<Other,AffineCompact,Options,Dim,HDim, Dim,HDim>
{
  typedef Transform<typename Other::Scalar,Dim,AffineCompact,Options> TransformType;
  typedef typename TransformType::MatrixType MatrixType;
  typedef TransformType ResultType;
  static ResultType run(const Other& other,const TransformType& tr)
  {
    ResultType res;
    res.matrix().noalias() = other.template block<Dim,Dim>(0,0) * tr.matrix();
    res.translation() += other.col(Dim);
    return res;
  }
};

template<typename Other,int Mode, int Options, int Dim, int HDim>
struct transform_left_product_impl<Other,Mode,Options,Dim,HDim, Dim,Dim>
{
  typedef Transform<typename Other::Scalar,Dim,Mode,Options> TransformType;
  typedef typename TransformType::MatrixType MatrixType;
  typedef TransformType ResultType;
  static ResultType run(const Other& other, const TransformType& tr)
  {
    TransformType res;
    if(Mode!=int(AffineCompact))
      res.matrix().row(Dim) = tr.matrix().row(Dim);
    res.matrix().template topRows<Dim>().noalias()
      = other * tr.matrix().template topRows<Dim>();
    return res;
  }
};

template<typename Scalar, int Dim, int LhsMode, int LhsOptions, int RhsMode, int RhsOptions>
struct transform_transform_product_impl<Transform<Scalar,Dim,LhsMode,LhsOptions>,Transform<Scalar,Dim,RhsMode,RhsOptions>,false >
{
  enum { ResultMode = transform_product_result<LhsMode,RhsMode>::Mode };
  typedef Transform<Scalar,Dim,LhsMode,LhsOptions> Lhs;
  typedef Transform<Scalar,Dim,RhsMode,RhsOptions> Rhs;
  typedef Transform<Scalar,Dim,ResultMode,LhsOptions> ResultType;
  static ResultType run(const Lhs& lhs, const Rhs& rhs)
  {
    ResultType res;
    res.linear() = lhs.linear() * rhs.linear();
    res.translation() = lhs.linear() * rhs.translation() + lhs.translation();
    res.makeAffine();
    return res;
  }
};

template<typename Scalar, int Dim, int LhsMode, int LhsOptions, int RhsMode, int RhsOptions>
struct transform_transform_product_impl<Transform<Scalar,Dim,LhsMode,LhsOptions>,Transform<Scalar,Dim,RhsMode,RhsOptions>,true >
{
  typedef Transform<Scalar,Dim,LhsMode,LhsOptions> Lhs;
  typedef Transform<Scalar,Dim,RhsMode,RhsOptions> Rhs;
  typedef Transform<Scalar,Dim,Projective> ResultType;
  static ResultType run(const Lhs& lhs, const Rhs& rhs)
  {
    return ResultType( lhs.matrix() * rhs.matrix() );
  }
};

template<typename Scalar, int Dim, int LhsOptions, int RhsOptions>
struct transform_transform_product_impl<Transform<Scalar,Dim,AffineCompact,LhsOptions>,Transform<Scalar,Dim,Projective,RhsOptions>,true >
{
  typedef Transform<Scalar,Dim,AffineCompact,LhsOptions> Lhs;
  typedef Transform<Scalar,Dim,Projective,RhsOptions> Rhs;
  typedef Transform<Scalar,Dim,Projective> ResultType;
  static ResultType run(const Lhs& lhs, const Rhs& rhs)
  {
    ResultType res;
    res.matrix().template topRows<Dim>() = lhs.matrix() * rhs.matrix();
    res.matrix().row(Dim) = rhs.matrix().row(Dim);
    return res;
  }
};

template<typename Scalar, int Dim, int LhsOptions, int RhsOptions>
struct transform_transform_product_impl<Transform<Scalar,Dim,Projective,LhsOptions>,Transform<Scalar,Dim,AffineCompact,RhsOptions>,true >
{
  typedef Transform<Scalar,Dim,Projective,LhsOptions> Lhs;
  typedef Transform<Scalar,Dim,AffineCompact,RhsOptions> Rhs;
  typedef Transform<Scalar,Dim,Projective> ResultType;
  static ResultType run(const Lhs& lhs, const Rhs& rhs)
  {
    ResultType res(lhs.matrix().template leftCols<Dim>() * rhs.matrix());
    res.matrix().col(Dim) += lhs.matrix().col(Dim);
    return res;
  }
};

} 

} 

#endif 
