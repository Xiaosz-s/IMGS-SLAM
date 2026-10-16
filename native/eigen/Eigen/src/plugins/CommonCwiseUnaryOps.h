#ifndef EIGEN_PARSED_BY_DOXYGEN

typedef typename internal::conditional<NumTraits<Scalar>::IsComplex,
                    const CwiseUnaryOp<internal::scalar_conjugate_op<Scalar>, const Derived>,
                    const Derived&
                  >::type ConjugateReturnType;

typedef typename internal::conditional<NumTraits<Scalar>::IsComplex,
                    const CwiseUnaryOp<internal::scalar_real_op<Scalar>, const Derived>,
                    const Derived&
                  >::type RealReturnType;

typedef typename internal::conditional<NumTraits<Scalar>::IsComplex,
                    CwiseUnaryView<internal::scalar_real_ref_op<Scalar>, Derived>,
                    Derived&
                  >::type NonConstRealReturnType;

typedef CwiseUnaryOp<internal::scalar_imag_op<Scalar>, const Derived> ImagReturnType;

typedef CwiseUnaryView<internal::scalar_imag_ref_op<Scalar>, Derived> NonConstImagReturnType;

typedef CwiseUnaryOp<internal::scalar_opposite_op<Scalar>, const Derived> NegativeReturnType;

#endif 

EIGEN_DOC_UNARY_ADDONS(operator-,opposite)

EIGEN_DEVICE_FUNC
inline const NegativeReturnType
operator-() const { return NegativeReturnType(derived()); }

template<class NewType> struct CastXpr { typedef typename internal::cast_return_type<Derived,const CwiseUnaryOp<internal::scalar_cast_op<Scalar, NewType>, const Derived> >::type Type; };

EIGEN_DOC_UNARY_ADDONS(cast,conversion function)

template<typename NewType>
EIGEN_DEVICE_FUNC
typename CastXpr<NewType>::Type
cast() const
{
  return typename CastXpr<NewType>::Type(derived());
}

EIGEN_DOC_UNARY_ADDONS(conjugate,complex conjugate)

EIGEN_DEVICE_FUNC
inline ConjugateReturnType
conjugate() const
{
  return ConjugateReturnType(derived());
}

EIGEN_DOC_UNARY_ADDONS(conjugate,complex conjugate)

template<bool Cond>
EIGEN_DEVICE_FUNC
inline typename internal::conditional<Cond,ConjugateReturnType,const Derived&>::type
conjugateIf() const
{
  typedef typename internal::conditional<Cond,ConjugateReturnType,const Derived&>::type ReturnType;
  return ReturnType(derived());
}

EIGEN_DOC_UNARY_ADDONS(real,real part function)

EIGEN_DEVICE_FUNC
inline RealReturnType
real() const { return RealReturnType(derived()); }

EIGEN_DOC_UNARY_ADDONS(imag,imaginary part function)

EIGEN_DEVICE_FUNC
inline const ImagReturnType
imag() const { return ImagReturnType(derived()); }

EIGEN_DOC_UNARY_ADDONS(unaryExpr,unary function)

template<typename CustomUnaryOp>
EIGEN_DEVICE_FUNC
inline const CwiseUnaryOp<CustomUnaryOp, const Derived>
unaryExpr(const CustomUnaryOp& func = CustomUnaryOp()) const
{
  return CwiseUnaryOp<CustomUnaryOp, const Derived>(derived(), func);
}

EIGEN_DOC_UNARY_ADDONS(unaryViewExpr,unary function)

template<typename CustomViewOp>
EIGEN_DEVICE_FUNC
inline const CwiseUnaryView<CustomViewOp, const Derived>
unaryViewExpr(const CustomViewOp& func = CustomViewOp()) const
{
  return CwiseUnaryView<CustomViewOp, const Derived>(derived(), func);
}

EIGEN_DOC_UNARY_ADDONS(real,real part function)

EIGEN_DEVICE_FUNC
inline NonConstRealReturnType
real() { return NonConstRealReturnType(derived()); }

EIGEN_DOC_UNARY_ADDONS(imag,imaginary part function)

EIGEN_DEVICE_FUNC
inline NonConstImagReturnType
imag() { return NonConstImagReturnType(derived()); }
