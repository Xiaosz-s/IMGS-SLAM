template<typename OtherDerived>
EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const EIGEN_CWISE_BINARY_RETURN_TYPE(Derived,OtherDerived,product)
operator*(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived> &other) const
{
  return EIGEN_CWISE_BINARY_RETURN_TYPE(Derived,OtherDerived,product)(derived(), other.derived());
}

template<typename OtherDerived>
EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_quotient_op<Scalar,typename OtherDerived::Scalar>, const Derived, const OtherDerived>
operator/(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived> &other) const
{
  return CwiseBinaryOp<internal::scalar_quotient_op<Scalar,typename OtherDerived::Scalar>, const Derived, const OtherDerived>(derived(), other.derived());
}

EIGEN_MAKE_CWISE_BINARY_OP(min,min)

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_min_op<Scalar,Scalar>, const Derived,
                                        const CwiseNullaryOp<internal::scalar_constant_op<Scalar>, PlainObject> >
#ifdef EIGEN_PARSED_BY_DOXYGEN
min
#else
(min)
#endif
(const Scalar &other) const
{
  return (min)(Derived::PlainObject::Constant(rows(), cols(), other));
}

EIGEN_MAKE_CWISE_BINARY_OP(max,max)

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_max_op<Scalar,Scalar>, const Derived,
                                        const CwiseNullaryOp<internal::scalar_constant_op<Scalar>, PlainObject> >
#ifdef EIGEN_PARSED_BY_DOXYGEN
max
#else
(max)
#endif
(const Scalar &other) const
{
  return (max)(Derived::PlainObject::Constant(rows(), cols(), other));
}

EIGEN_MAKE_CWISE_BINARY_OP(absolute_difference,absolute_difference)

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_absolute_difference_op<Scalar,Scalar>, const Derived,
                                        const CwiseNullaryOp<internal::scalar_constant_op<Scalar>, PlainObject> >
#ifdef EIGEN_PARSED_BY_DOXYGEN
absolute_difference
#else
(absolute_difference)
#endif
(const Scalar &other) const
{
  return (absolute_difference)(Derived::PlainObject::Constant(rows(), cols(), other));
}

EIGEN_MAKE_CWISE_BINARY_OP(pow,pow)

#ifndef EIGEN_PARSED_BY_DOXYGEN
EIGEN_MAKE_SCALAR_BINARY_OP_ONTHERIGHT(pow,pow)
#else

template<typename T>
const CwiseBinaryOp<internal::scalar_pow_op<Scalar,T>,Derived,Constant<T> > pow(const T& exponent) const;
#endif

#define EIGEN_MAKE_CWISE_COMP_OP(OP, COMPARATOR) \
template<typename OtherDerived> \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_cmp_op<Scalar, typename OtherDerived::Scalar, internal::cmp_ ## COMPARATOR>, const Derived, const OtherDerived> \
OP(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived> &other) const \
{ \
  return CwiseBinaryOp<internal::scalar_cmp_op<Scalar, typename OtherDerived::Scalar, internal::cmp_ ## COMPARATOR>, const Derived, const OtherDerived>(derived(), other.derived()); \
}\
typedef CwiseBinaryOp<internal::scalar_cmp_op<Scalar,Scalar, internal::cmp_ ## COMPARATOR>, const Derived, const CwiseNullaryOp<internal::scalar_constant_op<Scalar>, PlainObject> > Cmp ## COMPARATOR ## ReturnType; \
typedef CwiseBinaryOp<internal::scalar_cmp_op<Scalar,Scalar, internal::cmp_ ## COMPARATOR>, const CwiseNullaryOp<internal::scalar_constant_op<Scalar>, PlainObject>, const Derived > RCmp ## COMPARATOR ## ReturnType; \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Cmp ## COMPARATOR ## ReturnType \
OP(const Scalar& s) const { \
  return this->OP(Derived::PlainObject::Constant(rows(), cols(), s)); \
} \
EIGEN_DEVICE_FUNC friend EIGEN_STRONG_INLINE const RCmp ## COMPARATOR ## ReturnType \
OP(const Scalar& s, const EIGEN_CURRENT_STORAGE_BASE_CLASS<Derived>& d) { \
  return Derived::PlainObject::Constant(d.rows(), d.cols(), s).OP(d); \
}

#define EIGEN_MAKE_CWISE_COMP_R_OP(OP, R_OP, RCOMPARATOR) \
template<typename OtherDerived> \
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_cmp_op<typename OtherDerived::Scalar, Scalar, internal::cmp_##RCOMPARATOR>, const OtherDerived, const Derived> \
OP(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived> &other) const \
{ \
  return CwiseBinaryOp<internal::scalar_cmp_op<typename OtherDerived::Scalar, Scalar, internal::cmp_##RCOMPARATOR>, const OtherDerived, const Derived>(other.derived(), derived()); \
} \
EIGEN_DEVICE_FUNC \
inline const RCmp ## RCOMPARATOR ## ReturnType \
OP(const Scalar& s) const { \
  return Derived::PlainObject::Constant(rows(), cols(), s).R_OP(*this); \
} \
friend inline const Cmp ## RCOMPARATOR ## ReturnType \
OP(const Scalar& s, const Derived& d) { \
  return d.R_OP(Derived::PlainObject::Constant(d.rows(), d.cols(), s)); \
}

EIGEN_MAKE_CWISE_COMP_OP(operator<, LT)

EIGEN_MAKE_CWISE_COMP_OP(operator<=, LE)

EIGEN_MAKE_CWISE_COMP_R_OP(operator>, operator<, LT)

EIGEN_MAKE_CWISE_COMP_R_OP(operator>=, operator<=, LE)

EIGEN_MAKE_CWISE_COMP_OP(operator==, EQ)

EIGEN_MAKE_CWISE_COMP_OP(operator!=, NEQ)

#undef EIGEN_MAKE_CWISE_COMP_OP
#undef EIGEN_MAKE_CWISE_COMP_R_OP

#ifndef EIGEN_PARSED_BY_DOXYGEN
EIGEN_MAKE_SCALAR_BINARY_OP(operator+,sum)
#else

template<typename T>
const CwiseBinaryOp<internal::scalar_sum_op<Scalar,T>,Derived,Constant<T> > operator+(const T& scalar) const;

template<typename T> friend
const CwiseBinaryOp<internal::scalar_sum_op<T,Scalar>,Constant<T>,Derived> operator+(const T& scalar, const StorageBaseType& expr);
#endif

#ifndef EIGEN_PARSED_BY_DOXYGEN
EIGEN_MAKE_SCALAR_BINARY_OP(operator-,difference)
#else

template<typename T>
const CwiseBinaryOp<internal::scalar_difference_op<Scalar,T>,Derived,Constant<T> > operator-(const T& scalar) const;

template<typename T> friend
const CwiseBinaryOp<internal::scalar_difference_op<T,Scalar>,Constant<T>,Derived> operator-(const T& scalar, const StorageBaseType& expr);
#endif

#ifndef EIGEN_PARSED_BY_DOXYGEN
  EIGEN_MAKE_SCALAR_BINARY_OP_ONTHELEFT(operator/,quotient)
#else

  template<typename T> friend
  inline const CwiseBinaryOp<internal::scalar_quotient_op<T,Scalar>,Constant<T>,Derived>
  operator/(const T& s,const StorageBaseType& a);
#endif

template<typename OtherDerived>
EIGEN_DEVICE_FUNC
inline const CwiseBinaryOp<internal::scalar_boolean_xor_op, const Derived, const OtherDerived>
operator^(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived> &other) const
{
  EIGEN_STATIC_ASSERT((internal::is_same<bool,Scalar>::value && internal::is_same<bool,typename OtherDerived::Scalar>::value),
                      THIS_METHOD_IS_ONLY_FOR_EXPRESSIONS_OF_BOOL);
  return CwiseBinaryOp<internal::scalar_boolean_xor_op, const Derived, const OtherDerived>(derived(),other.derived());
}

#if 0

template<typename DerivedN>
inline const CwiseBinaryOp<internal::scalar_polygamma_op<Scalar>, const DerivedN, const Derived>
polygamma(const EIGEN_CURRENT_STORAGE_BASE_CLASS<DerivedN> &n) const
{
  return CwiseBinaryOp<internal::scalar_polygamma_op<Scalar>, const DerivedN, const Derived>(n.derived(), this->derived());
}
#endif

template<typename DerivedQ>
inline const CwiseBinaryOp<internal::scalar_zeta_op<Scalar>, const Derived, const DerivedQ>
zeta(const EIGEN_CURRENT_STORAGE_BASE_CLASS<DerivedQ> &q) const
{
  return CwiseBinaryOp<internal::scalar_zeta_op<Scalar>, const Derived, const DerivedQ>(this->derived(), q.derived());
}
