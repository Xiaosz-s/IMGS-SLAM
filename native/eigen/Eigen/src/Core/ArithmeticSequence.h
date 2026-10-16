#ifndef EIGEN_ARITHMETIC_SEQUENCE_H
#define EIGEN_ARITHMETIC_SEQUENCE_H

namespace Eigen {

namespace internal {

#if (!EIGEN_HAS_CXX11) || !((!EIGEN_COMP_GNUC) || EIGEN_COMP_GNUC>=48)
template<typename T> struct aseq_negate {};

template<> struct aseq_negate<Index> {
  typedef Index type;
};

template<int N> struct aseq_negate<FixedInt<N> > {
  typedef FixedInt<-N> type;
};

template<> struct aseq_negate<FixedInt<DynamicIndex> > {};

template<typename FirstType,typename SizeType,typename IncrType,
         bool FirstIsSymbolic=symbolic::is_symbolic<FirstType>::value,
         bool SizeIsSymbolic =symbolic::is_symbolic<SizeType>::value>
struct aseq_reverse_first_type {
  typedef Index type;
};

template<typename FirstType,typename SizeType,typename IncrType>
struct aseq_reverse_first_type<FirstType,SizeType,IncrType,true,true> {
  typedef symbolic::AddExpr<FirstType,
                            symbolic::ProductExpr<symbolic::AddExpr<SizeType,symbolic::ValueExpr<FixedInt<-1> > >,
                                                  symbolic::ValueExpr<IncrType> >
                           > type;
};

template<typename SizeType,typename IncrType,typename EnableIf = void>
struct aseq_reverse_first_type_aux {
  typedef Index type;
};

template<typename SizeType,typename IncrType>
struct aseq_reverse_first_type_aux<SizeType,IncrType,typename internal::enable_if<bool((SizeType::value+IncrType::value)|0x1)>::type> {
  typedef FixedInt<(SizeType::value-1)*IncrType::value> type;
};

template<typename FirstType,typename SizeType,typename IncrType>
struct aseq_reverse_first_type<FirstType,SizeType,IncrType,true,false> {
  typedef typename aseq_reverse_first_type_aux<SizeType,IncrType>::type Aux;
  typedef symbolic::AddExpr<FirstType,symbolic::ValueExpr<Aux> > type;
};

template<typename FirstType,typename SizeType,typename IncrType>
struct aseq_reverse_first_type<FirstType,SizeType,IncrType,false,true> {
  typedef symbolic::AddExpr<symbolic::ProductExpr<symbolic::AddExpr<SizeType,symbolic::ValueExpr<FixedInt<-1> > >,
                                                  symbolic::ValueExpr<IncrType> >,
                            symbolic::ValueExpr<> > type;
};
#endif

template<typename T> struct cleanup_seq_incr {
  typedef typename cleanup_index_type<T,DynamicIndex>::type type;
};

}

template<typename FirstType=Index,typename SizeType=Index,typename IncrType=internal::FixedInt<1> >
class ArithmeticSequence;

template<typename FirstType,typename SizeType,typename IncrType>
ArithmeticSequence<typename internal::cleanup_index_type<FirstType>::type,
                   typename internal::cleanup_index_type<SizeType>::type,
                   typename internal::cleanup_seq_incr<IncrType>::type >
seqN(FirstType first, SizeType size, IncrType incr);

template<typename FirstType,typename SizeType,typename IncrType>
class ArithmeticSequence
{
public:
  ArithmeticSequence(FirstType first, SizeType size) : m_first(first), m_size(size) {}
  ArithmeticSequence(FirstType first, SizeType size, IncrType incr) : m_first(first), m_size(size), m_incr(incr) {}

  enum {
    SizeAtCompileTime = internal::get_fixed_value<SizeType>::value,
    IncrAtCompileTime = internal::get_fixed_value<IncrType,DynamicIndex>::value
  };

  Index size()  const { return m_size; }

  Index first()  const { return m_first; }

  Index operator[](Index i) const { return m_first + i * m_incr; }

  const FirstType& firstObject() const { return m_first; }
  const SizeType&  sizeObject()  const { return m_size; }
  const IncrType&  incrObject()  const { return m_incr; }

protected:
  FirstType m_first;
  SizeType  m_size;
  IncrType  m_incr;

public:

#if EIGEN_HAS_CXX11 && ((!EIGEN_COMP_GNUC) || EIGEN_COMP_GNUC>=48)
  auto reverse() const -> decltype(Eigen::seqN(m_first+(m_size+fix<-1>())*m_incr,m_size,-m_incr)) {
    return seqN(m_first+(m_size+fix<-1>())*m_incr,m_size,-m_incr);
  }
#else
protected:
  typedef typename internal::aseq_negate<IncrType>::type ReverseIncrType;
  typedef typename internal::aseq_reverse_first_type<FirstType,SizeType,IncrType>::type ReverseFirstType;
public:
  ArithmeticSequence<ReverseFirstType,SizeType,ReverseIncrType>
  reverse() const {
    return seqN(m_first+(m_size+fix<-1>())*m_incr,m_size,-m_incr);
  }
#endif
};

template<typename FirstType,typename SizeType,typename IncrType>
ArithmeticSequence<typename internal::cleanup_index_type<FirstType>::type,typename internal::cleanup_index_type<SizeType>::type,typename internal::cleanup_seq_incr<IncrType>::type >
seqN(FirstType first, SizeType size, IncrType incr)  {
  return ArithmeticSequence<typename internal::cleanup_index_type<FirstType>::type,typename internal::cleanup_index_type<SizeType>::type,typename internal::cleanup_seq_incr<IncrType>::type>(first,size,incr);
}

template<typename FirstType,typename SizeType>
ArithmeticSequence<typename internal::cleanup_index_type<FirstType>::type,typename internal::cleanup_index_type<SizeType>::type >
seqN(FirstType first, SizeType size)  {
  return ArithmeticSequence<typename internal::cleanup_index_type<FirstType>::type,typename internal::cleanup_index_type<SizeType>::type>(first,size);
}

#ifdef EIGEN_PARSED_BY_DOXYGEN

template<typename FirstType,typename LastType, typename IncrType>
auto seq(FirstType f, LastType l, IncrType incr);

template<typename FirstType,typename LastType>
auto seq(FirstType f, LastType l);

#else 

#if EIGEN_HAS_CXX11
template<typename FirstType,typename LastType>
auto seq(FirstType f, LastType l) -> decltype(seqN(typename internal::cleanup_index_type<FirstType>::type(f),
                                                   (  typename internal::cleanup_index_type<LastType>::type(l)
                                                    - typename internal::cleanup_index_type<FirstType>::type(f)+fix<1>())))
{
  return seqN(typename internal::cleanup_index_type<FirstType>::type(f),
              (typename internal::cleanup_index_type<LastType>::type(l)
               -typename internal::cleanup_index_type<FirstType>::type(f)+fix<1>()));
}

template<typename FirstType,typename LastType, typename IncrType>
auto seq(FirstType f, LastType l, IncrType incr)
  -> decltype(seqN(typename internal::cleanup_index_type<FirstType>::type(f),
                   (   typename internal::cleanup_index_type<LastType>::type(l)
                     - typename internal::cleanup_index_type<FirstType>::type(f)+typename internal::cleanup_seq_incr<IncrType>::type(incr)
                   ) / typename internal::cleanup_seq_incr<IncrType>::type(incr),
                   typename internal::cleanup_seq_incr<IncrType>::type(incr)))
{
  typedef typename internal::cleanup_seq_incr<IncrType>::type CleanedIncrType;
  return seqN(typename internal::cleanup_index_type<FirstType>::type(f),
              ( typename internal::cleanup_index_type<LastType>::type(l)
               -typename internal::cleanup_index_type<FirstType>::type(f)+CleanedIncrType(incr)) / CleanedIncrType(incr),
              CleanedIncrType(incr));
}

#else 

template<typename FirstType,typename LastType>
typename internal::enable_if<!(symbolic::is_symbolic<FirstType>::value || symbolic::is_symbolic<LastType>::value),
                             ArithmeticSequence<typename internal::cleanup_index_type<FirstType>::type,Index> >::type
seq(FirstType f, LastType l)
{
  return seqN(typename internal::cleanup_index_type<FirstType>::type(f),
              Index((typename internal::cleanup_index_type<LastType>::type(l)-typename internal::cleanup_index_type<FirstType>::type(f)+fix<1>())));
}

template<typename FirstTypeDerived,typename LastType>
typename internal::enable_if<!symbolic::is_symbolic<LastType>::value,
    ArithmeticSequence<FirstTypeDerived, symbolic::AddExpr<symbolic::AddExpr<symbolic::NegateExpr<FirstTypeDerived>,symbolic::ValueExpr<> >,
                                                            symbolic::ValueExpr<internal::FixedInt<1> > > > >::type
seq(const symbolic::BaseExpr<FirstTypeDerived> &f, LastType l)
{
  return seqN(f.derived(),(typename internal::cleanup_index_type<LastType>::type(l)-f.derived()+fix<1>()));
}

template<typename FirstType,typename LastTypeDerived>
typename internal::enable_if<!symbolic::is_symbolic<FirstType>::value,
    ArithmeticSequence<typename internal::cleanup_index_type<FirstType>::type,
                        symbolic::AddExpr<symbolic::AddExpr<LastTypeDerived,symbolic::ValueExpr<> >,
                                          symbolic::ValueExpr<internal::FixedInt<1> > > > >::type
seq(FirstType f, const symbolic::BaseExpr<LastTypeDerived> &l)
{
  return seqN(typename internal::cleanup_index_type<FirstType>::type(f),(l.derived()-typename internal::cleanup_index_type<FirstType>::type(f)+fix<1>()));
}

template<typename FirstTypeDerived,typename LastTypeDerived>
ArithmeticSequence<FirstTypeDerived,
                    symbolic::AddExpr<symbolic::AddExpr<LastTypeDerived,symbolic::NegateExpr<FirstTypeDerived> >,symbolic::ValueExpr<internal::FixedInt<1> > > >
seq(const symbolic::BaseExpr<FirstTypeDerived> &f, const symbolic::BaseExpr<LastTypeDerived> &l)
{
  return seqN(f.derived(),(l.derived()-f.derived()+fix<1>()));
}

template<typename FirstType,typename LastType, typename IncrType>
typename internal::enable_if<!(symbolic::is_symbolic<FirstType>::value || symbolic::is_symbolic<LastType>::value),
    ArithmeticSequence<typename internal::cleanup_index_type<FirstType>::type,Index,typename internal::cleanup_seq_incr<IncrType>::type> >::type
seq(FirstType f, LastType l, IncrType incr)
{
  typedef typename internal::cleanup_seq_incr<IncrType>::type CleanedIncrType;
  return seqN(typename internal::cleanup_index_type<FirstType>::type(f),
              Index((typename internal::cleanup_index_type<LastType>::type(l)-typename internal::cleanup_index_type<FirstType>::type(f)+CleanedIncrType(incr))/CleanedIncrType(incr)), incr);
}

template<typename FirstTypeDerived,typename LastType, typename IncrType>
typename internal::enable_if<!symbolic::is_symbolic<LastType>::value,
    ArithmeticSequence<FirstTypeDerived,
                        symbolic::QuotientExpr<symbolic::AddExpr<symbolic::AddExpr<symbolic::NegateExpr<FirstTypeDerived>,
                                                                                   symbolic::ValueExpr<> >,
                                                                 symbolic::ValueExpr<typename internal::cleanup_seq_incr<IncrType>::type> >,
                                              symbolic::ValueExpr<typename internal::cleanup_seq_incr<IncrType>::type> >,
                        typename internal::cleanup_seq_incr<IncrType>::type> >::type
seq(const symbolic::BaseExpr<FirstTypeDerived> &f, LastType l, IncrType incr)
{
  typedef typename internal::cleanup_seq_incr<IncrType>::type CleanedIncrType;
  return seqN(f.derived(),(typename internal::cleanup_index_type<LastType>::type(l)-f.derived()+CleanedIncrType(incr))/CleanedIncrType(incr), incr);
}

template<typename FirstType,typename LastTypeDerived, typename IncrType>
typename internal::enable_if<!symbolic::is_symbolic<FirstType>::value,
    ArithmeticSequence<typename internal::cleanup_index_type<FirstType>::type,
                        symbolic::QuotientExpr<symbolic::AddExpr<symbolic::AddExpr<LastTypeDerived,symbolic::ValueExpr<> >,
                                                                 symbolic::ValueExpr<typename internal::cleanup_seq_incr<IncrType>::type> >,
                                               symbolic::ValueExpr<typename internal::cleanup_seq_incr<IncrType>::type> >,
                        typename internal::cleanup_seq_incr<IncrType>::type> >::type
seq(FirstType f, const symbolic::BaseExpr<LastTypeDerived> &l, IncrType incr)
{
  typedef typename internal::cleanup_seq_incr<IncrType>::type CleanedIncrType;
  return seqN(typename internal::cleanup_index_type<FirstType>::type(f),
              (l.derived()-typename internal::cleanup_index_type<FirstType>::type(f)+CleanedIncrType(incr))/CleanedIncrType(incr), incr);
}

template<typename FirstTypeDerived,typename LastTypeDerived, typename IncrType>
ArithmeticSequence<FirstTypeDerived,
                    symbolic::QuotientExpr<symbolic::AddExpr<symbolic::AddExpr<LastTypeDerived,
                                                                               symbolic::NegateExpr<FirstTypeDerived> >,
                                                             symbolic::ValueExpr<typename internal::cleanup_seq_incr<IncrType>::type> >,
                                          symbolic::ValueExpr<typename internal::cleanup_seq_incr<IncrType>::type> >,
                    typename internal::cleanup_seq_incr<IncrType>::type>
seq(const symbolic::BaseExpr<FirstTypeDerived> &f, const symbolic::BaseExpr<LastTypeDerived> &l, IncrType incr)
{
  typedef typename internal::cleanup_seq_incr<IncrType>::type CleanedIncrType;
  return seqN(f.derived(),(l.derived()-f.derived()+CleanedIncrType(incr))/CleanedIncrType(incr), incr);
}
#endif 

#endif 

#if EIGEN_HAS_CXX11 || defined(EIGEN_PARSED_BY_DOXYGEN)

template<typename SizeType,typename IncrType>
auto lastN(SizeType size, IncrType incr)
-> decltype(seqN(Eigen::last-(size-fix<1>())*incr, size, incr))
{
  return seqN(Eigen::last-(size-fix<1>())*incr, size, incr);
}

template<typename SizeType>
auto lastN(SizeType size)
-> decltype(seqN(Eigen::last+fix<1>()-size, size))
{
  return seqN(Eigen::last+fix<1>()-size, size);
}
#endif

namespace internal {

template<typename T>
struct make_size_type {
  typedef typename internal::conditional<symbolic::is_symbolic<T>::value, Index, T>::type type;
};

template<typename FirstType,typename SizeType,typename IncrType,int XprSize>
struct IndexedViewCompatibleType<ArithmeticSequence<FirstType,SizeType,IncrType>, XprSize> {
  typedef ArithmeticSequence<Index,typename make_size_type<SizeType>::type,IncrType> type;
};

template<typename FirstType,typename SizeType,typename IncrType>
ArithmeticSequence<Index,typename make_size_type<SizeType>::type,IncrType>
makeIndexedViewCompatible(const ArithmeticSequence<FirstType,SizeType,IncrType>& ids, Index size,SpecializedType) {
  return ArithmeticSequence<Index,typename make_size_type<SizeType>::type,IncrType>(
            eval_expr_given_size(ids.firstObject(),size),eval_expr_given_size(ids.sizeObject(),size),ids.incrObject());
}

template<typename FirstType,typename SizeType,typename IncrType>
struct get_compile_time_incr<ArithmeticSequence<FirstType,SizeType,IncrType> > {
  enum { value = get_fixed_value<IncrType,DynamicIndex>::value };
};

} 

namespace indexing {
  using Eigen::all;
  using Eigen::seq;
  using Eigen::seqN;
  #if EIGEN_HAS_CXX11
  using Eigen::lastN;
  #endif
  using Eigen::last;
  using Eigen::lastp1;
  using Eigen::fix;
}

} 

#endif 
