#ifndef EIGEN_COREITERATORS_H
#define EIGEN_COREITERATORS_H

namespace Eigen { 

namespace internal {

template<typename XprType, typename EvaluatorKind>
class inner_iterator_selector;

}

template<typename XprType>
class InnerIterator
{
protected:
  typedef internal::inner_iterator_selector<XprType, typename internal::evaluator_traits<XprType>::Kind> IteratorType;
  typedef internal::evaluator<XprType> EvaluatorType;
  typedef typename internal::traits<XprType>::Scalar Scalar;
public:
  
  InnerIterator(const XprType &xpr, const Index &outerId)
    : m_eval(xpr), m_iter(m_eval, outerId, xpr.innerSize())
  {}

  EIGEN_STRONG_INLINE Scalar value() const          { return m_iter.value(); }

  EIGEN_STRONG_INLINE InnerIterator& operator++()   { m_iter.operator++(); return *this; }
  EIGEN_STRONG_INLINE InnerIterator& operator+=(Index i) { m_iter.operator+=(i); return *this; }
  EIGEN_STRONG_INLINE InnerIterator operator+(Index i) 
  { InnerIterator result(*this); result+=i; return result; }

  EIGEN_STRONG_INLINE Index index() const           { return m_iter.index(); }
  
  EIGEN_STRONG_INLINE Index row() const             { return m_iter.row(); }
  
  EIGEN_STRONG_INLINE Index col() const             { return m_iter.col(); }
  
  EIGEN_STRONG_INLINE operator bool() const         { return m_iter; }
  
protected:
  EvaluatorType m_eval;
  IteratorType m_iter;
private:

  template<typename T> InnerIterator(const EigenBase<T>&,Index outer);
};

namespace internal {

template<typename XprType>
class inner_iterator_selector<XprType, IndexBased>
{
protected:
  typedef evaluator<XprType> EvaluatorType;
  typedef typename traits<XprType>::Scalar Scalar;
  enum { IsRowMajor = (XprType::Flags&RowMajorBit)==RowMajorBit };
  
public:
  EIGEN_STRONG_INLINE inner_iterator_selector(const EvaluatorType &eval, const Index &outerId, const Index &innerSize)
    : m_eval(eval), m_inner(0), m_outer(outerId), m_end(innerSize)
  {}

  EIGEN_STRONG_INLINE Scalar value() const
  {
    return (IsRowMajor) ? m_eval.coeff(m_outer, m_inner)
                        : m_eval.coeff(m_inner, m_outer);
  }

  EIGEN_STRONG_INLINE inner_iterator_selector& operator++() { m_inner++; return *this; }

  EIGEN_STRONG_INLINE Index index() const { return m_inner; }
  inline Index row() const { return IsRowMajor ? m_outer : index(); }
  inline Index col() const { return IsRowMajor ? index() : m_outer; }

  EIGEN_STRONG_INLINE operator bool() const { return m_inner < m_end && m_inner>=0; }

protected:
  const EvaluatorType& m_eval;
  Index m_inner;
  const Index m_outer;
  const Index m_end;
};

template<typename XprType>
class inner_iterator_selector<XprType, IteratorBased>
 : public evaluator<XprType>::InnerIterator
{
protected:
  typedef typename evaluator<XprType>::InnerIterator Base;
  typedef evaluator<XprType> EvaluatorType;
  
public:
  EIGEN_STRONG_INLINE inner_iterator_selector(const EvaluatorType &eval, const Index &outerId, const Index &)
    : Base(eval, outerId)
  {}  
};

} 

} 

#endif 
