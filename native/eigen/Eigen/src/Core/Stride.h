#ifndef EIGEN_STRIDE_H
#define EIGEN_STRIDE_H

namespace Eigen {

template<int OuterStrideAtCompileTime_, int InnerStrideAtCompileTime_>
class Stride
{
  public:
    typedef Eigen::Index Index; 
    enum {
      InnerStrideAtCompileTime = InnerStrideAtCompileTime_,
      OuterStrideAtCompileTime = OuterStrideAtCompileTime_
    };

    EIGEN_DEVICE_FUNC
    Stride()
      : m_outer(OuterStrideAtCompileTime), m_inner(InnerStrideAtCompileTime)
    {

      eigen_assert(InnerStrideAtCompileTime != Dynamic && OuterStrideAtCompileTime != Dynamic);
    }

    EIGEN_DEVICE_FUNC
    Stride(Index outerStride, Index innerStride)
      : m_outer(outerStride), m_inner(innerStride)
    {
    }

    EIGEN_DEVICE_FUNC
    Stride(const Stride& other)
      : m_outer(other.outer()), m_inner(other.inner())
    {}

    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index outer() const { return m_outer.value(); }
    
    EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
    inline Index inner() const { return m_inner.value(); }

  protected:
    internal::variable_if_dynamic<Index, OuterStrideAtCompileTime> m_outer;
    internal::variable_if_dynamic<Index, InnerStrideAtCompileTime> m_inner;
};

template<int Value>
class InnerStride : public Stride<0, Value>
{
    typedef Stride<0, Value> Base;
  public:
    EIGEN_DEVICE_FUNC InnerStride() : Base() {}
    EIGEN_DEVICE_FUNC InnerStride(Index v) : Base(0, v) {} 
};

template<int Value>
class OuterStride : public Stride<Value, 0>
{
    typedef Stride<Value, 0> Base;
  public:
    EIGEN_DEVICE_FUNC OuterStride() : Base() {}
    EIGEN_DEVICE_FUNC OuterStride(Index v) : Base(v,0) {} 
};

} 

#endif 
