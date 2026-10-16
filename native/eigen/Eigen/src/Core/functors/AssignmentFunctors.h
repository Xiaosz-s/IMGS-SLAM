#ifndef EIGEN_ASSIGNMENT_FUNCTORS_H
#define EIGEN_ASSIGNMENT_FUNCTORS_H

namespace Eigen {

namespace internal {

template<typename DstScalar,typename SrcScalar> struct assign_op {

  EIGEN_EMPTY_STRUCT_CTOR(assign_op)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void assignCoeff(DstScalar& a, const SrcScalar& b) const { a = b; }
  
  template<int Alignment, typename Packet>
  EIGEN_STRONG_INLINE void assignPacket(DstScalar* a, const Packet& b) const
  { internal::pstoret<DstScalar,Packet,Alignment>(a,b); }
};

template<typename DstScalar> struct assign_op<DstScalar,void> {};

template<typename DstScalar,typename SrcScalar>
struct functor_traits<assign_op<DstScalar,SrcScalar> > {
  enum {
    Cost = NumTraits<DstScalar>::ReadCost,
    PacketAccess = is_same<DstScalar,SrcScalar>::value && packet_traits<DstScalar>::Vectorizable && packet_traits<SrcScalar>::Vectorizable
  };
};

template<typename DstScalar,typename SrcScalar> struct add_assign_op {

  EIGEN_EMPTY_STRUCT_CTOR(add_assign_op)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void assignCoeff(DstScalar& a, const SrcScalar& b) const { a += b; }
  
  template<int Alignment, typename Packet>
  EIGEN_STRONG_INLINE void assignPacket(DstScalar* a, const Packet& b) const
  { internal::pstoret<DstScalar,Packet,Alignment>(a,internal::padd(internal::ploadt<Packet,Alignment>(a),b)); }
};
template<typename DstScalar,typename SrcScalar>
struct functor_traits<add_assign_op<DstScalar,SrcScalar> > {
  enum {
    Cost = NumTraits<DstScalar>::ReadCost + NumTraits<DstScalar>::AddCost,
    PacketAccess = is_same<DstScalar,SrcScalar>::value && packet_traits<DstScalar>::HasAdd
  };
};

template<typename DstScalar,typename SrcScalar> struct sub_assign_op {

  EIGEN_EMPTY_STRUCT_CTOR(sub_assign_op)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void assignCoeff(DstScalar& a, const SrcScalar& b) const { a -= b; }
  
  template<int Alignment, typename Packet>
  EIGEN_STRONG_INLINE void assignPacket(DstScalar* a, const Packet& b) const
  { internal::pstoret<DstScalar,Packet,Alignment>(a,internal::psub(internal::ploadt<Packet,Alignment>(a),b)); }
};
template<typename DstScalar,typename SrcScalar>
struct functor_traits<sub_assign_op<DstScalar,SrcScalar> > {
  enum {
    Cost = NumTraits<DstScalar>::ReadCost + NumTraits<DstScalar>::AddCost,
    PacketAccess = is_same<DstScalar,SrcScalar>::value && packet_traits<DstScalar>::HasSub
  };
};

template<typename DstScalar, typename SrcScalar=DstScalar>
struct mul_assign_op {

  EIGEN_EMPTY_STRUCT_CTOR(mul_assign_op)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void assignCoeff(DstScalar& a, const SrcScalar& b) const { a *= b; }
  
  template<int Alignment, typename Packet>
  EIGEN_STRONG_INLINE void assignPacket(DstScalar* a, const Packet& b) const
  { internal::pstoret<DstScalar,Packet,Alignment>(a,internal::pmul(internal::ploadt<Packet,Alignment>(a),b)); }
};
template<typename DstScalar, typename SrcScalar>
struct functor_traits<mul_assign_op<DstScalar,SrcScalar> > {
  enum {
    Cost = NumTraits<DstScalar>::ReadCost + NumTraits<DstScalar>::MulCost,
    PacketAccess = is_same<DstScalar,SrcScalar>::value && packet_traits<DstScalar>::HasMul
  };
};

template<typename DstScalar, typename SrcScalar=DstScalar> struct div_assign_op {

  EIGEN_EMPTY_STRUCT_CTOR(div_assign_op)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void assignCoeff(DstScalar& a, const SrcScalar& b) const { a /= b; }
  
  template<int Alignment, typename Packet>
  EIGEN_STRONG_INLINE void assignPacket(DstScalar* a, const Packet& b) const
  { internal::pstoret<DstScalar,Packet,Alignment>(a,internal::pdiv(internal::ploadt<Packet,Alignment>(a),b)); }
};
template<typename DstScalar, typename SrcScalar>
struct functor_traits<div_assign_op<DstScalar,SrcScalar> > {
  enum {
    Cost = NumTraits<DstScalar>::ReadCost + NumTraits<DstScalar>::MulCost,
    PacketAccess = is_same<DstScalar,SrcScalar>::value && packet_traits<DstScalar>::HasDiv
  };
};

template<typename Scalar> struct swap_assign_op {

  EIGEN_EMPTY_STRUCT_CTOR(swap_assign_op)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void assignCoeff(Scalar& a, const Scalar& b) const
  {
#ifdef EIGEN_GPUCC
    
    Scalar t=b; const_cast<Scalar&>(b)=a; a=t;
#else
    using std::swap;
    swap(a,const_cast<Scalar&>(b));
#endif
  }
};
template<typename Scalar>
struct functor_traits<swap_assign_op<Scalar> > {
  enum {
    Cost = 3 * NumTraits<Scalar>::ReadCost,
    PacketAccess = 
    #if defined(EIGEN_VECTORIZE_AVX) && EIGEN_COMP_CLANG && (EIGEN_COMP_CLANG<800 || defined(__apple_build_version__))

    0
    #else
    packet_traits<Scalar>::Vectorizable
    #endif
  };
};

} 

} 

#endif 
