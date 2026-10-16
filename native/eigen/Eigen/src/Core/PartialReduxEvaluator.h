#ifndef EIGEN_PARTIALREDUX_H
#define EIGEN_PARTIALREDUX_H

namespace Eigen { 

namespace internal {

template<typename Func, typename Evaluator>
struct packetwise_redux_traits
{
  enum {
    OuterSize = int(Evaluator::IsRowMajor) ? Evaluator::RowsAtCompileTime : Evaluator::ColsAtCompileTime,
    Cost = OuterSize == Dynamic ? HugeCost
         : OuterSize * Evaluator::CoeffReadCost + (OuterSize-1) * functor_traits<Func>::Cost,
    Unrolling = Cost <= EIGEN_UNROLLING_LIMIT ? CompleteUnrolling : NoUnrolling
  };

};

template<typename PacketType,typename Func>
EIGEN_DEVICE_FUNC
PacketType packetwise_redux_empty_value(const Func& ) { return pset1<PacketType>(0); }

template<typename PacketType,typename Scalar>
EIGEN_DEVICE_FUNC
PacketType packetwise_redux_empty_value(const scalar_product_op<Scalar,Scalar>& ) { return pset1<PacketType>(1); }

template<typename Func, typename Evaluator,
         int Unrolling = packetwise_redux_traits<Func, Evaluator>::Unrolling
>
struct packetwise_redux_impl;

template<typename Func, typename Evaluator>
struct packetwise_redux_impl<Func, Evaluator, CompleteUnrolling>
{
  typedef redux_novec_unroller<Func,Evaluator, 0, Evaluator::SizeAtCompileTime> Base;
  typedef typename Evaluator::Scalar Scalar;

  template<typename PacketType>
  EIGEN_DEVICE_FUNC static EIGEN_STRONG_INLINE
  PacketType run(const Evaluator &eval, const Func& func, Index )
  {
    return redux_vec_unroller<Func, Evaluator, 0, packetwise_redux_traits<Func, Evaluator>::OuterSize>::template run<PacketType>(eval,func);
  }
};

template<typename Func, typename Evaluator, int Start>
struct redux_vec_unroller<Func, Evaluator, Start, 0>
{
  template<typename PacketType>
  EIGEN_DEVICE_FUNC
  static EIGEN_STRONG_INLINE PacketType run(const Evaluator &, const Func& f)
  {
    return packetwise_redux_empty_value<PacketType>(f);
  }
};

template<typename Func, typename Evaluator>
struct packetwise_redux_impl<Func, Evaluator, NoUnrolling>
{
  typedef typename Evaluator::Scalar Scalar;
  typedef typename redux_traits<Func, Evaluator>::PacketType PacketScalar;

  template<typename PacketType>
  EIGEN_DEVICE_FUNC
  static PacketType run(const Evaluator &eval, const Func& func, Index size)
  {
    if(size==0)
      return packetwise_redux_empty_value<PacketType>(func);
    
    const Index size4 = (size-1)&(~3);
    PacketType p = eval.template packetByOuterInner<Unaligned,PacketType>(0,0);
    Index i = 1;

    for(; i<size4; i+=4)
      p = func.packetOp(p,
            func.packetOp(
              func.packetOp(eval.template packetByOuterInner<Unaligned,PacketType>(i+0,0),eval.template packetByOuterInner<Unaligned,PacketType>(i+1,0)),
              func.packetOp(eval.template packetByOuterInner<Unaligned,PacketType>(i+2,0),eval.template packetByOuterInner<Unaligned,PacketType>(i+3,0))));
    for(; i<size; ++i)
      p = func.packetOp(p, eval.template packetByOuterInner<Unaligned,PacketType>(i,0));
    return p;
  }
};

template< typename ArgType, typename MemberOp, int Direction>
struct evaluator<PartialReduxExpr<ArgType, MemberOp, Direction> >
  : evaluator_base<PartialReduxExpr<ArgType, MemberOp, Direction> >
{
  typedef PartialReduxExpr<ArgType, MemberOp, Direction> XprType;
  typedef typename internal::nested_eval<ArgType,1>::type ArgTypeNested;
  typedef typename internal::add_const_on_value_type<ArgTypeNested>::type ConstArgTypeNested;
  typedef typename internal::remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
  typedef typename ArgType::Scalar InputScalar;
  typedef typename XprType::Scalar Scalar;
  enum {
    TraversalSize = Direction==int(Vertical) ? int(ArgType::RowsAtCompileTime) :  int(ArgType::ColsAtCompileTime)
  };
  typedef typename MemberOp::template Cost<int(TraversalSize)> CostOpType;
  enum {
    CoeffReadCost = TraversalSize==Dynamic ? HugeCost
                  : TraversalSize==0 ? 1
                  : int(TraversalSize) * int(evaluator<ArgType>::CoeffReadCost) + int(CostOpType::value),
    
    _ArgFlags = evaluator<ArgType>::Flags,

    _Vectorizable =  bool(int(_ArgFlags)&PacketAccessBit)
                  && bool(MemberOp::Vectorizable)
                  && (Direction==int(Vertical) ? bool(_ArgFlags&RowMajorBit) : (_ArgFlags&RowMajorBit)==0)
                  && (TraversalSize!=0),
                  
    Flags = (traits<XprType>::Flags&RowMajorBit)
          | (evaluator<ArgType>::Flags&(HereditaryBits&(~RowMajorBit)))
          | (_Vectorizable ? PacketAccessBit : 0)
          | LinearAccessBit,
    
    Alignment = 0 
  };

  EIGEN_DEVICE_FUNC explicit evaluator(const XprType xpr)
    : m_arg(xpr.nestedExpression()), m_functor(xpr.functor())
  {
    EIGEN_INTERNAL_CHECK_COST_VALUE(TraversalSize==Dynamic ? HugeCost : (TraversalSize==0 ? 1 : int(CostOpType::value)));
    EIGEN_INTERNAL_CHECK_COST_VALUE(CoeffReadCost);
  }

  typedef typename XprType::CoeffReturnType CoeffReturnType;

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  const Scalar coeff(Index i, Index j) const
  {
    return coeff(Direction==Vertical ? j : i);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  const Scalar coeff(Index index) const
  {
    return m_functor(m_arg.template subVector<DirectionType(Direction)>(index));
  }

  template<int LoadMode,typename PacketType>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  PacketType packet(Index i, Index j) const
  {
    return packet<LoadMode,PacketType>(Direction==Vertical ? j : i);
  }
  
  template<int LoadMode,typename PacketType>
  EIGEN_STRONG_INLINE EIGEN_DEVICE_FUNC
  PacketType packet(Index idx) const
  {
    enum { PacketSize = internal::unpacket_traits<PacketType>::size };
    typedef Block<const ArgTypeNestedCleaned,
                  Direction==Vertical ? int(ArgType::RowsAtCompileTime) : int(PacketSize),
                  Direction==Vertical ? int(PacketSize) : int(ArgType::ColsAtCompileTime),
                  true > PanelType;
    
    PanelType panel(m_arg,
                    Direction==Vertical ? 0 : idx,
                    Direction==Vertical ? idx : 0,
                    Direction==Vertical ? m_arg.rows() : Index(PacketSize),
                    Direction==Vertical ? Index(PacketSize) : m_arg.cols());

    if(PacketSize==1)
      return internal::pset1<PacketType>(coeff(idx));
    
    typedef typename internal::redux_evaluator<PanelType> PanelEvaluator;
    PanelEvaluator panel_eval(panel);
    typedef typename MemberOp::BinaryOp BinaryOp;
    PacketType p = internal::packetwise_redux_impl<BinaryOp,PanelEvaluator>::template run<PacketType>(panel_eval,m_functor.binaryFunc(),m_arg.outerSize());
    return p;
  }

protected:
  ConstArgTypeNested m_arg;
  const MemberOp m_functor;
};

} 

} 

#endif 
