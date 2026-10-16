#ifndef EIGEN_SWAP_H
#define EIGEN_SWAP_H

namespace Eigen { 

namespace internal {

template<typename DstEvaluatorTypeT, typename SrcEvaluatorTypeT>
class generic_dense_assignment_kernel<DstEvaluatorTypeT, SrcEvaluatorTypeT, swap_assign_op<typename DstEvaluatorTypeT::Scalar>, Specialized>
 : public generic_dense_assignment_kernel<DstEvaluatorTypeT, SrcEvaluatorTypeT, swap_assign_op<typename DstEvaluatorTypeT::Scalar>, BuiltIn>
{
protected:
  typedef generic_dense_assignment_kernel<DstEvaluatorTypeT, SrcEvaluatorTypeT, swap_assign_op<typename DstEvaluatorTypeT::Scalar>, BuiltIn> Base;
  using Base::m_dst;
  using Base::m_src;
  using Base::m_functor;
  
public:
  typedef typename Base::Scalar Scalar;
  typedef typename Base::DstXprType DstXprType;
  typedef swap_assign_op<Scalar> Functor;
  
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  generic_dense_assignment_kernel(DstEvaluatorTypeT &dst, const SrcEvaluatorTypeT &src, const Functor &func, DstXprType& dstExpr)
    : Base(dst, src, func, dstExpr)
  {}
  
  template<int StoreMode, int LoadMode, typename PacketType>
  EIGEN_STRONG_INLINE void assignPacket(Index row, Index col)
  {
    PacketType tmp = m_src.template packet<LoadMode,PacketType>(row,col);
    const_cast<SrcEvaluatorTypeT&>(m_src).template writePacket<LoadMode>(row,col, m_dst.template packet<StoreMode,PacketType>(row,col));
    m_dst.template writePacket<StoreMode>(row,col,tmp);
  }
  
  template<int StoreMode, int LoadMode, typename PacketType>
  EIGEN_STRONG_INLINE void assignPacket(Index index)
  {
    PacketType tmp = m_src.template packet<LoadMode,PacketType>(index);
    const_cast<SrcEvaluatorTypeT&>(m_src).template writePacket<LoadMode>(index, m_dst.template packet<StoreMode,PacketType>(index));
    m_dst.template writePacket<StoreMode>(index,tmp);
  }

  template<int StoreMode, int LoadMode, typename PacketType>
  EIGEN_STRONG_INLINE void assignPacketByOuterInner(Index outer, Index inner)
  {
    Index row = Base::rowIndexByOuterInner(outer, inner); 
    Index col = Base::colIndexByOuterInner(outer, inner);
    assignPacket<StoreMode,LoadMode,PacketType>(row, col);
  }
};

} 

} 

#endif 
