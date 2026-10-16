#include <torch/extension.h>
#include <cuda_runtime.h>
#include <c10/cuda/CUDAGuard.h>
#include <c10/cuda/CUDAStream.h>
#include <c10/cuda/CUDAException.h>
#include <cstdint>
template<class T> T* take(char*& p,size_t count){
 uintptr_t q=(reinterpret_cast<uintptr_t>(p)+127)&~uintptr_t(127);
 T* out=reinterpret_cast<T*>(q);p=reinterpret_cast<char*>(out+count);return out;
}
__global__ void score_kernel(const float2* xy,const float4* conic,const uint32_t* ids,
 const uint2* ranges,const int* metric,int* counts,int width,int height){
 int x=blockIdx.x*16+threadIdx.x,y=blockIdx.y*16+threadIdx.y;
 if(x>=width||y>=height||metric[y*width+x]==0)return;
 uint2 range=ranges[blockIdx.y*((width+15)/16)+blockIdx.x];float T=1.f;
 for(uint32_t k=range.x;k<range.y;++k){
  int id=ids[k];float dx=xy[id].x-x,dy=xy[id].y-y;float4 c=conic[id];
  float power=-.5f*(c.x*dx*dx+c.z*dy*dy)-c.y*dx*dy;
  if(power>0)continue;
  float alpha=fminf(.99f,c.w*expf(power));if(alpha<1.f/255.f)continue;
  float next=T*(1-alpha);if(next<.0001f)break;
  atomicAdd(counts+id,1);T=next;
 }
}
torch::Tensor ScoreGaussiansCUDA(torch::Tensor geom,torch::Tensor binning,torch::Tensor image,
 torch::Tensor metric,int points,int instances){
 TORCH_CHECK(metric.is_cuda()&&geom.is_cuda()&&binning.is_cuda()&&image.is_cuda(),"CUDA buffers required");
 TORCH_CHECK(metric.dim()==2&&metric.scalar_type()==torch::kInt32&&metric.is_contiguous(),"metric must be contiguous int32 HxW");
 TORCH_CHECK(geom.scalar_type()==torch::kUInt8&&binning.scalar_type()==torch::kUInt8&&image.scalar_type()==torch::kUInt8,"byte buffers required");
 c10::cuda::CUDAGuard guard(metric.device());
 auto counts=torch::zeros({points},metric.options());if(points==0||instances==0)return counts;
 char* p=reinterpret_cast<char*>(geom.data_ptr());
 take<float>(p,points);take<float2>(p,points);take<float>(p,points);take<bool>(p,points*3);take<int>(p,points);
 auto xy=take<float2>(p,points);take<float>(p,points*3);take<float>(p,points*6);auto conic=take<float4>(p,points);
 p=reinterpret_cast<char*>(binning.data_ptr());auto ids=take<uint32_t>(p,instances);
 int h=metric.size(0),w=metric.size(1);p=reinterpret_cast<char*>(image.data_ptr());take<uint32_t>(p,w*h);auto ranges=take<uint2>(p,w*h);
 score_kernel<<<dim3((w+15)/16,(h+15)/16),dim3(16,16),0,c10::cuda::getCurrentCUDAStream()>>>(xy,conic,ids,ranges,metric.data_ptr<int>(),counts.data_ptr<int>(),w,h);
 C10_CUDA_KERNEL_LAUNCH_CHECK();return counts;
}
