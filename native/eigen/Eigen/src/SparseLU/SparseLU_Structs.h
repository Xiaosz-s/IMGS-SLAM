#ifndef EIGEN_LU_STRUCTS
#define EIGEN_LU_STRUCTS
namespace Eigen {
namespace internal {
  
typedef enum {LUSUP, UCOL, LSUB, USUB, LLVL, ULVL} MemType; 

template <typename IndexVector, typename ScalarVector>
struct LU_GlobalLU_t {
  typedef typename IndexVector::Scalar StorageIndex; 
  IndexVector xsup; 
  IndexVector supno; 
  ScalarVector  lusup; 
  IndexVector lsub; 
  IndexVector xlusup; 
  IndexVector xlsub; 
  Index   nzlmax; 
  Index   nzlumax; 
  ScalarVector  ucol; 
  IndexVector usub; 
  IndexVector xusub; 
  Index   nzumax; 
  Index   n; 
  Index   num_expansions; 
};

struct perfvalues {
  Index panel_size; 
  Index relax; 

  Index maxsuper; 
  Index rowblk; 
  Index colblk; 
  Index fillfactor; 
}; 

} 

} 
#endif 
