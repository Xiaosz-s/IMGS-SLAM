#ifndef EIGEN_CONSTANTS_H
#define EIGEN_CONSTANTS_H

namespace Eigen {

const int Dynamic = -1;

const int DynamicIndex = 0xffffff;

const int UndefinedIncr = 0xfffffe;

const int Infinity = -1;

const int HugeCost = 10000;

const unsigned int RowMajorBit = 0x1;

const unsigned int EvalBeforeNestingBit = 0x2;

EIGEN_DEPRECATED
const unsigned int EvalBeforeAssigningBit = 0x4; 

const unsigned int PacketAccessBit = 0x8;

#ifdef EIGEN_VECTORIZE

const unsigned int ActualPacketAccessBit = PacketAccessBit;
#else
const unsigned int ActualPacketAccessBit = 0x0;
#endif

const unsigned int LinearAccessBit = 0x10;

const unsigned int LvalueBit = 0x20;

const unsigned int DirectAccessBit = 0x40;

EIGEN_DEPRECATED const unsigned int AlignedBit = 0x80;

const unsigned int NestByRefBit = 0x100;

const unsigned int NoPreferredStorageOrderBit = 0x200;

const unsigned int CompressedAccessBit = 0x400;

const unsigned int HereditaryBits = RowMajorBit
                                  | EvalBeforeNestingBit;

enum UpLoType {
  
  Lower=0x1,                      
  
  Upper=0x2,                      
  
  UnitDiag=0x4, 
  
  ZeroDiag=0x8,
  
  UnitLower=UnitDiag|Lower, 
  
  UnitUpper=UnitDiag|Upper,
  
  StrictlyLower=ZeroDiag|Lower, 
  
  StrictlyUpper=ZeroDiag|Upper,
  
  SelfAdjoint=0x10,
  
  Symmetric=0x20
};

enum AlignmentType {
  Unaligned=0,        
  Aligned8=8,         
  Aligned16=16,       
  Aligned32=32,       
  Aligned64=64,       
  Aligned128=128,     
  AlignedMask=255,
  Aligned=16,         
#if EIGEN_MAX_ALIGN_BYTES==128
  AlignedMax = Aligned128
#elif EIGEN_MAX_ALIGN_BYTES==64
  AlignedMax = Aligned64
#elif EIGEN_MAX_ALIGN_BYTES==32
  AlignedMax = Aligned32
#elif EIGEN_MAX_ALIGN_BYTES==16
  AlignedMax = Aligned16
#elif EIGEN_MAX_ALIGN_BYTES==8
  AlignedMax = Aligned8
#elif EIGEN_MAX_ALIGN_BYTES==0
  AlignedMax = Unaligned
#else
#error Invalid value for EIGEN_MAX_ALIGN_BYTES
#endif
};

enum DirectionType { 

  Vertical, 

  Horizontal, 

  BothDirections 
};

enum TraversalType {
  
  DefaultTraversal,
  
  LinearTraversal,

  InnerVectorizedTraversal,

  LinearVectorizedTraversal,

  SliceVectorizedTraversal,
  
  InvalidTraversal,
  
  AllAtOnceTraversal
};

enum UnrollingType {
  
  NoUnrolling,
  
  InnerUnrolling,

  CompleteUnrolling
};

enum SpecializedType {
  Specialized,
  BuiltIn
};

enum StorageOptions {
  
  ColMajor = 0,
  
  RowMajor = 0x1,  
  
  AutoAlign = 0,
   
  DontAlign = 0x2
};

enum SideType {
  
  OnTheLeft = 1,
  
  OnTheRight = 2
};

enum NaNPropagationOptions {
  
  PropagateFast = 0,
  
  PropagateNaN,
  
  PropagateNumbers
};

enum NoChange_t   { NoChange };
enum Sequential_t { Sequential };
enum Default_t    { Default };

enum AmbiVectorMode {
  IsDense         = 0,
  IsSparse
};

enum AccessorLevels {
  
  ReadOnlyAccessors, 
  
  WriteAccessors, 
  
  DirectAccessors, 
  
  DirectWriteAccessors
};

enum DecompositionOptions {
  
  Pivoting            = 0x01, 
  
  NoPivoting          = 0x02, 
  
  ComputeFullU        = 0x04,
  
  ComputeThinU        = 0x08,
  
  ComputeFullV        = 0x10,
  
  ComputeThinV        = 0x20,

  EigenvaluesOnly     = 0x40,

  ComputeEigenvectors = 0x80,
  
  EigVecMask = EigenvaluesOnly | ComputeEigenvectors,

  Ax_lBx              = 0x100,

  ABx_lx              = 0x200,

  BAx_lx              = 0x400,
  
  GenEigMask = Ax_lBx | ABx_lx | BAx_lx
};

enum QRPreconditioners {
  
  NoQRPreconditioner,
  
  HouseholderQRPreconditioner,
  
  ColPivHouseholderQRPreconditioner,
  
  FullPivHouseholderQRPreconditioner
};

#ifdef Success
#error The preprocessor symbol 'Success' is defined, possibly by the X11 header file X.h
#endif

enum ComputationInfo {
  
  Success = 0,        
  
  NumericalIssue = 1, 
  
  NoConvergence = 2,

  InvalidInput = 3
};

enum TransformTraits {
  
  Isometry      = 0x1,

  Affine        = 0x2,
  
  AffineCompact = 0x10 | Affine,
  
  Projective    = 0x20
};

namespace Architecture
{
  enum Type {
    Generic = 0x0,
    SSE = 0x1,
    AltiVec = 0x2,
    VSX = 0x3,
    NEON = 0x4,
    MSA = 0x5,
    SVE = 0x6,
#if defined EIGEN_VECTORIZE_SSE
    Target = SSE
#elif defined EIGEN_VECTORIZE_ALTIVEC
    Target = AltiVec
#elif defined EIGEN_VECTORIZE_VSX
    Target = VSX
#elif defined EIGEN_VECTORIZE_NEON
    Target = NEON
#elif defined EIGEN_VECTORIZE_SVE
    Target = SVE
#elif defined EIGEN_VECTORIZE_MSA
    Target = MSA
#else
    Target = Generic
#endif
  };
}

enum ProductImplType
{ DefaultProduct=0, LazyProduct, AliasFreeProduct, CoeffBasedProductMode, LazyCoeffBasedProductMode, OuterProduct, InnerProduct, GemvProduct, GemmProduct };

enum Action {GetAction, SetAction};

struct Dense {};

struct Sparse {};

struct SolverStorage {};

struct PermutationStorage {};

struct TranspositionsStorage {};

struct MatrixXpr {};

struct ArrayXpr {};

struct DenseShape             { static std::string debugName() { return "DenseShape"; } };
struct SolverShape            { static std::string debugName() { return "SolverShape"; } };
struct HomogeneousShape       { static std::string debugName() { return "HomogeneousShape"; } };
struct DiagonalShape          { static std::string debugName() { return "DiagonalShape"; } };
struct BandShape              { static std::string debugName() { return "BandShape"; } };
struct TriangularShape        { static std::string debugName() { return "TriangularShape"; } };
struct SelfAdjointShape       { static std::string debugName() { return "SelfAdjointShape"; } };
struct PermutationShape       { static std::string debugName() { return "PermutationShape"; } };
struct TranspositionsShape    { static std::string debugName() { return "TranspositionsShape"; } };
struct SparseShape            { static std::string debugName() { return "SparseShape"; } };

namespace internal {

struct IndexBased {};

struct IteratorBased {};

enum ComparisonName {
  cmp_EQ = 0,
  cmp_LT = 1,
  cmp_LE = 2,
  cmp_UNORD = 3,
  cmp_NEQ = 4,
  cmp_GT = 5,
  cmp_GE = 6
};
} 

} 

#endif 
