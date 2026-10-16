#ifndef EIGEN_CONFIGURE_VECTORIZATION_H
#define EIGEN_CONFIGURE_VECTORIZATION_H

#if (defined EIGEN_CUDACC)
  #define EIGEN_ALIGN_TO_BOUNDARY(n) __align__(n)
  #define EIGEN_ALIGNOF(x) __alignof(x)
#elif EIGEN_HAS_ALIGNAS
  #define EIGEN_ALIGN_TO_BOUNDARY(n) alignas(n)
  #define EIGEN_ALIGNOF(x) alignof(x)
#elif EIGEN_COMP_GNUC || EIGEN_COMP_PGI || EIGEN_COMP_IBM || EIGEN_COMP_ARM
  #define EIGEN_ALIGN_TO_BOUNDARY(n) __attribute__((aligned(n)))
  #define EIGEN_ALIGNOF(x) __alignof(x)
#elif EIGEN_COMP_MSVC
  #define EIGEN_ALIGN_TO_BOUNDARY(n) __declspec(align(n))
  #define EIGEN_ALIGNOF(x) __alignof(x)
#elif EIGEN_COMP_SUNCC
  
  #define EIGEN_ALIGN_TO_BOUNDARY(n) __attribute__((aligned(n)))
  #define EIGEN_ALIGNOF(x) __alignof(x)
#else
  #error Please tell me what is the equivalent of alignas(n) and alignof(x) for your compiler
#endif

#if defined(EIGEN_DONT_VECTORIZE)
  #if defined(EIGEN_GPUCC)

    #define EIGEN_IDEAL_MAX_ALIGN_BYTES 16
  #else
    #define EIGEN_IDEAL_MAX_ALIGN_BYTES 0
  #endif
#elif defined(__AVX512F__)
  
  #define EIGEN_IDEAL_MAX_ALIGN_BYTES 64
#elif defined(__AVX__)
  
  #define EIGEN_IDEAL_MAX_ALIGN_BYTES 32
#else
  #define EIGEN_IDEAL_MAX_ALIGN_BYTES 16
#endif

#define EIGEN_MIN_ALIGN_BYTES 16

#if (defined(EIGEN_DONT_ALIGN_STATICALLY) || defined(EIGEN_DONT_ALIGN))  && defined(EIGEN_MAX_STATIC_ALIGN_BYTES) && EIGEN_MAX_STATIC_ALIGN_BYTES>0
#error EIGEN_MAX_STATIC_ALIGN_BYTES and EIGEN_DONT_ALIGN[_STATICALLY] are both defined with EIGEN_MAX_STATIC_ALIGN_BYTES!=0. Use EIGEN_MAX_STATIC_ALIGN_BYTES=0 as a synonym of EIGEN_DONT_ALIGN_STATICALLY.
#endif

#if defined(EIGEN_DONT_ALIGN_STATICALLY) || defined(EIGEN_DONT_ALIGN)
  #ifdef EIGEN_MAX_STATIC_ALIGN_BYTES
    #undef EIGEN_MAX_STATIC_ALIGN_BYTES
  #endif
  #define EIGEN_MAX_STATIC_ALIGN_BYTES 0
#endif

#ifndef EIGEN_MAX_STATIC_ALIGN_BYTES

  #if EIGEN_COMP_GNUC && !(EIGEN_ARCH_i386_OR_x86_64 || EIGEN_ARCH_ARM_OR_ARM64 || EIGEN_ARCH_PPC || EIGEN_ARCH_IA64 || EIGEN_ARCH_MIPS)
  #define EIGEN_GCC_AND_ARCH_DOESNT_WANT_STACK_ALIGNMENT 1
  #elif EIGEN_ARCH_ARM_OR_ARM64 && EIGEN_COMP_GNUC_STRICT && EIGEN_GNUC_AT_MOST(4, 6)

  #define EIGEN_GCC_AND_ARCH_DOESNT_WANT_STACK_ALIGNMENT 1
  #else
  #define EIGEN_GCC_AND_ARCH_DOESNT_WANT_STACK_ALIGNMENT 0
  #endif

  #if !EIGEN_GCC_AND_ARCH_DOESNT_WANT_STACK_ALIGNMENT \
  && !EIGEN_GCC3_OR_OLDER \
  && !EIGEN_COMP_SUNCC \
  && !EIGEN_OS_QNX
    #define EIGEN_ARCH_WANTS_STACK_ALIGNMENT 1
  #else
    #define EIGEN_ARCH_WANTS_STACK_ALIGNMENT 0
  #endif

  #if EIGEN_ARCH_WANTS_STACK_ALIGNMENT
    #define EIGEN_MAX_STATIC_ALIGN_BYTES EIGEN_IDEAL_MAX_ALIGN_BYTES
  #else
    #define EIGEN_MAX_STATIC_ALIGN_BYTES 0
  #endif

#endif

#if defined(EIGEN_MAX_ALIGN_BYTES) && EIGEN_MAX_ALIGN_BYTES<EIGEN_MAX_STATIC_ALIGN_BYTES
#undef EIGEN_MAX_STATIC_ALIGN_BYTES
#define EIGEN_MAX_STATIC_ALIGN_BYTES EIGEN_MAX_ALIGN_BYTES
#endif

#if EIGEN_MAX_STATIC_ALIGN_BYTES==0 && !defined(EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT)
  #define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#endif

#define EIGEN_ALIGN8  EIGEN_ALIGN_TO_BOUNDARY(8)
#define EIGEN_ALIGN16 EIGEN_ALIGN_TO_BOUNDARY(16)
#define EIGEN_ALIGN32 EIGEN_ALIGN_TO_BOUNDARY(32)
#define EIGEN_ALIGN64 EIGEN_ALIGN_TO_BOUNDARY(64)
#if EIGEN_MAX_STATIC_ALIGN_BYTES>0
#define EIGEN_ALIGN_MAX EIGEN_ALIGN_TO_BOUNDARY(EIGEN_MAX_STATIC_ALIGN_BYTES)
#else
#define EIGEN_ALIGN_MAX
#endif

#if defined(EIGEN_DONT_ALIGN) && defined(EIGEN_MAX_ALIGN_BYTES) && EIGEN_MAX_ALIGN_BYTES>0
#error EIGEN_MAX_ALIGN_BYTES and EIGEN_DONT_ALIGN are both defined with EIGEN_MAX_ALIGN_BYTES!=0. Use EIGEN_MAX_ALIGN_BYTES=0 as a synonym of EIGEN_DONT_ALIGN.
#endif

#ifdef EIGEN_DONT_ALIGN
  #ifdef EIGEN_MAX_ALIGN_BYTES
    #undef EIGEN_MAX_ALIGN_BYTES
  #endif
  #define EIGEN_MAX_ALIGN_BYTES 0
#elif !defined(EIGEN_MAX_ALIGN_BYTES)
  #define EIGEN_MAX_ALIGN_BYTES EIGEN_IDEAL_MAX_ALIGN_BYTES
#endif

#if EIGEN_IDEAL_MAX_ALIGN_BYTES > EIGEN_MAX_ALIGN_BYTES
#define EIGEN_DEFAULT_ALIGN_BYTES EIGEN_IDEAL_MAX_ALIGN_BYTES
#else
#define EIGEN_DEFAULT_ALIGN_BYTES EIGEN_MAX_ALIGN_BYTES
#endif

#ifndef EIGEN_UNALIGNED_VECTORIZE
#define EIGEN_UNALIGNED_VECTORIZE 1
#endif

#if EIGEN_MAX_ALIGN_BYTES==0
  #ifndef EIGEN_DONT_VECTORIZE
    #define EIGEN_DONT_VECTORIZE
  #endif
#endif

#if EIGEN_COMP_MSVC
  #include <malloc.h> 
  #if (EIGEN_COMP_MSVC >= 1500) 
    
    #if (defined(_M_IX86_FP) && (_M_IX86_FP >= 2)) || EIGEN_ARCH_x86_64
      #define EIGEN_SSE2_ON_MSVC_2008_OR_LATER
    #endif
  #endif
#else
  #if (defined __SSE2__) && ( (!EIGEN_COMP_GNUC) || EIGEN_COMP_ICC || EIGEN_GNUC_AT_LEAST(4,2) )
    #define EIGEN_SSE2_ON_NON_MSVC_BUT_NOT_OLD_GCC
  #endif
#endif

#if !(defined(EIGEN_DONT_VECTORIZE) || defined(EIGEN_GPUCC))

  #if defined (EIGEN_SSE2_ON_NON_MSVC_BUT_NOT_OLD_GCC) || defined(EIGEN_SSE2_ON_MSVC_2008_OR_LATER)

    #define EIGEN_VECTORIZE
    #define EIGEN_VECTORIZE_SSE
    #define EIGEN_VECTORIZE_SSE2

    #ifdef __SSE3__
      #define EIGEN_VECTORIZE_SSE3
    #endif
    #ifdef __SSSE3__
      #define EIGEN_VECTORIZE_SSSE3
    #endif
    #ifdef __SSE4_1__
      #define EIGEN_VECTORIZE_SSE4_1
    #endif
    #ifdef __SSE4_2__
      #define EIGEN_VECTORIZE_SSE4_2
    #endif
    #ifdef __AVX__
      #ifndef EIGEN_USE_SYCL 
        #define EIGEN_VECTORIZE_AVX
      #endif
      #define EIGEN_VECTORIZE_SSE3
      #define EIGEN_VECTORIZE_SSSE3
      #define EIGEN_VECTORIZE_SSE4_1
      #define EIGEN_VECTORIZE_SSE4_2
    #endif
    #ifdef __AVX2__
      #ifndef EIGEN_USE_SYCL 
        #define EIGEN_VECTORIZE_AVX2
        #define EIGEN_VECTORIZE_AVX
      #endif
      #define EIGEN_VECTORIZE_SSE3
      #define EIGEN_VECTORIZE_SSSE3
      #define EIGEN_VECTORIZE_SSE4_1
      #define EIGEN_VECTORIZE_SSE4_2
    #endif
    #if defined(__FMA__) || (EIGEN_COMP_MSVC && defined(__AVX2__))

      #define EIGEN_VECTORIZE_FMA
    #endif
    #if defined(__AVX512F__)
      #ifndef EIGEN_VECTORIZE_FMA
      #if EIGEN_COMP_GNUC
      #error Please add -mfma to your compiler flags: compiling with -mavx512f alone without SSE/AVX FMA is not supported (bug 1638).
      #else
      #error Please enable FMA in your compiler flags (e.g. -mfma): compiling with AVX512 alone without SSE/AVX FMA is not supported (bug 1638).
      #endif
      #endif
      #ifndef EIGEN_USE_SYCL
        #define EIGEN_VECTORIZE_AVX512
        #define EIGEN_VECTORIZE_AVX2
        #define EIGEN_VECTORIZE_AVX
      #endif
      #define EIGEN_VECTORIZE_FMA
      #define EIGEN_VECTORIZE_SSE3
      #define EIGEN_VECTORIZE_SSSE3
      #define EIGEN_VECTORIZE_SSE4_1
      #define EIGEN_VECTORIZE_SSE4_2
      #ifndef EIGEN_USE_SYCL
        #ifdef __AVX512DQ__
          #define EIGEN_VECTORIZE_AVX512DQ
        #endif
        #ifdef __AVX512ER__
          #define EIGEN_VECTORIZE_AVX512ER
        #endif
        #ifdef __AVX512BF16__
          #define EIGEN_VECTORIZE_AVX512BF16
        #endif
      #endif
    #endif

    #if defined(__apple_build_version__) && (__apple_build_version__ == 11000033 ) && ( __MAC_OS_X_VERSION_MIN_REQUIRED == 101500 )

      #ifdef EIGEN_VECTORIZE_AVX
        #undef EIGEN_VECTORIZE_AVX
        #warning "Disabling AVX support: clang compiler shipped with XCode 11.[012] generates broken assembly with -macosx-version-min=10.15 and AVX enabled. "
        #ifdef EIGEN_VECTORIZE_AVX2
          #undef EIGEN_VECTORIZE_AVX2
        #endif
        #ifdef EIGEN_VECTORIZE_FMA
          #undef EIGEN_VECTORIZE_FMA
        #endif
        #ifdef EIGEN_VECTORIZE_AVX512
          #undef EIGEN_VECTORIZE_AVX512
        #endif
        #ifdef EIGEN_VECTORIZE_AVX512DQ
          #undef EIGEN_VECTORIZE_AVX512DQ
        #endif
        #ifdef EIGEN_VECTORIZE_AVX512ER
          #undef EIGEN_VECTORIZE_AVX512ER
        #endif
      #endif

    #endif

    extern "C" {

      #if EIGEN_COMP_ICC >= 1110
        #include <immintrin.h>
      #else
        #include <mmintrin.h>
        #include <emmintrin.h>
        #include <xmmintrin.h>
        #ifdef  EIGEN_VECTORIZE_SSE3
        #include <pmmintrin.h>
        #endif
        #ifdef EIGEN_VECTORIZE_SSSE3
        #include <tmmintrin.h>
        #endif
        #ifdef EIGEN_VECTORIZE_SSE4_1
        #include <smmintrin.h>
        #endif
        #ifdef EIGEN_VECTORIZE_SSE4_2
        #include <nmmintrin.h>
        #endif
        #if defined(EIGEN_VECTORIZE_AVX) || defined(EIGEN_VECTORIZE_AVX512)
        #include <immintrin.h>
        #endif
      #endif
    } 

  #elif defined __VSX__

    #define EIGEN_VECTORIZE
    #define EIGEN_VECTORIZE_VSX
    #include <altivec.h>

    #undef bool
    #undef vector
    #undef pixel

  #elif defined __ALTIVEC__

    #define EIGEN_VECTORIZE
    #define EIGEN_VECTORIZE_ALTIVEC
    #include <altivec.h>

    #undef bool
    #undef vector
    #undef pixel

  #elif ((defined  __ARM_NEON) || (defined __ARM_NEON__)) && !(defined EIGEN_ARM64_USE_SVE)

    #define EIGEN_VECTORIZE
    #define EIGEN_VECTORIZE_NEON
    #include <arm_neon.h>

  #elif (defined __ARM_FEATURE_SVE) && (defined EIGEN_ARM64_USE_SVE)

    #define EIGEN_VECTORIZE
    #define EIGEN_VECTORIZE_SVE
    #include <arm_sve.h>

    #if defined __ARM_FEATURE_SVE_BITS
      #define EIGEN_ARM64_SVE_VL __ARM_FEATURE_SVE_BITS
    #else
#error "Eigen requires a fixed SVE lector length but EIGEN_ARM64_SVE_VL is not set."
#endif

#elif (defined __s390x__ && defined __VEC__)

#define EIGEN_VECTORIZE
#define EIGEN_VECTORIZE_ZVECTOR
#include <vecintrin.h>

#elif defined __mips_msa

#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
#if defined(__LP64__)
#define EIGEN_MIPS_64
#else
#define EIGEN_MIPS_32
#endif
#define EIGEN_VECTORIZE
#define EIGEN_VECTORIZE_MSA
#include <msa.h>
#endif

#endif
#endif

#if defined(EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC)
  #include <arm_fp16.h>
#endif

#if defined(__F16C__) && (!defined(EIGEN_GPUCC) && (!defined(EIGEN_COMP_CLANG) || EIGEN_COMP_CLANG>=380))
  
  #define EIGEN_HAS_FP16_C

  #if defined(EIGEN_COMP_CLANG)

    #include <immintrin.h>
  #endif
#endif

#if defined EIGEN_CUDACC
  #define EIGEN_VECTORIZE_GPU
  #include <vector_types.h>
  #if EIGEN_CUDA_SDK_VER >= 70500
    #define EIGEN_HAS_CUDA_FP16
  #endif
#endif

#if defined(EIGEN_HAS_CUDA_FP16)
  #include <cuda_runtime_api.h>
  #include <cuda_fp16.h>
#endif

#if defined(EIGEN_HIPCC)
  #define EIGEN_VECTORIZE_GPU
  #include <hip/hip_vector_types.h>
  #define EIGEN_HAS_HIP_FP16
  #include <hip/hip_fp16.h>
#endif

namespace Eigen {

inline static const char *SimdInstructionSetsInUse(void) {
#if defined(EIGEN_VECTORIZE_AVX512)
  return "AVX512, FMA, AVX2, AVX, SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2";
#elif defined(EIGEN_VECTORIZE_AVX)
  return "AVX SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2";
#elif defined(EIGEN_VECTORIZE_SSE4_2)
  return "SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2";
#elif defined(EIGEN_VECTORIZE_SSE4_1)
  return "SSE, SSE2, SSE3, SSSE3, SSE4.1";
#elif defined(EIGEN_VECTORIZE_SSSE3)
  return "SSE, SSE2, SSE3, SSSE3";
#elif defined(EIGEN_VECTORIZE_SSE3)
  return "SSE, SSE2, SSE3";
#elif defined(EIGEN_VECTORIZE_SSE2)
  return "SSE, SSE2";
#elif defined(EIGEN_VECTORIZE_ALTIVEC)
  return "AltiVec";
#elif defined(EIGEN_VECTORIZE_VSX)
  return "VSX";
#elif defined(EIGEN_VECTORIZE_NEON)
  return "ARM NEON";
#elif defined(EIGEN_VECTORIZE_SVE)
  return "ARM SVE";
#elif defined(EIGEN_VECTORIZE_ZVECTOR)
  return "S390X ZVECTOR";
#elif defined(EIGEN_VECTORIZE_MSA)
  return "MIPS MSA";
#else
  return "None";
#endif
}

} 

#endif 
