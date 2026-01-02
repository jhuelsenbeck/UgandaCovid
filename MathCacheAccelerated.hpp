#ifndef MathCacheAccelerated_hpp
#define MathCacheAccelerated_hpp

#include "Container.hpp"
#include "MathCache.hpp"

// -------------------------------------------------------------------
// Platform detection for BLAS/LAPACK
// -------------------------------------------------------------------
// - macOS: Use Apple's Accelerate framework (highly optimized for Apple Silicon)
// - Linux: Use OpenBLAS or system BLAS if available
// - Other: Fall back to pure C++ implementation
// -------------------------------------------------------------------

#if defined(__APPLE__)
    #define USE_ACCELERATE 1
    #ifndef ACCELERATE_NEW_LAPACK
    #define ACCELERATE_NEW_LAPACK
    #endif
    #include <Accelerate/Accelerate.h>
    typedef __LAPACK_int lapack_int;
#elif defined(USE_OPENBLAS) || defined(USE_BLAS)
    #define USE_CBLAS 1
    extern "C" {
        #include <cblas.h>
    }
    // OpenBLAS LAPACK declarations
    extern "C" {
        void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
        void dgetrs_(char* trans, int* n, int* nrhs, double* a, int* lda,
                     int* ipiv, double* b, int* ldb, int* info);
    }
    typedef int lapack_int;
#else
    #define USE_PORTABLE 1
    typedef int lapack_int;
#endif

#include <vector>



// -------------------------------------------------------------------
// MathCacheAccelerated
// -------------------------------------------------------------------
// Extends MathCache with optimized matrix operations.
// Uses Padé[13/13] approximant with scaling and squaring for
// matrix exponentials - the standard high-accuracy method.
//
// Key optimizations:
// - BLAS for matrix multiplication when available
// - LAPACK for LU solve when available
// - Pre-allocated working buffers to avoid heap fragmentation
// -------------------------------------------------------------------

class MathCacheAccelerated : public MathCache {

    public:
                            MathCacheAccelerated(void);
                           ~MathCacheAccelerated(void);
        
        // Compute P = exp(Q * t) using Padé[13/13] with scaling and squaring
        void                computeMatrixExponential(const DoubleMatrix& Q, double t, DoubleMatrix& P);
        
        // Matrix multiply: C = A * B (uses BLAS when available)
        void                optimizedMultiply(double* A, double* B, double* C, int n);
        
        // Compute infinity norm (max row sum of absolute values)
        double              infinityNorm(double* A, size_t n);
        
    private:
        // Working buffers for matrix exponential
        std::vector<double> workA, workA2, workA4, workA6;
        std::vector<double> workU, workV, workTemp1, workTemp2;
        std::vector<double> workTemp1T, workTemp2T;
        std::vector<lapack_int> workPivot;
        size_t              workSize;
        
        void                ensureWorkBuffers(size_t n);
        
        // Portable matrix multiply fallback
        void                portableMultiply(double* A, double* B, double* C, int n);
        
        // Portable LU solve fallback  
        void                portableLUSolve(double* A, double* B, int n);
};

#endif
