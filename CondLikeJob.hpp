#ifndef CondLikeJob_hpp
#define CondLikeJob_hpp

#include <mutex>
#include <vector>
#include "Threads.hpp"

// -------------------------------------------------------------------
// Platform detection for BLAS and vectorized math
// -------------------------------------------------------------------
#if defined(__APPLE__)
    #define CONDLIKE_USE_ACCELERATE 1
    #ifndef ACCELERATE_NEW_LAPACK
    #define ACCELERATE_NEW_LAPACK
    #endif
    #include <Accelerate/Accelerate.h>
#elif defined(USE_OPENBLAS) || defined(USE_CBLAS)
    #define CONDLIKE_USE_CBLAS 1
    extern "C" {
        #include <cblas.h>
    }
#endif

class CondLikeJobMngr;
class Node;
class Tree;



// -------------------------------------------------------------------
// CondLikeJob
// -------------------------------------------------------------------
// Computes conditional likelihoods for a set of nodes in the tree.
// This is the core computation in phylogenetic likelihood calculation.
//
// Key optimizations:
// - BLAS dgemv for matrix-vector multiply (P * CL)
// - Vectorized log via vvlog (Accelerate) or manual SIMD
// - Pre-allocated working buffers
// - Minimized memory passes
// -------------------------------------------------------------------

class CondLikeJob : public ThreadTask {

    public:
                                    CondLikeJob(void) = delete;
                                    CondLikeJob(ThreadPool* tp, int na);
                                   ~CondLikeJob(void);
        void                        addDependency(CondLikeJob* ) { numDependencies++; }
        void                        addNode(Node* p) { nodes.push_back(p); }
        int                         getJobId(void) { return jobId; }
        std::vector<Node*>&         getJobNodes(void) { return nodes; }
        double                      getLnScaler(void) { return lnScaler; }
        int                         getNumDependencies(void) { return numDependencies; }
        int                         numNodesInJob(void) { return (int)nodes.size(); }
        void                        print(void);
        void                        resolveDependency(void);
        virtual void                run(void);
        void                        setJobId(int x) { jobId = x; }
        void                        setNumDepencies(int x) { numDependencies = x; }
        void                        setNumResolvedDependencies(int x) { numResolvedDependencies = x; }
    
    private:
        void                        conditionalLikelihood(void);
        void                        conditionalLikelihoodOptimized(void);
        void                        conditionalLikelihoodPortable(void);
        
                                    // working buffers (pre-allocated, reused)
        double*                     clSum;          // accumulated log-likelihoods
        double*                     clSumEnd;
        double*                     tempResult;     // result of P * CL (for BLAS)
        
        // Job management
        ThreadPool*                 threadPool;
        std::mutex                  mtx;
        double                      lnScaler;
        int                         jobId;
        std::vector<Node*>          nodes;
        int                         numDependencies;
        int                         numResolvedDependencies;
        int                         numStates;
};

#endif
