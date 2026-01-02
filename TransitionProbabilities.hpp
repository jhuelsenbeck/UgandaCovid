#ifndef TransitionProbabilities_hpp
#define TransitionProbabilities_hpp

#include <cmath>
#include "Container.hpp"
#include "MathCacheAccelerated.hpp"
#include "Node.hpp"
#include "Threads.hpp"
#include "Tree.hpp"

// Legacy functions - kept for compatibility
double factorial(int x);
int logBase2Plus1(double x);
int setQvalue(double tol);


// -------------------------------------------------------------------
// TransitionProbabilitiesTask
// -------------------------------------------------------------------
// Computes P = exp(Q * brlen) using Padé[13/13] approximant with
// scaling and squaring. Uses Apple's Accelerate BLAS/LAPACK for
// high-performance matrix operations.
//
// Key optimizations:
// - Thread-local MathCacheAccelerated with persistent working buffers
// - cblas_dgemm for matrix multiplication (vectorized, cache-blocked)
// - dgetrf/dgetrs for LU solve (optimized pivoting)
// - No heap allocations in hot path after first use
// -------------------------------------------------------------------

class TransitionProbabilitiesTask : public ThreadTask {

    public:
    
        // Thread-local accelerated cache - each worker gets its own
        // Persists across task executions within the same thread
        static MathCacheAccelerated& getCache() {
            thread_local MathCacheAccelerated cache;
            return cache;
        }
    
        TransitionProbabilitiesTask(void) {
        
            taskId     = 0;
            numStates  = 0;
            brlen      = 0.0;
            Q          = nullptr;
            P          = nullptr;
        }
        
        void init(int i, int n, double v, DoubleMatrix* q, DoubleMatrix* p) {
        
            taskId     = i;
            numStates  = n;
            brlen      = v;
            Q          = q;
            P          = p;
        }

        void run(void) override {
            
            MathCacheAccelerated& cache = getCache();
            cache.computeMatrixExponential(*Q, brlen, *P);
        }

    private:
        int             taskId;
        int             numStates;
        double          brlen;
        DoubleMatrix*     Q;
        DoubleMatrix*     P;
};



class TransitionProbabilities : public DoubleMatrix {

    public:
                    TransitionProbabilities(void) = delete;
                    TransitionProbabilities(const TransitionProbabilities& m);
                    TransitionProbabilities(size_t d);
        virtual    ~TransitionProbabilities(void);
        size_t      dim(void) { return numRows; }
        int         getBrlen(void) { return brlen; }
        double      getCalculatedBrlen(void) { return calculatedBrlen; }
        void        setBrlen(int x) { brlen = x; }
        void        setCalculatedBrlen(double x) { calculatedBrlen = x; }
    
    private:
        int         brlen;
        double      calculatedBrlen;

    friend std::ostream& operator<<(std::ostream& os, const TransitionProbabilities& m);
};

inline std::ostream& operator<<(std::ostream& os, const TransitionProbabilities& m) {

    os << std::fixed << std::setprecision(3) << std::scientific;
    for (size_t i=0; i<m.numRows; i++)
        {
        for (size_t j=0; j<m.numCols; j++)
            {
            if (m(i,j) > 0)
                os << " ";
            os << m(i,j) << " ";
            }
        os << std::endl;
        }
    return os;
}

#endif
