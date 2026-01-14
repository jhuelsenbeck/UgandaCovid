#ifndef TransitionProbabilities_hpp
#define TransitionProbabilities_hpp

#include <atomic>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include "Container.hpp"
#include "MathCacheAccelerated.hpp"
#include "Node.hpp"
#include "Threads.hpp"
#include "Tree.hpp"

// Legacy functions - kept for compatibility
double factorial(int x);
int logBase2Plus1(double x);
int setQvalue(double tol);

// Substitution model used for transition probabilities
//
//  jc69        : equal stationary frequencies
//  f81         : stationary frequencies are parameters (estimated)
//  custom_f81  : stationary frequencies are fixed (set at initialization)
//  gtr         : general time-reversible model (pi + exchangeabilities)
//
// NOTE: For f81/custom_f81, we compute P analytically (rank-1 update),
//       which is substantially faster than a full matrix exponential.
//       For gtr we use the existing Padé[13/13] + scaling/squaring.

enum SubModel { jc69, f81, custom_f81, gtr };


// -------------------------------------------------------------------
// TransitionProbabilitiesTask
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
            pi         = nullptr;
            r          = nullptr;
            k          = 1.0;
            ugandaIdx  = 0;
            modelType  = gtr;
        }

        #if defined(__GNUC__) || defined(__clang__)
        __attribute__((noinline))
        #endif
        void init(int i,
                  int n,
                  double v,
                  DoubleMatrix* q,
                  DoubleMatrix* p,
                  const std::vector<double>* piPtr,
                  const std::vector<double>* rPtr,
                  const double kappa,
                  size_t ugi) {

            std::cerr << "    init() entered, this=" << (void*)this << std::endl;
            std::cerr << "    params: i=" << i << " n=" << n << " v=" << v << std::endl;
            std::cerr << "    params: q=" << (void*)q << " p=" << (void*)p << std::endl;
            std::cerr << "    params: piPtr=" << (void*)piPtr << " rPtr=" << (void*)rPtr << std::endl;
            std::cerr << "    params: kappa=" << kappa << " ugi=" << ugi << std::endl;
            
            std::cerr << "    assigning taskId..." << std::endl;
            this->taskId     = i;
            std::cerr << "    assigning numStates..." << std::endl;
            this->numStates  = n;
            std::cerr << "    assigning brlen..." << std::endl;
            this->brlen      = v;
            std::cerr << "    assigning Q..." << std::endl;
            this->Q          = q;
            std::cerr << "    assigning P..." << std::endl;
            this->P          = p;
            std::cerr << "    assigning pi..." << std::endl;
            this->pi         = piPtr;
            std::cerr << "    assigning r..." << std::endl;
            this->r          = rPtr;
            std::cerr << "    assigning k..." << std::endl;
            this->k          = kappa;
            std::cerr << "    assigning ugandaIdx..." << std::endl;
            this->ugandaIdx  = ugi;
            
            std::cerr << "    init() done" << std::endl;
            // Memory fence: ensure all writes above are visible to other threads
            // before this task is pushed to the thread pool
            std::atomic_thread_fence(std::memory_order_release);
        }

        void run(void) override {

            // Memory fence: ensure all writes from init() are visible
            std::atomic_thread_fence(std::memory_order_acquire);

            // DEBUG: Print state at start of run
            std::cerr << "DEBUG run() starting: taskId=" << taskId 
                      << " numStates=" << numStates 
                      << " modelType=" << modelType 
                      << " pi=" << (void*)pi 
                      << " P=" << (void*)P << std::endl;

            if (modelType == jc69)
                {
                // k-state Jukes-Cantor, scaled so that the average rate is 1.0.
                // With q_ij = 1/(k-1), the non-zero eigenvalue is -k/(k-1).
                const double k = (double)numStates;
                const double e = std::exp(-(k / (k - 1.0)) * brlen);

                const double a = 1.0 / k;
                const double diag = a + (1.0 - a) * e;
                const double off  = a - a * e;

                // fill entire matrix (compiler should do something clever here)
                double* p = P->begin();
                double* end = p + (numStates * numStates);
                while (p < end)
                    *p++ = off;

                // walk the diagonal
                p = P->begin();
                size_t stride = numStates + 1;
                for (size_t i=0; i<numStates; i++)
                    {
                    *p = diag;
                    p += stride;
                    }
                }
            else if (modelType == f81)
                {
                // F81: P_ij(t) = pi_j + (delta_ij - pi_j) * exp(-mu t)
                // where mu is chosen so that the average rate is 1.0:
                //   average rate = mu * (1 - sum_i pi_i^2)  => mu = 1/(1 - sum pi^2)
                double s = 0.0;
                for (size_t i=0; i<numStates; i++)
                    s += (*pi)[i] * (*pi)[i];
                double mu = 1.0 / (1.0 - s);

                const double e = std::exp(-mu * brlen);

                for (int i=0; i<numStates; i++)
                    {
                    for (int j=0; j<numStates; j++)
                        {
                        const double pij = (*pi)[j];
                        const double delta = (i == j ? 1.0 : 0.0);
                        (*P)(i,j) = pij + (delta - pij) * e;
                        }
                    }
                }
            else if (modelType == custom_f81)
                {
                // custom F81
                double p = (*pi)[ugandaIdx];
                double q = 1.0 - p;
                double qInv = 1.0 / q;
                double s = 0.0;
                for (int i = 0; i < numStates; i++)
                    s += (*pi)[i] * (*pi)[i];
                double r = 1.0 - s + 2.0 * (k - 1.0) * p * q;
                double lambda1 = (-1.0 + (1.0 - k) * p) / r;
                double lambda2 = -k / r;
                double exp1 = exp(lambda1 * brlen);
                double exp2 = exp(lambda2 * brlen);

                // Precompute common terms
                double one_minus_exp2 = 1.0 - exp2;
                double q_exp2 = q * exp2;
                double p_qInv = p * qInv;

                // Case 1: i != ugandaIdx, j != ugandaIdx, off-diagonal (i != j)
                for (int i = 0; i < numStates; i++) 
                    {
                    if (i == ugandaIdx) 
                        continue;
                    for (int j = 0; j < numStates; j++) 
                        {
                        if (j == ugandaIdx || i == j) 
                            continue;
                        (*P)(i, j) = (*pi)[j] - (*pi)[j] * qInv * exp1 + (*pi)[j] * exp2 * p_qInv;
                        }
                    }

                // Case 2: i != ugandaIdx, j != ugandaIdx, diagonal (i == j)
                for (int i = 0; i < numStates; i++) 
                    {
                    if (i == ugandaIdx) 
                        continue;
                    (*P)(i, i) = (*pi)[i] + exp1 * (1.0 - (*pi)[i] * qInv) + (*pi)[i] * exp2 * p_qInv;
                    }

                // Case 3: i != ugandaIdx, j == ugandaIdx
                double pi_k_one_minus_exp2 = (*pi)[ugandaIdx] * one_minus_exp2;
                for (int i = 0; i < numStates; i++) 
                    {
                    if (i == ugandaIdx) 
                        continue;
                    (*P)(i, ugandaIdx) = pi_k_one_minus_exp2;
                    }

                // Case 4: i == ugandaIdx, j != ugandaIdx
                for (int j = 0; j < numStates; j++) 
                    {
                    if (j == ugandaIdx) 
                        continue;
                    (*P)(ugandaIdx, j) = (*pi)[j] * one_minus_exp2;
                    }

                // Case 5: i == ugandaIdx, j == ugandaIdx (diagonal)
                (*P)(ugandaIdx, ugandaIdx) = (*pi)[ugandaIdx] + q_exp2;
                }
            else // gtr
                {
                MathCacheAccelerated& cache = getCache();
                cache.computeMatrixExponential(*Q, brlen, *P);
                }
        }

        void setModelType(SubModel mt) { modelType = mt; }

    private:
        int                         taskId;
        int                         numStates;
        double                      brlen;
        DoubleMatrix*               Q;
        DoubleMatrix*               P;
        const std::vector<double>*  pi;
        const std::vector<double>*  r;
        double                      k;
        size_t                      ugandaIdx;
        SubModel                    modelType;
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
