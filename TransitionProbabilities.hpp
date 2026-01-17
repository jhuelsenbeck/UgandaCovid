#ifndef TransitionProbabilities_hpp
#define TransitionProbabilities_hpp

#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>
#include "Container.hpp"
#include "MathCacheAccelerated.hpp"
#include "Node.hpp"
#include "Threads.hpp"
#include "Tree.hpp"

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
            k2         = 1.0;
            ugandaIdx  = 0;
            modelType  = gtr;
        }

        // Declaration only - definition in TransitionProbabilities.cpp
        void init(int i,
                  int n,
                  double v,
                  double v2,
                  double v3,
                  DoubleMatrix* q,
                  DoubleMatrix* p,
                  const std::vector<double>* piPtr,
                  const std::vector<double>* rPtr,
                  const double kappa,
                  const double kappa2,
                  size_t ugi);

        void run(void) override {

            if (modelType == jc69)
                {
                tiProbsJC69();
                }
            else if (modelType == f81)
                {
                tiProbsF81();
                }
            else if (modelType == custom_f81)
                {
                tiProbsF81_Custom();
                }
            else
                {
                MathCacheAccelerated& cache = getCache();
                cache.computeMatrixExponential(*Q, brlen, *P);
                }
        }
        void tiProbsJC69(void);
        void tiProbsF81(void);
        void tiProbsF81_Custom(void);
        void tiProbsF81_Custom_Variable(void);
        void setModelType(SubModel mt) { modelType = mt; }
        void setCustomModelParms(double p, double q, double lambda1, double lambda2, double lambda12, double lambda22);

    private:
        int                         taskId;
        int                         numStates;
        double                      brlen;
        double                      brlen2;
        double                      brlen3;
        DoubleMatrix*               Q;
        DoubleMatrix*               P;
        const std::vector<double>*  pi;
        const std::vector<double>*  r;
        double                      k;
        double                      k2;
        double                      p, q, lambda1, lambda2, lambda12, lambda22;
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
