#ifndef Model_hpp
#define Model_hpp

#include <iostream>
#include "Threads.hpp"
class RateMatrix;
class ThreadPool;
class TransitionProbabilitiesMngr;
class Tree;


class Model {

    public:
                                        Model(void) = delete;
                                        Model(Tree* tp, RateMatrix* m, ThreadPool* thp);
                                       ~Model(void);
        RateMatrix*                     getRateMatrix(void) { return q[activeRateMatrix]; }
        double                          lnLikelihood(void);
        double                          lnPriorProbability(void);
    
    private:
        void                            initializeConditionalLikelihoods(void);
        Tree*                           tree;
        RateMatrix*                     q[2];
        TransitionProbabilitiesMngr*    tiMngr;
        ThreadPool*                     threadPool;
        int                             activeRateMatrix;
        int                             activeSubstitutionRate;
        int                             activePi;
        int                             activeR;
        double                          substitutionRate[2];
        std::vector<double>             pi[2];
        std::vector<double>             r[2];
        size_t                          numStates;
};

class ConditionalLikelihoodTask : public ThreadTask {

    public:
        ConditionalLikelihoodTask(void) {
        
            numStates     = 0;
            clBegin       = NULL;
            tiBegin       = NULL;
        }
        
        void init(int n, double* l, double* b, double* p) {
        
            numStates     = n;
            clBegin       = l;
            tiBegin       = p;
        }

        void condLike(MathCache& cache) {
        
            double* P = tiBegin;
            double* clEnd = clBegin + numStates;
            auto lnSum  = cache.pushArray(numStates);
            double* sStart = lnSum->begin();
            double* sEnd = lnSum->end();
            
            for (double* s=sStart; s != sEnd; s++)
                {
                double sum = 0.0;
                for (double* L=clBegin; L != clEnd; L++)
                    {
                    sum += (*P) * (*L);
                    P++;
                    }
                *s = log(sum);
                }
        }

        virtual void Run(MathCache& cache) {
            
            //std::cout << "Running task " << taskId << "  << std::endl;
            condLike(cache);
            }

    private:
        int             numStates;
        double*         clBegin;
        double*         tiBegin;
};


#endif
