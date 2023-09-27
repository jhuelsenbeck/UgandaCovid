#ifndef Model_hpp
#define Model_hpp

#include <cmath>
#include <iostream>
#include "jph.hpp"
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
        real                          lnLikelihood(void);
        real                          lnPriorProbability(void);
    
    private:
        void                            initializeConditionalLikelihoods(void);
        Tree*                           tree;
        RateMatrix*                     q[2];
        TransitionProbabilitiesMngr*    tiMngr;
        ThreadPool*                     threadPool;
        real*                         condLikes;
        int                             activeRateMatrix;
        int                             activeSubstitutionRate;
        int                             activePi;
        int                             activeR;
        real                            substitutionRate[2];
        std::vector<real>             pi[2];
        std::vector<real>             r[2];
        size_t                          numStates;
};

class ConditionalLikelihoodTask : public ThreadTask {

    public:
        ConditionalLikelihoodTask(void) {
        
            numStates     = 0;
            clBegin       = NULL;
            tiBegin       = NULL;
        }
        
        void init(int n, real* l, real* b, real* p) {
        
            numStates     = n;
            clBegin       = l;
            tiBegin       = p;
        }

        void condLike(MathCache& cache) {
        
            real* P = tiBegin;
            real* clEnd = clBegin + numStates;
            auto lnSum  = cache.pushArray(numStates);
            real* sStart = lnSum->begin();
            real* sEnd = lnSum->end();
            
            for (real* s=sStart; s != sEnd; s++)
                {
                real sum = 0.0;
                for (real* L=clBegin; L != clEnd; L++)
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
        real*         clBegin;
        real*         tiBegin;
};


#endif
