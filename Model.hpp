#ifndef Model_hpp
#define Model_hpp

#include <cmath>
#include <iostream>
#include <string>
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
        void                            accept(void);
        std::vector<double>&            getPi(void) { return pi[1]; }
        std::vector<double>&            getR(void) { return r[1]; }
        RateMatrix*                     getRateMatrix(void) { return q[activeRateMatrix]; }
        double                          getSubstitutionRate(void) { return substitutionRate[1]; }
        std::string                     getUpdateType(void) { return updateType; }
        double                          lnLikelihood(void);
        double                          lnPriorProbability(void);
        void                            reject(void);
        double                          update(void);
    
    private:
        void                            checkConditionalLikelihoods(void);
        void                            initializeConditionalLikelihoods(void);
        void                            switchActiveRateMatrix(void);
        double                          updateSubstitutionRate(void);
        double                          updatePi(void);
        double                          updateR(void);
        void                            updateRateMatrix(void);
        double                          updateSimplex(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0);
        double                          updateSimplex(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, int k);
        Tree*                           tree;
        int                             activeRateMatrix;
        RateMatrix*                     q[2];
        TransitionProbabilitiesMngr*    tiMngr;
        ThreadPool*                     threadPool;
        double*                         condLikes;
        double                          substitutionRate[2];
        std::vector<double>             pi[2];
        std::vector<double>             r[2];
        size_t                          numStates;
        std::string                     updateType;
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
