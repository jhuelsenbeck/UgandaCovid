#ifndef Model_hpp
#define Model_hpp

#include <cmath>
#include <iostream>
#include <string>
#include "Threads.hpp"
class CondLikeJobMngr;
class RateMatrix;
class ThreadPool;
class TransitionProbabilitiesMngr;
class Tree;


class Model {

    public:
                                        Model(void) = delete;
                                        Model(Tree* tp, RateMatrix* m, ThreadPool* thp, CondLikeJobMngr* mngr);
                                       ~Model(void);
        void                            accept(void);
        std::vector<double>&            getPi(void) { return pi[1]; }
        std::vector<double>&            getR(void) { return r[1]; }
        RateMatrix*                     getRateMatrix(void) { return q[activeRateMatrix]; }
        double                          getSubstitutionRate(void) { return substitutionRate[1]; }
        std::string                     getUpdateType(void) { return updateType; }
        double                          lnLikelihood(void);
        double                          lnPriorProbability(void);
        void                            map(void);
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
        CondLikeJobMngr*                clManager;
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


#endif
