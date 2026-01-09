#ifndef Model_hpp
#define Model_hpp

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "History.hpp"
#include "Threads.hpp"
class CondLikeJobMngr;
class History;
class MetaData;
class Node;
class RandomVariable;
class RateMatrix;
class Samples;
class ThreadPool;
class TransitionProbabilitiesMngr;
class Tree;

#define MAX_NUM_CHANGES     10

class Model {

    public:
                                        Model(void) = delete;
                                        Model(RandomVariable* r, Tree* tp, MetaData* md, ThreadPool* thp, CondLikeJobMngr* mngr);
                                       ~Model(void);
        void                            accept(void);
        int***                          getIntervalTransitions(void) { return intervalTransitions; }
        int                             getNumIntervals(void) { return numIntervals; }
        int                             getNumStates(void) { return (int)numStates; }
        std::vector<double>&            getPi(void) { return pi[1]; }
        std::vector<double>&            getR(void) { return r[1]; }
        RateMatrix*                     getRateMatrix(void) { return q[activeRateMatrix]; }
        double                          getSubstitutionRate(void) { return substitutionRate[1]; }
        Tree*                           getTree(void) { return tree; }
        std::string                     getUpdateType(void) { return updateType; }
        double                          lnLikelihood(void);
        double                          lnPriorProbability(void);
        void                            map(void);
        void                            reject(void);
        double                          update(void);
    
    private:
        void                            assignNodeStates(RandomVariable* rng);
        void                            checkConditionalLikelihoods(void);
        void                            deleteHistories(void);
        void                            incrementDwellTimes(History* h, Node* p, Change* change);
        void                            initializeConditionalLikelihoods(void);
        void                            initializeHistories(void);
        void                            initializeMatrixPowers(size_t num);
        void                            initializeParameters(Tree* tp, RateMatrix* m);
        bool                            isValidSimplex(const std::vector<double>& x, double eps=1e-14);
        void                            printMatrixPowers(void);
        std::vector<Samples*>           readParameterFile(std::string fn);
        int                             sampleHistoriesUsingRejectionSamplign(RandomVariable* rng);
        int                             sampleHistoriesUsingUniformization(RandomVariable* rng);
        void                            switchActiveRateMatrix(void);
        double                          updateSubstitutionRate(void);
        double                          updatePi(void);
        double                          updateR(void);
        void                            updateRateMatrix(void);
        double                          updateSimplexTransfer(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal);
        double                          updateSimplex(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal);
        double                          updateSimplex(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, size_t k, double minVal);
        double                          updateSimplexALRMVN(std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal);
        CondLikeJobMngr*                clManager;
        Tree*                           tree;
        MetaData*                       metaData;
        RandomVariable*                 rng;
        RateMatrix*                     q[2];
        RateMatrix*                     uniformizedRateMatrix;
        TransitionProbabilitiesMngr*    tiMngr;
        ThreadPool*                     threadPool;
        size_t                          numStates;
        double*                         condLikes;
        double**                        intervalDwellTimes;
        int***                          intervalTransitions;
        double                          substitutionRate[2];
        std::vector<double>             pi[2];
        std::vector<double>             r[2];
        std::vector<RateMatrix*>        matrixPowers;
        int                             activeRateMatrix;
        int                             numIntervals;
        double                          nFactorial[MAX_NUM_CHANGES];
        std::string                     updateType;
};


#endif
