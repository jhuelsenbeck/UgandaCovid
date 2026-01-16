#ifndef Model_hpp
#define Model_hpp

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "History.hpp"
#include "Threads.hpp"
#include "TransitionProbabilities.hpp"
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
        void                            checkPoint(std::string fileName);
        int***                          getIntervalTransitions(void) { return intervalTransitions; }
        int                             getNumIntervals(void) { return numIntervals; }
        int                             getNumStates(void) { return (int)numStates; }
        std::vector<double>&            getPi(void) { return pi[1]; }
        std::vector<double>&            getR(void) { return r[1]; }
        double                          getKappa(void) { return kappa[1]; }
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
        void                            loadCheckpont(void);
        void                            printMatrixPowers(void);
        std::vector<Samples*>           readParameterFile(std::string fn);
        int                             sampleHistoriesUsingRejectionSamplign(RandomVariable* rng);
        int                             sampleHistoriesUsingUniformization(RandomVariable* rng);
        void                            switchActiveRateMatrix(void);
        double                          updateSubstitutionRate(void);
        double                          updatePi(void);
        double                          updateR(void);
        double                          updateKappa(void);
        double                          updateKappaLockdown(void);
        double                          updateKappaNoLockdown(void);
        void                            updateRateMatrix(void);
        double                          updateSimplex(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal);
        double                          updateSimplex(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, size_t k, double minVal);
        double                          updateSimplexTransfer(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal);
        double                          updateSimplexALRMVN(std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal);
        double                          updateSimplexALRMVN(const std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal, size_t blockSize);
        double                          updateSimplexPrior(const std::vector<double>& oldVec, std::vector<double>& newVec, std::vector<double>& alpha, double minVal);
        CondLikeJobMngr*                clManager;
        Tree*                           tree;
        MetaData*                       metaData;
        RandomVariable*                 rng;
        RateMatrix*                     q[2];
        RateMatrix*                     uniformizedRateMatrix;
        TransitionProbabilitiesMngr*    tiMngr;
        ThreadPool*                     threadPool;
        SubModel                        modelType;
        bool                            variableUgandaRate;
        size_t                          numStates;
        double*                         condLikes;
        double**                        intervalDwellTimes;
        int***                          intervalTransitions;
        double                          substitutionRate[2];
        std::vector<double>             pi[2];
        std::vector<double>             r[2];
        std::vector<double>             alphaPi;
        std::vector<double>             alphaR;
        double                          kappa[2];
        double                          kappaLockdown[2];
        double                          kappaNoLockdown[2];
        size_t                          ugandaIdx;
        std::vector<RateMatrix*>        matrixPowers;
        int                             activeRateMatrix;
        int                             numIntervals;
        double                          nFactorial[MAX_NUM_CHANGES];
        std::string                     updateType;
};


#endif
