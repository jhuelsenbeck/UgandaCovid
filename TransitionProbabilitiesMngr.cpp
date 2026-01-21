#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include "GPUMatrixExponentialBatch.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "RateMatrix.hpp"
#include "Threads.hpp"
#include "TransitionProbabilitiesMngr.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#define MIN_BRLEN 1.0e-5



TransitionProbabilitiesMngr::TransitionProbabilitiesMngr(Model* m, Tree* t, size_t d, ThreadPool* tp, size_t ugi) :
    modelPtr(m), dim(d), threadPool(tp), ugandaIdx(ugi) {

    // set the transition-probability model from user settings
    UserSettings& settings = UserSettings::getUserSettings();
    switch (settings.getLikelihoodModel())
        {
        case LikelihoodModel::JC69:       modelType = jc69;       break;
        case LikelihoodModel::F81:        modelType = f81;        break;
        case LikelihoodModel::CUSTOM_F81: modelType = custom_f81; break;
        case LikelihoodModel::GTR:        modelType = gtr;        break;
        }
    isUgandaRateVariable = settings.isUgandaRateVariable();

    // initialize GPU batcher singleton
    gpuBatcher = &GPUMatrixExponentialBatch::getInstance();
    
    // default settings
    computeBackend = TiProbComputeBackend::Auto;
    batchThreshold = 8;  // Use batched approach for >= 8 matrices
    
    // for large matrices, the batched approach is always better
    if (dim >= 64)
        batchThreshold = 4;
    
    // allocate matrices for the unique branch lengths
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (int i=0, n=(int)dpSeq.size(); i<n; i++)
        {
        Node* p = dpSeq[i];
        if (p == t->getRoot())
            continue;
            
        int v = p->getBrlen();
        
        if (modelType == custom_f81 && isUgandaRateVariable == true)
            {
            double t0=0.0, t1=0.0, t2=0.0;
            if (p->getAncestor() != nullptr)
                modelPtr->computeDwellTimes(p->getTime(), p->getAncestor()->getTime(), &t0, &t1, &t2);
            key = std::make_tuple((int)t0, (int)t1, (int)t2);
            }
        else 
            {
            key = std::make_tuple(v,0,0);
            }
            
        ti_map::iterator it = tiMap.find(key);
        if (it == tiMap.end())
            {
            TransitionProbabilities* ti = new TransitionProbabilities(dim);
            ti->setBrlen(v);
            tiMap.insert( std::make_pair(key,ti) );
            }
        }
            
    // pre-allocate vectors for batch operations
    batchBranchLengths.reserve(tiMap.size());
    batchOutputs.reserve(tiMap.size());

    // preallocate the transition probability tasks for threaded calculations
    tasks = new TransitionProbabilitiesTask[tiMap.size()];
    for (size_t i=0; i<tiMap.size(); i++)
        tasks[i].setModelType(modelType);
    
    std::cout << "     Tree has " << tiMap.size() << " unique branch lengths" << std::endl;
    std::cout << "     GPU batch: " << (gpuBatcher->isAvailable() ? gpuBatcher->getDeviceName() : "Not available") << std::endl;
}

TransitionProbabilitiesMngr::~TransitionProbabilitiesMngr(void) {

    for (auto p : tiMap)
        {
        TransitionProbabilities* ti = p.second;
        delete ti;
        }
    tiMap.clear();
    delete [] tasks;
}

void TransitionProbabilitiesMngr::checkTiProbs(void) {

    for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
        {
        TransitionProbabilities* p = it->second;
        size_t n = p->dim();
        std::cout << "matrix " << std::get<0>(it->first) << " " << std::get<1>(it->first) << " " << std::get<2>(it->first) << std::endl;
        
        bool hasIssues = false;
        for (size_t i=0; i<n; i++)
            {
            double sum = 0.0;
            for (size_t j=0; j<n; j++)
                {
                double x = (*p)(i,j);
                sum += x;
                if (x < 0.0) {
                    Msg::warning("Negative transition probability!!!!!");
                    hasIssues = true;
                }
                if (std::isnan(x)) {
                    Msg::warning("NaN transition probability!!!!!");
                    hasIssues = true;
                }
                if (std::isinf(x)) {
                    Msg::warning("Infinite transition probability!!!!!");
                    hasIssues = true;
                }
                }
            std::cout << "Row " << i << " sum: " << sum << std::endl;
            if (std::abs(sum - 1.0) > 1e-10) {
                std::cout << "WARNING: Row " << i << " sum deviates significantly from 1.0: " << sum << std::endl;
                hasIssues = true;
            }
            }
            
        if (hasIssues) {
            std::cout << "Problematic transition matrix:" << std::endl;
            for (size_t i=0; i<n; i++) {
                for (size_t j=0; j<n; j++) {
                    std::cout << std::setw(12) << (*p)(i,j) << " ";
                }
                std::cout << std::endl;
            }
        }
        }
}

TransitionProbabilities* TransitionProbabilitiesMngr::getTiProb(int brlen) {

    key = std::make_tuple(brlen,0,0);
    ti_map::iterator it = tiMap.find(key);
    if (it != tiMap.end())
        return it->second;
    return nullptr;
}

TransitionProbabilities* TransitionProbabilitiesMngr::getTiProb(Node* node) {
    
    if (modelType == custom_f81 && isUgandaRateVariable == true) {
        // For custom_f81 with variable Uganda rates, compute dwell times
        double t0=0.0, t1=0.0, t2=0.0;
        if (node->getAncestor() != nullptr) {
            modelPtr->computeDwellTimes(node->getTime(), node->getAncestor()->getTime(), &t0, &t1, &t2);
        }
        key = std::make_tuple((int)t0, (int)t1, (int)t2);
        ti_map::iterator it = tiMap.find(key);
        if (it != tiMap.end())
            return it->second;
        return nullptr;
    } else {
        // For other models, fall back to the simple branch length lookup
        return getTiProb(node->getBrlen());
    }
}

bool TransitionProbabilitiesMngr::isGPUAvailable(void) const {

    return gpuBatcher != nullptr && gpuBatcher->isAvailable();
}

const char* TransitionProbabilitiesMngr::getGPUDeviceName(void) const {

    if (gpuBatcher != nullptr)
        return gpuBatcher->getDeviceName();
    return "N/A";
}

void TransitionProbabilitiesMngr::printMap(void) {

    for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
        {
        std::cout << std::get<0>(it->first) << ":" << std::get<1>(it->first) << " (" << it->second << ")" << std::endl;
        }
}

void TransitionProbabilitiesMngr::updateTransitionProbabilities(void) {

    RateMatrix* Q = modelPtr->getRateMatrix();

    // For JC69/F81/custom_F81 we use the fast analytic forms (threaded tasks).
    // The batched/GPU matrix-exponential path is only relevant for the general (GTR) model.
    if (modelType != gtr)
        {
        updateThreaded(Q);
        return;
        }

    size_t numMatrices = tiMap.size();

    // select compute path based on backend setting
    bool useBatched = false;

    // analytic models never use the matrix-exponential backends
    if (modelType != gtr)
        useBatched = false;

    switch (computeBackend)
        {
        case TiProbComputeBackend::Auto:
            useBatched = (numMatrices >= batchThreshold);
            break;
        case TiProbComputeBackend::ThreadedTasks:
            useBatched = false;
            break;
        case TiProbComputeBackend::BatchedAccelerate:
            useBatched = true;
            break;
        }
        
    if (useBatched)
        updateBatched(Q);
    else
        updateThreaded(Q);
    
#   if 0
    checkTiProbs();
#   endif
}

/* Original threaded approach: one task per branch length
   Each task creates exp(Q * brlen) independently */
void TransitionProbabilitiesMngr::updateThreaded(DoubleMatrix* Q) {

    if (modelType == custom_f81 && isUgandaRateVariable == true)
        {
        // variable rates of Uganda transmission
        TransitionProbabilitiesTask* task = tasks;

        // For analytic models we need the current stationary frequencies.
        // (For GTR this pointer is ignored by TransitionProbabilitiesTask.)
        const std::vector<double>& pi = modelPtr->getPi();
        const std::vector<double>& exch = modelPtr->getR();
        double r = modelPtr->getSubstitutionRate();
        double kappaLockdown = modelPtr->getKappaLockdown();
        double kappaNoLockdown = modelPtr->getKappaNoLockdown();

        double p = pi[ugandaIdx];
        double q = 1.0 - p;
        double s = 0.0;
        for (size_t i = 0, numStates=pi.size(); i < numStates; ++i)
            s += pi[i] * pi[i];
        double r1 = 1.0 - s + 2.0 * (kappaNoLockdown - 1.0) * p * q;
        double lambda1 = (-1.0 + (1.0 - kappaNoLockdown) * p) / r1;
        double lambda2 = -kappaNoLockdown / r1;
        double r2 = 1.0 - s + 2.0 * (kappaLockdown - 1.0) * p * q;
        double lambda12 = (-1.0 + (1.0 - kappaLockdown) * p) / r2;
        double lambda22 = -kappaLockdown / r2;
            
        int id = 0;
        for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
            {
            int sum = std::get<0>(it->first) + std::get<1>(it->first) + std::get<2>(it->first);
            double v1 = (double)(std::get<0>(it->first)) * r;
            double v2 = (double)(std::get<1>(it->first)) * r;
            double v3 = (double)(std::get<2>(it->first)) * r;
            if (sum == 0)
                {
                v1 = MIN_BRLEN/3.0;
                v2 = MIN_BRLEN/3.0;
                v3 = MIN_BRLEN/3.0;
                }
            it->second->setCalculatedBrlen(v1 + v2 + v3);
            task->init(id, (int)dim, v1, v2, v3, (DoubleMatrix*)Q, it->second, &pi, &exch, kappaNoLockdown, kappaLockdown, ugandaIdx);
            task->setCustomModelParms(p, q, lambda1, lambda2, lambda12, lambda22);
                
            threadPool->pushTask(task);
            ++task;
            ++id;
            }

        threadPool->wait();
        }
    else 
        {
        // time-homogeneous models
        TransitionProbabilitiesTask* task = tasks;

        // For analytic models we need the current stationary frequencies.
        // (For GTR this pointer is ignored by TransitionProbabilitiesTask.)
        const std::vector<double>& pi = modelPtr->getPi();
        const std::vector<double>& exch = modelPtr->getR();
        double kappa = modelPtr->getKappa();
        double r = modelPtr->getSubstitutionRate();
            
        int id = 0;
        for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
            {
            double v = (double)(std::get<0>(it->first)) * r;
            if (std::get<0>(it->first) == 0)
                v = MIN_BRLEN;
            it->second->setCalculatedBrlen(v);
            task->init(id, (int)dim, v, 0.0, 0.0, (DoubleMatrix*)Q, it->second, &pi, &exch, kappa, kappa, ugandaIdx);

            threadPool->pushTask(task);
            ++task;
            ++id;
            }

        threadPool->wait();
        }
}

/* Batched approach: all branch lengths processed together
   Uses GPUMatrixExponentialBatch which:
     - Shares the same Q matrix across all computations
     - Uses static buffer pools to avoid heap fragmentation
     - Leverages ThreadPool for parallel Padé[13/13] computation
     - Each worker thread has its own working buffers

   This is more efficient than the threaded approach when:
     - There are many unique branch lengths (>= batchThreshold)
     - The rate matrix Q is large (>= 64x64)*/

void TransitionProbabilitiesMngr::updateBatched(DoubleMatrix* Q) {

    if (modelType != gtr)
        Msg::error("Batched transition-probability computation is only valid for gtr");

    // clear and populate the batch vectors
    batchBranchLengths.clear();
    batchOutputs.clear();
    
    double r = modelPtr->getSubstitutionRate();
    
    for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
        {
        double v = (double)(std::get<0>(it->first)) * r;
        if (std::get<0>(it->first) == 0)
            v = MIN_BRLEN;
        it->second->setCalculatedBrlen(v);
        
        batchBranchLengths.push_back(v);
        batchOutputs.push_back(it->second);  // TransitionProbabilities inherits from DoubleMatrix
        }
    
    // Compute all matrix exponentials in parallel
    // GPUMatrixExponentialBatch::computeBatch uses the ThreadPool internally
    gpuBatcher->computeBatch(*Q, batchBranchLengths, batchOutputs, threadPool);
}
