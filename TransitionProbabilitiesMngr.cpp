#include "GPUMatrixExponentialBatch.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "RateMatrix.hpp"
#include "Threads.hpp"
#include "TransitionProbabilities.hpp"
#include "TransitionProbabilitiesMngr.hpp"
#include "Tree.hpp"

#define MIN_BRLEN 10e-5



TransitionProbabilitiesMngr::TransitionProbabilitiesMngr(Model* m, Tree* t, size_t d, ThreadPool* tp) {

    // pointer to the model object
    modelPtr = m;
    
    threadPool = tp;
    
    // set the dimensions of the rate matrices
    dim = d;
    
    // Initialize GPU batcher singleton
    gpuBatcher = &GPUMatrixExponentialBatch::getInstance();
    
    // Default settings
    computeBackend = TiProbComputeBackend::Auto;
    batchThreshold = 8;  // Use batched approach for >= 8 matrices
    
    // For large matrices, the batched approach is always better
    if (dim >= 64)
        batchThreshold = 4;
    
    // allocate matrices for the unique branch lengths
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (int i=0, n=(int)dpSeq.size(); i<n; i++)
        {
        Node* p = dpSeq[i];
        int v = p->getBrlen();
        key.first = v;
        key.second = 0;
        ti_map::iterator it = tiMap.find(key);
        if (it == tiMap.end())
            {
            TransitionProbabilities* ti = new TransitionProbabilities(dim);
            ti->setBrlen(v);
            tiMap.insert( std::make_pair(key,ti) );
            }
        }
    
    // Pre-allocate vectors for batch operations
    batchBranchLengths.reserve(tiMap.size());
    batchOutputs.reserve(tiMap.size());
    
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
}

void TransitionProbabilitiesMngr::checkTiProbs(void) {

    for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
        {
        TransitionProbabilities* p = it->second;
        size_t n = p->dim();
        std::cout << "matrix " << it->first.first << " " << it->first.second << std::endl;
        for (size_t i=0; i<n; i++)
            {
            double sum = 0.0;
            for (size_t j=0; j<n; j++)
                {
                double x = (*p)(i,j);
                sum += x;
                if (x < 0.0)
                    Msg::warning("Negative transition probability!!!!!");
                }
            std::cout << sum << std::endl;
            }
        }
}

TransitionProbabilities* TransitionProbabilitiesMngr::getTiProb(int brlen) {

    key.first = brlen;
    key.second = 0;
    ti_map::iterator it = tiMap.find(key);
    if (it != tiMap.end())
        return it->second;
    return nullptr;
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
        std::cout << it->first.first << ":" << it->first.second << " (" << it->second << ")" << std::endl;
        }
}

void TransitionProbabilitiesMngr::updateTransitionProbabilities(void) {

    RateMatrix* Q = modelPtr->getRateMatrix();
    size_t numMatrices = tiMap.size();
    
    // Select compute path based on backend setting
    bool useBatched = false;
    
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

// -------------------------------------------------------------------
// Original threaded approach: one task per branch length
// Each task creates exp(Q * brlen) independently
// -------------------------------------------------------------------

void TransitionProbabilitiesMngr::updateThreaded(DoubleMatrix* Q) {

    auto tasks = new TransitionProbabilitiesTask[tiMap.size()];
    auto task = tasks;
    
    int id = 0;
    double r = modelPtr->getSubstitutionRate();
    for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
        {
        double v = (double)it->first.first * r;
        if (it->first.first == 0)
            v = MIN_BRLEN;
        it->second->setCalculatedBrlen(v);
        task->init(id, (int)dim, v, (DoubleMatrix*)Q, it->second);
        threadPool->pushTask(task);
        ++task;
        ++id;
        }
        
    threadPool->wait();
    delete [] tasks;
}

// -------------------------------------------------------------------
// Batched approach: all branch lengths processed together
// Uses GPUMatrixExponentialBatch which:
//   - Shares the same Q matrix across all computations
//   - Uses static buffer pools to avoid heap fragmentation
//   - Leverages ThreadPool for parallel Padé[13/13] computation
//   - Each worker thread has its own working buffers
//
// This is more efficient than the threaded approach when:
//   - There are many unique branch lengths (>= batchThreshold)
//   - The rate matrix Q is large (>= 64x64)
// -------------------------------------------------------------------

void TransitionProbabilitiesMngr::updateBatched(DoubleMatrix* Q) {

    // Clear and populate the batch vectors
    batchBranchLengths.clear();
    batchOutputs.clear();
    
    double r = modelPtr->getSubstitutionRate();
    
    for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
        {
        double v = (double)it->first.first * r;
        if (it->first.first == 0)
            v = MIN_BRLEN;
        it->second->setCalculatedBrlen(v);
        
        batchBranchLengths.push_back(v);
        batchOutputs.push_back(it->second);  // TransitionProbabilities inherits from DoubleMatrix
        }
    
    // Compute all matrix exponentials in parallel
    // GPUMatrixExponentialBatch::computeBatch uses the ThreadPool internally
    gpuBatcher->computeBatch(*Q, batchBranchLengths, batchOutputs, threadPool);
}
