#ifndef TransitionProbabilitiesMngr_hpp
#define TransitionProbabilitiesMngr_hpp

#include <iostream>
#include <map>
#include <vector>
#include "Container.hpp"
#include "TransitionProbabilities.hpp"
class GPUMatrixExponentialBatch;
class Model;
class ThreadPool;
class TransitionProbabilities;
class Tree;

struct TransitionProbabilityPair {

    TransitionProbabilities*    tipr[2];
};

typedef std::map<std::pair<int,int>,TransitionProbabilities*> ti_map;


// -------------------------------------------------------------------
// Compute backend selection for transition probability calculation
// -------------------------------------------------------------------
// Auto: Chooses batched for large batch sizes (>= batchThreshold)
// ThreadedTasks: Original approach - one task per branch length
// BatchedAccelerate: Uses GPUMatrixExponentialBatch with ThreadPool
//                    (same Q matrix, multiple branch lengths in parallel)
// -------------------------------------------------------------------

enum class TiProbComputeBackend {
    Auto,              // Automatically choose based on batch size
    ThreadedTasks,     // Original: individual tasks per branch length
    BatchedAccelerate  // Batched: uses GPUMatrixExponentialBatch
};



class TransitionProbabilitiesMngr {

    public:
                                        TransitionProbabilitiesMngr(void) = delete;
                                        TransitionProbabilitiesMngr(const TransitionProbabilitiesMngr& m) = delete;
                                        TransitionProbabilitiesMngr(Model* m, Tree* t, size_t d, ThreadPool* tp, size_t ugi);
                                       ~TransitionProbabilitiesMngr(void);
        TransitionProbabilities*        getTiProb(int brlen);
        void                            printMap(void);
        void                            updateTransitionProbabilities(void);
        
        // Backend selection
        void                            setComputeBackend(TiProbComputeBackend backend) { computeBackend = backend; }
        TiProbComputeBackend            getComputeBackend(void) const { return computeBackend; }
        void                            setBatchThreshold(size_t threshold) { batchThreshold = threshold; }
        size_t                          getBatchThreshold(void) const { return batchThreshold; }
        
        // Diagnostics
        bool                            isGPUAvailable(void) const;
        const char*                     getGPUDeviceName(void) const;
        size_t                          getNumUniqueBranchLengths(void) const { return tiMap.size(); }
                
    private:
        void                            checkTiProbs(void);
        void                            updateThreaded(DoubleMatrix* Q);
        void                            updateBatched(DoubleMatrix* Q);
        
        Model*                          modelPtr;
        ThreadPool*                     threadPool;
        GPUMatrixExponentialBatch*      gpuBatcher;
        TransitionProbabilitiesTask*    tasks;
        size_t                          dim;
        ti_map                          tiMap;
        std::pair<int,int>              key;
        SubModel                        modelType;
        size_t                          ugandaIdx;
        
        // Batch computation settings
        TiProbComputeBackend            computeBackend;
        size_t                          batchThreshold;
        
        // Pre-allocated vectors for batch operations (avoid repeated allocation)
        std::vector<double>             batchBranchLengths;
        std::vector<DoubleMatrix*>      batchOutputs;
};

#endif
