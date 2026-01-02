// -------------------------------------------------------------------
// GPUMatrixExponentialBatch.cpp
// -------------------------------------------------------------------
// Portable implementation of batch matrix exponential computation.
// Uses Padé[13/13] with scaling and squaring, parallelized via ThreadPool.
//
// On macOS: Uses Accelerate BLAS/LAPACK for optimal performance
// On Linux: Uses OpenBLAS if available, otherwise portable fallback
//
// The "GPU" in the name is aspirational - currently uses optimized CPU
// computation, but the interface is designed to support Metal compute
// shaders in the future.
// -------------------------------------------------------------------

#include "GPUMatrixExponentialBatch.hpp"
#include "MathCacheAccelerated.hpp"
#include <atomic>
#include <cmath>
#include <cstring>
#include <iostream>
#include <thread>
#include <vector>



// -------------------------------------------------------------------
// Static buffer pool for working memory
// -------------------------------------------------------------------
// Each worker thread acquires a MathCacheAccelerated from the pool,
// uses it for matrix exponential computation, then releases it.
// This eliminates heap fragmentation from repeated alloc/free cycles.
// -------------------------------------------------------------------

class CachePool {

    static constexpr size_t MAX_CACHES = 64;
    MathCacheAccelerated* caches[MAX_CACHES];
    std::atomic<bool> inUse[MAX_CACHES];
    
    public:
        CachePool(void) {
            for (size_t i = 0; i < MAX_CACHES; i++) {
                caches[i] = nullptr;
                inUse[i].store(false, std::memory_order_relaxed);
            }
        }
        
        ~CachePool(void) {
            for (size_t i = 0; i < MAX_CACHES; i++) {
                delete caches[i];
            }
        }
        
        MathCacheAccelerated* acquire(void) {
            // Try to find a free cache
            for (size_t i = 0; i < MAX_CACHES; i++) {
                bool expected = false;
                if (inUse[i].compare_exchange_strong(expected, true,
                        std::memory_order_acquire, std::memory_order_relaxed)) {
                    // Lazy initialization
                    if (caches[i] == nullptr) {
                        caches[i] = new MathCacheAccelerated();
                    }
                    return caches[i];
                }
            }
            // All caches in use - spin on first one
            while (true) {
                bool expected = false;
                if (inUse[0].compare_exchange_strong(expected, true,
                        std::memory_order_acquire, std::memory_order_relaxed)) {
                    if (caches[0] == nullptr) {
                        caches[0] = new MathCacheAccelerated();
                    }
                    return caches[0];
                }
                std::this_thread::yield();
            }
        }
        
        void release(MathCacheAccelerated* cache) {
            for (size_t i = 0; i < MAX_CACHES; i++) {
                if (caches[i] == cache) {
                    inUse[i].store(false, std::memory_order_release);
                    return;
                }
            }
        }
};

static CachePool& getCachePool(void) {
    static CachePool pool;
    return pool;
}

// -------------------------------------------------------------------
// Task for computing a single matrix exponential
// -------------------------------------------------------------------

class MatrixExpTask : public ThreadTask {

    public:
        MatrixExpTask(void) : Q(nullptr), n(0), t(0), P(nullptr) {}
        
        void set(const DoubleMatrix* qMatrix, double branchLength, DoubleMatrix* output) {
            Q = qMatrix;
            n = qMatrix->getNumRows();
            t = branchLength;
            P = output;
        }
        
        void run(void) override {
            CachePool& pool = getCachePool();
            MathCacheAccelerated* cache = pool.acquire();
            cache->computeMatrixExponential(*Q, t, *P);
            pool.release(cache);
        }
    
    private:
        const DoubleMatrix*   Q;
        size_t              n;
        double              t;
        DoubleMatrix*         P;
};

// -------------------------------------------------------------------
// Pool of reusable task objects
// -------------------------------------------------------------------

class TaskPool {

    static constexpr size_t MAX_TASKS = 2048;
    MatrixExpTask tasks[MAX_TASKS];
    std::atomic<size_t> nextIndex{0};
    
    public:
        void reset(void) {
            nextIndex.store(0, std::memory_order_relaxed);
        }
        
        MatrixExpTask* getTask(void) {
            size_t idx = nextIndex.fetch_add(1, std::memory_order_relaxed);
            if (idx >= MAX_TASKS) {
                // Overflow - shouldn't happen with reasonable batch sizes
                return &tasks[0];
            }
            return &tasks[idx];
        }
};

static TaskPool& getTaskPool(void) {
    static TaskPool pool;
    return pool;
}

// -------------------------------------------------------------------
// GPUMatrixExponentialBatch implementation
// -------------------------------------------------------------------

GPUMatrixExponentialBatch& GPUMatrixExponentialBatch::getInstance(void) {

    static GPUMatrixExponentialBatch instance;
    return instance;
}

GPUMatrixExponentialBatch::GPUMatrixExponentialBatch(void) {

    device = nullptr;
    commandQueue = nullptr;
    scaleMatrixPipeline = nullptr;
    matrixMultiplyPipeline = nullptr;
    matrixAddPipeline = nullptr;
    absoluteValuePipeline = nullptr;
    
    qBuffer = nullptr;
    branchLengthBuffer = nullptr;
    workBuffer = nullptr;
    qBufferCapacity = 0;
    branchLengthBufferCapacity = 0;
    workBufferCapacity = 0;
    
    gpuAvailable = false;
    deviceName[0] = '\0';
    minBatchSize = 8;
    minMatrixSize = 32;
    
    initialize();
}

GPUMatrixExponentialBatch::~GPUMatrixExponentialBatch(void) {

    // No Metal resources to clean up in portable version
}

void GPUMatrixExponentialBatch::initialize(void) {

    // In the portable version, we use optimized CPU computation
    // via BLAS/LAPACK (Accelerate on macOS, OpenBLAS on Linux)
    
#if defined(__APPLE__)
    gpuAvailable = true;  // Accelerate is always available on macOS
    strncpy(deviceName, "Apple Accelerate (CPU)", 255);
#elif defined(USE_OPENBLAS) || defined(USE_CBLAS)
    gpuAvailable = true;
    strncpy(deviceName, "OpenBLAS (CPU)", 255);
#else
    gpuAvailable = true;
    strncpy(deviceName, "Portable (CPU)", 255);
#endif
    
    deviceName[255] = '\0';
    
    std::cout << "   * Matrix exponential batch: " << deviceName << std::endl;
}

bool GPUMatrixExponentialBatch::isAvailable(void) const {

    return gpuAvailable;
}

const char* GPUMatrixExponentialBatch::getDeviceName(void) const {

    return deviceName;
}

bool GPUMatrixExponentialBatch::shouldUseGPU(size_t batchSize, size_t matrixSize) const {

    if (!gpuAvailable)
        return false;
    
    return batchSize >= minBatchSize && matrixSize >= minMatrixSize;
}

void GPUMatrixExponentialBatch::computeBatch(const DoubleMatrix& Q,
                                              const std::vector<double>& branchLengths,
                                              std::vector<DoubleMatrix*>& outputs,
                                              ThreadPool* pool) {

    size_t count = branchLengths.size();
    if (count == 0)
        return;
    
    // Single-threaded fallback if no pool or small batch
    if (pool == nullptr || count < 4) {
        computeSingleThreaded(Q.begin(), Q.getNumRows(), branchLengths, outputs);
        return;
    }
    
    // Reset task pool for this batch
    TaskPool& taskPool = getTaskPool();
    taskPool.reset();
    
    // Create and submit tasks
    for (size_t i = 0; i < count; i++) {
        MatrixExpTask* task = taskPool.getTask();
        task->set(&Q, branchLengths[i], outputs[i]);
        pool->pushTask(task);
    }
    
    // Wait for all tasks to complete
    pool->wait();
}

void GPUMatrixExponentialBatch::computeBatch(const DoubleMatrix& Q,
                                              const std::vector<double>& branchLengths,
                                              std::vector<DoubleMatrix*>& outputs) {

    // Legacy interface - single-threaded
    computeSingleThreaded(Q.begin(), Q.getNumRows(), branchLengths, outputs);
}

void GPUMatrixExponentialBatch::computeSingleThreaded(const double* qData, size_t n,
                                                       const std::vector<double>& branchLengths,
                                                       std::vector<DoubleMatrix*>& outputs) {

    CachePool& pool = getCachePool();
    MathCacheAccelerated* cache = pool.acquire();
    
    // Create a temporary matrix wrapper for Q
    DoubleMatrix Q(n, n);
    memcpy(Q.begin(), qData, n * n * sizeof(double));
    
    for (size_t i = 0; i < branchLengths.size(); i++) {
        cache->computeMatrixExponential(Q, branchLengths[i], *outputs[i]);
    }
    
    pool.release(cache);
}
