#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#ifndef ACCELERATE_NEW_LAPACK
#define ACCELERATE_NEW_LAPACK
#endif
#import <Accelerate/Accelerate.h>

#include "GPUMatrixExponentialBatch.hpp"
#include "MathCache.hpp"
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>
#include <atomic>



static const double PADE_B[14] = {
    64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
    1187353796428800.0, 129060195264000.0, 10559470521600.0,
    670442572800.0, 33522128640.0, 1323241920.0, 40840800.0,
    960960.0, 16380.0, 182.0, 1.0
};

static const double THETA_13 = 5.371920351148152;



/* Static buffer pool that persists for the program's lifetime.
   Each worker thread acquires a buffer set, uses it for all its
   matrix exponential computations, then releases it back to the pool.
   This eliminates heap fragmentation from repeated alloc/free cycles. */
struct WorkingBuffers {

    std::vector<double> A, A2, A4, A6;
    std::vector<double> U, V, temp1, temp2;
    std::vector<double> temp1T, temp2T;
    std::vector<__LAPACK_int> ipiv;
    size_t lastSize = 0;
    std::atomic<bool> inUse{false};
    
    void resize(size_t nn, size_t n) {
        if (lastSize != nn) 
            {
            A.resize(nn);
            A2.resize(nn);
            A4.resize(nn);
            A6.resize(nn);
            U.resize(nn);
            V.resize(nn);
            temp1.resize(nn);
            temp2.resize(nn);
            temp1T.resize(nn);
            temp2T.resize(nn);
            ipiv.resize(n);
            lastSize = nn;
            }
    }
};

class BufferPool {

    static constexpr size_t MAX_BUFFERS = 64;
    WorkingBuffers buffers[MAX_BUFFERS];
    
    public:
        WorkingBuffers* acquire(void) {
            // try to find a free buffer
            for (size_t i = 0; i < MAX_BUFFERS; i++) 
                {
                bool expected = false;
                if (buffers[i].inUse.compare_exchange_strong(expected, true, std::memory_order_acquire, std::memory_order_relaxed)) 
                    {
                    return &buffers[i];
                    }
                }
            // all buffers in use - spin on first buffer
            while (true) 
                {
                bool expected = false;
                if (buffers[0].inUse.compare_exchange_strong(expected, true, std::memory_order_acquire, std::memory_order_relaxed)) 
                    {
                    return &buffers[0];
                    }
                std::this_thread::yield();
                }
        }
        
        void release(WorkingBuffers* buf) {
            buf->inUse.store(false, std::memory_order_release);
        }
};

static BufferPool& getBufferPool(void) {
    static BufferPool pool;
    return pool;
}

static void computeSingleExponentialAccelerate(const double* Q, size_t n, double t, double* P, WorkingBuffers* buf) {

    size_t nn = n * n;
    int N = static_cast<int>(n);
    
    // resize buffers if needed (typically only on first use) 
    buf->resize(nn, n);
    
    // use references for cleaner code
    std::vector<double>& A = buf->A;
    std::vector<double>& A2 = buf->A2;
    std::vector<double>& A4 = buf->A4;
    std::vector<double>& A6 = buf->A6;
    std::vector<double>& U = buf->U;
    std::vector<double>& V = buf->V;
    std::vector<double>& temp1 = buf->temp1;
    std::vector<double>& temp2 = buf->temp2;
    std::vector<double>& temp1T = buf->temp1T;
    std::vector<double>& temp2T = buf->temp2T;
    std::vector<__LAPACK_int>& ipiv = buf->ipiv;
    
    // A = Q * t 
    for (size_t i = 0; i < nn; i++)
        A[i] = Q[i] * t;
    
    // compute infinity norm
    double normA = 0.0;
    for (size_t i = 0; i < n; i++)
        {
        double rowSum = 0.0;
        for (size_t j = 0; j < n; j++)
            rowSum += std::fabs(A[i * n + j]);
        if (rowSum > normA)
            normA = rowSum;
        }
    
    // compute scaling factor
    int s = 0;
    if (normA > THETA_13)
        {
        s = static_cast<int>(std::ceil(std::log2(normA / THETA_13)));
        double scale = 1.0 / (1L << s);
        cblas_dscal(static_cast<int>(nn), scale, A.data(), 1);
        }
    
    // compute A^2, A^4, A^6 using BLAS
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, A.data(), N, A.data(), N, 0.0, A2.data(), N);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, A2.data(), N, A2.data(), N, 0.0, A4.data(), N);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, A2.data(), N, A4.data(), N, 0.0, A6.data(), N);
    
    // build U polynomial
    std::fill(temp1.begin(), temp1.end(), 0.0);
    for (size_t i = 0; i < n; i++)
        temp1[i * n + i] = PADE_B[7];
    for (size_t i = 0; i < nn; i++)
        temp1[i] += PADE_B[9] * A2[i] + PADE_B[11] * A4[i] + PADE_B[13] * A6[i];
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, A6.data(), N, temp1.data(), N, 0.0, temp2.data(), N);
    
    for (size_t i = 0; i < nn; i++)
        temp2[i] += PADE_B[5] * A4[i] + PADE_B[3] * A2[i];
    for (size_t i = 0; i < n; i++)
        temp2[i * n + i] += PADE_B[1];
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, A.data(), N, temp2.data(), N, 0.0, U.data(), N);
    
    // build V polynomial 
    std::fill(temp1.begin(), temp1.end(), 0.0);
    for (size_t i = 0; i < n; i++)
        temp1[i * n + i] = PADE_B[6];
    for (size_t i = 0; i < nn; i++)
        temp1[i] += PADE_B[8] * A2[i] + PADE_B[10] * A4[i] + PADE_B[12] * A6[i];
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, A6.data(), N, temp1.data(), N, 0.0, V.data(), N);
    
    for (size_t i = 0; i < nn; i++)
        V[i] += PADE_B[4] * A4[i] + PADE_B[2] * A2[i];
    for (size_t i = 0; i < n; i++)
        V[i * n + i] += PADE_B[0];
    
    // solve (V - U) * P = (V + U) using LAPACK
    for (size_t i = 0; i < nn; i++)
        {
        temp1[i] = V[i] - U[i];
        temp2[i] = V[i] + U[i];
        }
    
    __LAPACK_int info;
    __LAPACK_int LN = static_cast<__LAPACK_int>(n);
    
    // transpose for column-major LAPACK
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            {
            temp1T[j * n + i] = temp1[i * n + j];
            temp2T[j * n + i] = temp2[i * n + j];
            }
    
    dgetrf_(&LN, &LN, temp1T.data(), &LN, ipiv.data(), &info);
    if (info == 0)
        dgetrs_("N", &LN, &LN, temp1T.data(), &LN, ipiv.data(), temp2T.data(), &LN, &info);
    
    // transpose back to row-major
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            P[i * n + j] = temp2T[j * n + i];
    
    // squaring phase
    for (int sq = 0; sq < s; sq++)
        {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    N, N, N, 1.0, P, N, P, N, 0.0, temp1.data(), N);
        memcpy(P, temp1.data(), nn * sizeof(double));
        }
    
    // absolute value cleanup
    for (size_t i = 0; i < nn; i++)
        P[i] = std::fabs(P[i]);
}

/* Task class for computing a single matrix exponential.
   Acquires a working buffer from the pool, computes exp(Q*t),
   and releases the buffer when done. */
class MatrixExponentialTask : public ThreadTask {

    public:
                            MatrixExponentialTask(void) : Q(nullptr), n(0), t(0), P(nullptr) {}
        
        void                set(const double* qMatrix, size_t dim, double branchLength, double* output) {
                                Q = qMatrix;
                                n = dim;
                                t = branchLength;
                                P = output;
                            }
        
        void                run(void) override {
                                BufferPool& pool = getBufferPool();
                                WorkingBuffers* buf = pool.acquire();
                                computeSingleExponentialAccelerate(Q, n, t, P, buf);
                                pool.release(buf);
                            }
    
    private:
        const double*       Q;
        size_t              n;
        double              t;
        double*             P;
};

/* Pool of reusable MatrixExponentialTask objects.
   Avoids allocating new task objects for every matrix exponential. */
class TaskPool {

    static constexpr size_t MAX_TASKS = 1024;
    MatrixExponentialTask tasks[MAX_TASKS];
    std::atomic<size_t> nextIndex{0};
    
    public:
        void reset(void) {
            nextIndex.store(0, std::memory_order_relaxed);
        }
        
        MatrixExponentialTask* getTask() {
            size_t idx = nextIndex.fetch_add(1, std::memory_order_relaxed);
            if (idx >= MAX_TASKS) {
                // overflow - this shouldn't happen with reasonable batch sizes
                return &tasks[0];
            }
            return &tasks[idx];
        }
};

static TaskPool& getTaskPool(void) {
    static TaskPool pool;
    return pool;
}

GPUMatrixExponentialBatch& GPUMatrixExponentialBatch::getInstance(void) {

    static GPUMatrixExponentialBatch instance;
    return instance;
}

GPUMatrixExponentialBatch::GPUMatrixExponentialBatch(void) {

    device = nil;
    commandQueue = nil;
    scaleMatrixPipeline = nil;
    matrixMultiplyPipeline = nil;
    matrixAddPipeline = nil;
    absoluteValuePipeline = nil;
    
    qBuffer = nil;
    branchLengthBuffer = nil;
    workBuffer = nil;
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

    @autoreleasepool
        {
        if (qBuffer) CFRelease(qBuffer);
        if (branchLengthBuffer) CFRelease(branchLengthBuffer);
        if (workBuffer) CFRelease(workBuffer);
        if (scaleMatrixPipeline) CFRelease(scaleMatrixPipeline);
        if (matrixMultiplyPipeline) CFRelease(matrixMultiplyPipeline);
        if (absoluteValuePipeline) CFRelease(absoluteValuePipeline);
        if (commandQueue) CFRelease(commandQueue);
        if (device) CFRelease(device);
        }
}

void GPUMatrixExponentialBatch::initialize(void) {

    @autoreleasepool
        {
        // get default GPU device - works on all Apple Silicon
        id<MTLDevice> mtlDevice = MTLCreateSystemDefaultDevice();
        if (!mtlDevice)
            {
            gpuAvailable = false;
            strcpy(deviceName, "No GPU available");
            return;
            }
        
        device = (__bridge_retained void*)mtlDevice;
        strncpy(deviceName, [[mtlDevice name] UTF8String], 255);
        deviceName[255] = '\0';
        
        id<MTLCommandQueue> queue = [mtlDevice newCommandQueue];
        if (!queue)
            {
            CFRelease(device);
            device = nil;
            gpuAvailable = false;
            return;
            }
        commandQueue = (__bridge_retained void*)queue;
        
        // GPU is available for future use
        // currently using optimized CPU path which is excellent on Apple Silicon 
        gpuAvailable = true;
        
        std::cout << "   * GPU batch: Initialized with device " <<  deviceName << std::endl;
        std::cout << "   * GPU batch: Using optimized Accelerate BLAS with ThreadPool" << std::endl;
        }
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

void GPUMatrixExponentialBatch::computeBatch(const DoubleMatrix& Q, const std::vector<double>& branchLengths, std::vector<DoubleMatrix*>& outputs, ThreadPool* pool) {

    size_t count = branchLengths.size();
    if (count == 0)
        return;
    
    size_t n = Q.getNumRows();
    const double* qData = Q.begin();
    
    /* Use the existing ThreadPool for parallel computation.
       This avoids creating/destroying threads on every call,
       which was causing heap fragmentation on Intel Macs.
       
       Each task acquires a working buffer from the static pool,
       computes one matrix exponential, and releases the buffer. */
    
    if (pool == nullptr || count < 4)
        {
        // single-threaded fallback
        computeSingleThreaded(qData, n, branchLengths, outputs);
        return;
        }
    
    // reset task pool for this batch
    TaskPool& taskPool = getTaskPool();
    taskPool.reset();
    
    // create and submit tasks
    for (size_t i = 0; i < count; i++)
        {
        MatrixExponentialTask* task = taskPool.getTask();
        task->set(qData, n, branchLengths[i], outputs[i]->begin());
        pool->pushTask(task);
        }
    
    // wait for all tasks to complete
    pool->wait();
}

void GPUMatrixExponentialBatch::computeBatch(const DoubleMatrix& Q, const std::vector<double>& branchLengths, std::vector<DoubleMatrix*>& outputs) {

    // legacy interface - single-threaded
    size_t n = Q.getNumRows();
    computeSingleThreaded(Q.begin(), n, branchLengths, outputs);
}

void GPUMatrixExponentialBatch::computeSingleThreaded(const double* Q, size_t n, const std::vector<double>& branchLengths, std::vector<DoubleMatrix*>& outputs) {

    BufferPool& pool = getBufferPool();
    WorkingBuffers* buf = pool.acquire();
    
    for (size_t i = 0; i < branchLengths.size(); i++)
        computeSingleExponentialAccelerate(Q, n, branchLengths[i], outputs[i]->begin(), buf);
    
    pool.release(buf);
}
