#ifndef GPUMatrixExponentialBatch_hpp
#define GPUMatrixExponentialBatch_hpp

#include "Container.hpp"
#include "Threads.hpp"
#include <vector>



class GPUMatrixExponentialBatch {

    public:
        static GPUMatrixExponentialBatch&   getInstance(void);
        
        bool                isAvailable(void) const;
        const char*         getDeviceName(void) const;
        
                            // compute exp(Q * t_i) for each branch length t_i
                            // Uses ThreadPool for parallel computation
        void                computeBatch(const DoubleMatrix& Q,
                                         const std::vector<double>& branchLengths,
                                         std::vector<DoubleMatrix*>& outputs,
                                         ThreadPool* pool);
        
                            // Legacy interface without ThreadPool (single-threaded)
        void                computeBatch(const DoubleMatrix& Q,
                                         const std::vector<double>& branchLengths,
                                         std::vector<DoubleMatrix*>& outputs);
        
                            // thresholds for when to use GPU vs CPU 
        void                setMinBatchSize(size_t size) { minBatchSize = size; }
        size_t              getMinBatchSize(void) const { return minBatchSize; }
        void                setMinMatrixSize(size_t size) { minMatrixSize = size; }
        size_t              getMinMatrixSize(void) const { return minMatrixSize; }
        
                            // returns true if GPU should be used for this workload 
        bool                shouldUseGPU(size_t batchSize, size_t matrixSize) const;
        
    private:
                            GPUMatrixExponentialBatch(void);
                           ~GPUMatrixExponentialBatch(void);
                            GPUMatrixExponentialBatch(const GPUMatrixExponentialBatch&) = delete;
        GPUMatrixExponentialBatch& operator=(const GPUMatrixExponentialBatch&) = delete;
        
        void                initialize(void);
        void                computeSingleThreaded(const double* Q, size_t n,
                                                  const std::vector<double>& branchLengths,
                                                  std::vector<DoubleMatrix*>& outputs);
        
                            // Metal objects (for GPU compute) 
        void*               device;
        void*               commandQueue;
        
                            // unused - reserved for Metal compute pipelines 
        void*               scaleMatrixPipeline;
        void*               matrixMultiplyPipeline;
        void*               matrixAddPipeline;
        void*               absoluteValuePipeline;
        
        void*               qBuffer;
        void*               branchLengthBuffer;
        void*               workBuffer;
        size_t              qBufferCapacity;
        size_t              branchLengthBufferCapacity;
        size_t              workBufferCapacity;
        
        bool                gpuAvailable;
        char                deviceName[256];
        size_t              minBatchSize;
        size_t              minMatrixSize;
};

#endif
