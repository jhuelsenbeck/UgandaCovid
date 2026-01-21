#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <mutex>
#include "CondLikeJob.hpp"
#include "CondLikeJobMngr.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "TransitionProbabilities.hpp"



CondLikeJob::CondLikeJob(ThreadPool* tp, int na) : ThreadTask(), threadPool(tp), numStates(na) {

    jobId = 0;
    numDependencies = 0;
    lnScaler = 0.0;
    numResolvedDependencies = 0;
    
    // Allocate working buffers with alignment for SIMD
    // Using aligned_alloc for better cache line alignment
    size_t alignment = 64;  // Cache line size
    size_t bufferSize = ((numStates * sizeof(double) + alignment - 1) / alignment) * alignment;
    
    clSum = (double*)aligned_alloc(alignment, bufferSize);
    clSumEnd = clSum + numStates;
    
    tempResult = (double*)aligned_alloc(alignment, bufferSize);
    
    if (clSum == nullptr || tempResult == nullptr) {
        Msg::error("Failed to allocate aligned memory for CondLikeJob buffers");
    }
}

CondLikeJob::~CondLikeJob(void) {

    free(clSum);
    free(tempResult);
}

void CondLikeJob::conditionalLikelihood(void) {

    // Dispatch to optimized or portable version based on platform
#   if defined(CONDLIKE_USE_ACCELERATE)
    conditionalLikelihoodOptimized();
#   elif defined(CONDLIKE_USE_CBLAS)
    conditionalLikelihoodOptimized();
#   else
    conditionalLikelihoodPortable();
#   endif
}

// optimized version using BLAS and vectorized math
void CondLikeJob::conditionalLikelihoodOptimized(void) {

#   if defined(CONDLIKE_USE_ACCELERATE) || defined(CONDLIKE_USE_CBLAS)
    lnScaler = 0.0;
    const int n = numStates;
    
    for (auto p = nodes.begin(); p != nodes.end(); p++) 
        {
        Node* parentNode = *p;
        
        // Initialize clSum to ones (multiplicative identity)
#       if defined(CONDLIKE_USE_ACCELERATE)
        double one = 1.0;
        vDSP_vfillD(&one, clSum, 1, (vDSP_Length)n);
#       else
        std::fill(clSum, clSum + n, 1.0);
#       endif
        
        int descendantCount = 0;
        
        // process each descendant
        for (Node* d = parentNode->getFirstChild(); d != nullptr; d = d->getNextSibling()) 
            {
            // prefetch next sibling's data while processing current
            Node* nextSibling = d->getNextSibling();
            if (nextSibling != nullptr) 
                {
                // prefetch transition probability matrix (n×n doubles)
                // temporal hint 1 (T1): keep in L2/L3, not L1
                double* nextP = nextSibling->getTransitionProbability()->begin();
                __builtin_prefetch(nextP, 0, 1);
                
                // prefetch conditional likelihood vector (n doubles)
                // temporal hint 0 (T0): keep in all cache levels
                double* nextCL = nextSibling->getConditionalLikelihood();
                __builtin_prefetch(nextCL, 0, 0);
                
                // for larger state spaces, prefetch additional cache lines
                // (each cache line is 64 bytes = 8 doubles)
                if (n > 8)
                    {
                    __builtin_prefetch(nextP + 8, 0, 1);
                    __builtin_prefetch(nextP + 16, 0, 1);
                    __builtin_prefetch(nextCL + 8, 0, 0);
                    }
                if (n > 16)
                    {
                    __builtin_prefetch(nextCL + 16, 0, 0);
                    }
                }
            
            // get transition probability matrix P (n x n, row-major)
            double* P = d->getTransitionProbability()->begin();
            
            // get descendants's conditional likelihood vector (length n)
            double* CL = d->getConditionalLikelihood();
            
            // compute tempResult = P * CL using BLAS
            cblas_dgemv(CblasRowMajor, CblasNoTrans,
                        n, n,           // matrix dimensions
                        1.0,            // alpha
                        P, n,           // matrix P with leading dimension n
                        CL, 1,          // vector CL with stride 1
                        0.0,            // beta
                        tempResult, 1); // output vector with stride 1
            
            // element-wise multiply: clSum = clSum .* tempResult
#           if defined(CONDLIKE_USE_ACCELERATE)
            vDSP_vmulD(clSum, 1, tempResult, 1, clSum, 1, (vDSP_Length)n);
#           else
            for (int i=0; i<n; i++)
                clSum[i] *= tempResult[i];
#           endif
            
            descendantCount++;
            
            /* Periodic rescaling to prevent underflow.
               With typical conditional likelihoods, we can safely multiply
               ~100-200 times before risking underflow near 10^-308 */
            if ((descendantCount & (rescaleFrequency - 1)) == 0) 
                {
                double maxVal;
#               if defined(CONDLIKE_USE_ACCELERATE)
                vDSP_maxvD(clSum, 1, &maxVal, (vDSP_Length)n);
#               else
                // cblas_idamax returns 0-based index of max |value|
                int maxIdx = cblas_idamax(n, clSum, 1);
                maxVal = clSum[maxIdx];
#               endif
                
                if (maxVal > 0.0) 
                    {
                    // Scale so largest element becomes 1.0
                    double invMax = 1.0 / maxVal;
                    cblas_dscal(n, invMax, clSum, 1);
                    lnScaler += log(maxVal);
                    }
                }
            }
        
        // final rescaling and store result
        double maxVal;
#       if defined(CONDLIKE_USE_ACCELERATE)
        vDSP_maxvD(clSum, 1, &maxVal, (vDSP_Length)n);
#       else
        int maxIdx = cblas_idamax(n, clSum, 1);
        maxVal = clSum[maxIdx];
#       endif
        
        double* clP = parentNode->getConditionalLikelihood();
        
        if (maxVal > 0.0) 
            {
            double invMax = 1.0 / maxVal;
#           if defined(CONDLIKE_USE_ACCELERATE)
            // vectorized multiply-and-copy: clP = clSum * invMax
            vDSP_vsmulD(clSum, 1, &invMax, clP, 1, (vDSP_Length)n);
#           else
            // copy then scale in place
            memcpy(clP, clSum, n * sizeof(double));
            cblas_dscal(n, invMax, clP, 1);
#           endif
            lnScaler += log(maxVal);
            } 
        else 
            {
            // all zeros (shouldn't happen with valid data)
            memcpy(clP, clSum, n * sizeof(double));
            }
        
        // resolve dependency for dependent job
        CondLikeJob* dependentJob = parentNode->getDependentJob();
        if (dependentJob != nullptr)
            dependentJob->resolveDependency();
        }
#   else
    // fallback if no BLAS available
    conditionalLikelihoodPortable();
#   endif
}

// portable version (no BLAS dependency)
void CondLikeJob::conditionalLikelihoodPortable(void) {

    lnScaler = 0.0;
    const int n = numStates;
    
    for (auto p = nodes.begin(); p != nodes.end(); p++)
        {
        Node* parentNode = *p;
        
        // initialize clSum to ones
        for (double* v=clSum; v != clSumEnd; v++)
            (*v) = 1.0;
        
        int descendantCount = 0;
            
        // process each descendant
        for (Node* d = parentNode->getFirstChild(); d != nullptr; d = d->getNextSibling())
            {
            // prefetch next sibling's data while processing current
            Node* nextSibling = d->getNextSibling();
            if (nextSibling != nullptr) 
                {
                double* nextP = nextSibling->getTransitionProbability()->begin();
                double* nextCL = nextSibling->getConditionalLikelihood();
                __builtin_prefetch(nextP, 0, 1);
                __builtin_prefetch(nextCL, 0, 0);
                if (n > 8)
                    {
                    __builtin_prefetch(nextP + 8, 0, 1);
                    __builtin_prefetch(nextP + 16, 0, 1);
                    __builtin_prefetch(nextCL + 8, 0, 0);
                    }
                if (n > 16)
                    {
                    __builtin_prefetch(nextCL + 16, 0, 0);
                    }
                }
            
            double* P = d->getTransitionProbability()->begin();
            double* CL = d->getConditionalLikelihood();
            
            // matrix-vector multiply: tempResult = P * CL
            // then multiply into clSum
            for (int i = 0; i < n; i++)
                {
                double sum = 0.0;
                double* Prow = P + i * n;
                
                for (int j = 0; j < n; j++)
                    {
                    sum += Prow[j] * CL[j];
                    }
                
                clSum[i] *= sum;
                }
            
            descendantCount++;
            
            // periodic rescaling to prevent underflow
            if ((descendantCount & (rescaleFrequency - 1)) == 0) 
                {
                double maxVal = clSum[0];
                for (double* v=clSum; v != clSumEnd; v++)
                    {
                    if (*v > maxVal)
                        maxVal = *v;
                    }
                
                if (maxVal > 0.0)
                    {
                    double invMax = 1.0 / maxVal;
                    for (double* v=clSum; v != clSumEnd; v++)
                        (*v) *= invMax;
                    lnScaler += log(maxVal);
                    }
                }
            }

        // final rescaling
        double maxVal = clSum[0];
        for (double* v=clSum; v != clSumEnd; v++)
            {
            if (*v > maxVal)
                maxVal = *v;
            }
        
        // store rescaled result
        double* clP = parentNode->getConditionalLikelihood();
        if (maxVal > 0.0)
            {
            double invMax = 1.0 / maxVal;
            for (int i = 0; i < n; i++)
                clP[i] = clSum[i] * invMax;
            lnScaler += log(maxVal);
            }
        else
            {
            for (int i = 0; i < n; i++)
                clP[i] = clSum[i];
            }
                    
        // resolve dependency
        CondLikeJob* dependentJob = parentNode->getDependentJob();
        if (dependentJob != nullptr)
            dependentJob->resolveDependency();
        }
}

void CondLikeJob::print(void) {

    std::cout << "Job " << jobId << std::endl;
    std::cout << "   Number of nodes : " << nodes.size() << std::endl;
    std::cout << "   Nodes           : ";
    for (auto p = nodes.begin(); p != nodes.end(); p++)
        std::cout << (*p)->getIndex() << " ";
    std::cout << std::endl;
    std::cout << "   No. Dependencies: " << numDependencies << std::endl;
}

void CondLikeJob::resolveDependency(void) {

    mtx.lock();
    numResolvedDependencies++;
    if (numResolvedDependencies == numDependencies)
        threadPool->pushTask(this);
    mtx.unlock();
}

void CondLikeJob::run(void) {

    conditionalLikelihood();
}
