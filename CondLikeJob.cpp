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



CondLikeJob::CondLikeJob(ThreadPool* tp, int na) : 
    ThreadTask(), threadPool(tp), numStates(na) {

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

// -------------------------------------------------------------------
// Optimized version using BLAS and vectorized math
// -------------------------------------------------------------------
// Uses:
// - cblas_dgemv for matrix-vector multiply (P * CL)
// - vvlog (Accelerate only) for vectorized logarithm
// - Efficient memory access patterns
// -------------------------------------------------------------------

void CondLikeJob::conditionalLikelihoodOptimized(void) {

#   if defined(CONDLIKE_USE_ACCELERATE) || defined(CONDLIKE_USE_CBLAS)
    lnScaler = 0.0;
    const int n = numStates;
    
    for (auto p = nodes.begin(); p != nodes.end(); p++) {
        Node* parentNode = *p;
        
        // Initialize clSum to zeros
        memset(clSum, 0, n * sizeof(double));
        
        // Process each child
        for (Node* d = parentNode->getFirstChild(); d != nullptr; d = d->getNextSibling()) {
            // Get transition probability matrix P (n x n, row-major)
            double* P = d->getTransitionProbability()->begin();
            
            // Get child's conditional likelihood vector (length n)
            double* CL = d->getConditionalLikelihood();
            
            // Compute tempResult = P * CL using BLAS
            // This replaces the manual double loop with optimized SIMD code
            cblas_dgemv(CblasRowMajor, CblasNoTrans,
                        n, n,           // matrix dimensions
                        1.0,            // alpha
                        P, n,           // matrix P with leading dimension n
                        CL, 1,          // vector CL with stride 1
                        0.0,            // beta
                        tempResult, 1); // output vector with stride 1
            
            // Add log(tempResult) to clSum
#           if defined(CONDLIKE_USE_ACCELERATE)
            // Use Accelerate's vectorized log for all n values at once
            // vvlog computes log(x) for a vector - much faster than scalar loop
            int nn = n;
            vvlog(tempResult, tempResult, &nn);  // in-place log
            
            // Vectorized addition: clSum += tempResult
            cblas_daxpy(n, 1.0, tempResult, 1, clSum, 1);
#           else
            // OpenBLAS path: scalar log (still benefits from dgemv)
            for (int i = 0; i < n; i++) {
                clSum[i] += log(tempResult[i]);
            }
#           endif
        }
        
        // Find maximum for rescaling (prevents underflow)
        double largestLogProb = clSum[0];
        for (int i = 1; i < n; i++) {
            if (clSum[i] > largestLogProb)
                largestLogProb = clSum[i];
        }
        
        // Rescale and convert back to probabilities
        double* clP = parentNode->getConditionalLikelihood();
        
#       if defined(CONDLIKE_USE_ACCELERATE)
        // Vectorized: clSum = clSum - largestLogProb, then exp
        double negLargest = -largestLogProb;
        int nn = n;
        
        // clSum += (-largestLogProb) * ones  =>  clSum[i] -= largestLogProb
        // Use vDSP for this
        vDSP_vsaddD(clSum, 1, &negLargest, clP, 1, n);
        
        // Vectorized exp
        vvexp(clP, clP, &nn);
#       else
        // Scalar fallback
        for (int i = 0; i < n; i++) {
            clP[i] = exp(clSum[i] - largestLogProb);
        }
#       endif
        
        lnScaler += largestLogProb;
        
        // Resolve dependency for dependent job
        CondLikeJob* dependentJob = parentNode->getDependentJob();
        if (dependentJob != nullptr)
            dependentJob->resolveDependency();
    }
#   else
    // Fallback if no BLAS available
    conditionalLikelihoodPortable();
#   endif
}

// -------------------------------------------------------------------
// Portable version (no BLAS dependency)
// -------------------------------------------------------------------
// Still optimized with:
// - Better memory access patterns (row-major traversal)
// - Loop structure that helps auto-vectorization
// - Reduced pointer arithmetic in inner loop
// -------------------------------------------------------------------

void CondLikeJob::conditionalLikelihoodPortable(void) {

    lnScaler = 0.0;
    const int n = numStates;
    
    for (auto p = nodes.begin(); p != nodes.end(); p++) 
        {
        Node* parentNode = *p;
        
        // initialize clSum to zeros
        for (int i=0; i<n; i++)
            clSum[i] = 0.0;
            
        // process each descendant
        for (Node* d = parentNode->getFirstChild(); d != nullptr; d = d->getNextSibling()) 
            {
            double* P = d->getTransitionProbability()->begin();
            double* CL = d->getConditionalLikelihood();
            
            // matrix-vector multiply: tempResult = P * CL
            // written to help auto-vectorization
            for (int i=0; i<n; i++) 
                {
                double sum = 0.0;
                double* Prow = P + i * n;  // Pointer to row i of P
                
                // inner loop - compiler can auto-vectorize this
                for (int j = 0; j < n; j++) 
                    {
                    sum += Prow[j] * CL[j];
                    }
                
                clSum[i] += log(sum);
                }
            }

        // find maximum for rescaling
        double largestLogProb = clSum[0];
        for (int i = 1; i < n; i++) 
            {
            if (clSum[i] > largestLogProb)
                largestLogProb = clSum[i];
            }
        
        // Rescale and convert to probabilities
        double* clP = parentNode->getConditionalLikelihood();
        for (int i = 0; i < n; i++) 
            {
            clP[i] = exp(clSum[i] - largestLogProb);
            }
        
        lnScaler += largestLogProb;
                    
        // Resolve dependency
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
