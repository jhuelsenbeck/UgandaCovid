#ifndef ThreadPool_hpp
#define ThreadPool_hpp

#include <thread>
#include <mutex> 
#include <condition_variable>
#include <atomic>
#include <iostream>
#include "MathCache.hpp"



class ThreadTask {

    public:
                                ThreadTask(void);
        virtual                ~ThreadTask(void);
        virtual void            run(void) = 0;
};

class LikelihoodTask : public ThreadTask {

    public: 
        virtual                ~LikelihoodTask(void);
        
        void                    run(void) override {
                                    result = computeLnLikelihood();
                                }
        
        double                  getResult(void) const { return result; }
        
    protected:
        virtual double          computeLnLikelihood(void) = 0;
        
    private:
        double                  result;
};

class TransitionProbabilityTask : public ThreadTask {

    public:
        virtual                ~TransitionProbabilityTask(void);
        void                    run(void) override {
                                    computeTransitionProbabilities();
                                }
    
    protected:
        virtual void            computeTransitionProbabilities(void) = 0;
};

class ThreadPool {

    public:
        explicit                ThreadPool(void);
        explicit                ThreadPool(int numThreads);
                               ~ThreadPool(void);
        void                    pushTask(ThreadTask* task);
        void                    wait(void);
        int                     getThreadCount(void) const { return threadCount; }

    private:
        void                    worker(void);
        ThreadTask*             popTask(void);
        static constexpr size_t queueCapacity = 65536;          // cache-friendly circular buffer for tasks
        static constexpr size_t queueMask = queueCapacity - 1;  // note the size must be power of 2 for fast modulo via bitwise AND

                                // accessed on every push/pop operation
                                // group atomics together for cache coherence
        std::atomic<size_t>     tasksInFlight;                  // tasks queued + tasks being executed
        std::atomic<bool>       running;                        // 1 byte
        
                                // queue management (accessed under lock)
        size_t                  queueHead;                      // next slot to read from
        size_t                  queueTail;                      // next slot to write to
        size_t                  queueSize;                      // current number of tasks in queue
        
                                // accessed occasionally
        int                     threadCount;                    // 4 bytes
        std::thread*            threads;
        
                                // synchronization primitives
        std::mutex              mutex;                          // single mutex protects the queue
        std::condition_variable taskAvailable;                  // signaled when queue non-empty
        std::condition_variable allComplete;                    // signaled when tasksInFlight hits 0
        
                                // large array at end (8KB)
        ThreadTask*             taskQueue[queueCapacity];
};


#endif
