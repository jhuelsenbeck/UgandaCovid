#ifndef Threads_hpp
#define Threads_hpp

#include <condition_variable>
#include <mutex>
#include <queue>
#include <thread>
#include "MathCache.hpp"



class ThreadTask {

    public:
                                ThreadTask(void);
        virtual                ~ThreadTask(void);
        virtual void            run(MathCache& cache) = 0;
};

class ThreadPool {

    public:
        explicit                ThreadPool(void);
        explicit                ThreadPool(int n);
                               ~ThreadPool(void);
        void                    pushTask(ThreadTask* task);
        void                    wait(void);
        int                     threadCount;

    private:
        void                    worker(void);
        ThreadTask*             popTask(void);
        std::atomic<size_t>     taskCount;
        std::atomic<bool>       running;
        std::mutex              taskMutex,
                                waitMutex,
                                checkMutex;
        std::condition_variable waitCondition,
                                checkCondition;
        std::queue<ThreadTask*> tasks;
        std::thread*            threads;
};


#endif
