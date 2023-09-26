#ifndef threads_hpp
#define threads_hpp

#include <thread>
#include <mutex> 
#include <queue>
#include "MathCache.hpp"

class ThreadTask {

    public:
        ThreadTask(void);
        virtual ~ThreadTask(void);
        virtual void Run(MathCache& cache) = 0;
};


class ThreadPool {

    public:
        int  ThreadCount;

             explicit ThreadPool(void);
             ~ThreadPool(void);

        void PushTask(ThreadTask* task);
        void Wait();

    private:
        std::atomic<size_t>     TaskCount;
        std::atomic<bool>       Running;
        std::mutex              TaskMutex,
                                WaitMutex,
                                CheckMutex;
        std::condition_variable WaitCondition,
                                CheckCondition;
        std::queue<ThreadTask*> Tasks;
        std::thread*            Threads;

        void                    Worker(void);
        ThreadTask*             PopTask(void);
};


#endif
