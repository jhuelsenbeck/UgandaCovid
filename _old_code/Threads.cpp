#include "Threads.hpp"

ThreadTask::ThreadTask(void) {

}

ThreadTask::~ThreadTask(void) {

}


#if 1
// Threaded version

ThreadPool::ThreadPool(void):
    threadCount(std::thread::hardware_concurrency()),
    taskCount(0),
    running(true),
    threads(new std::thread[threadCount]) {
    
    for (int i = 0; i < threadCount; i++)
        threads[i] = std::thread(&ThreadPool::worker, this);
}

ThreadPool::ThreadPool(int n):
    threadCount(n),
    taskCount(0),
    running(true),
    threads(new std::thread[threadCount]) {
    
    for (int i = 0; i < threadCount; i++)
        threads[i] = std::thread(&ThreadPool::worker, this);
}

ThreadPool::~ThreadPool(void) {

    wait();
    if (threads)
        {
        running = false;
        taskCount = 1;
        checkCondition.notify_all();
        for (auto* t = threads; t < threads + threadCount; ++t)
            t->join();
        delete[] threads;
        }
}

void ThreadPool::pushTask(ThreadTask* task) {

    {
    std::lock_guard<std::mutex> lock(taskMutex);
    ++taskCount;
    tasks.push(task);
    }
    std::unique_lock mlock(checkMutex);
    checkCondition.notify_one();
}
    
ThreadTask* ThreadPool::popTask(void) {

    {
    std::unique_lock mlock(checkMutex);
    checkCondition.wait(mlock, [this]{return taskCount > 0;});
    }

    std::lock_guard<std::mutex> tasklock(taskMutex);
    if (tasks.empty())
        return NULL;
    else 
        {
        auto task = tasks.front();
        tasks.pop();
        return task;
        }
}

void ThreadPool::wait(void) {

    for (;;) 
        {
        if (taskCount == 0)
            break;
        else
            std::this_thread::yield();
        }
}

void ThreadPool::worker(void) {

    MathCache cache;
    while (running)
        {
        ThreadTask* task = popTask();
        if (task)
            {
            task->run(cache);
            --taskCount;
            }
        std::this_thread::yield();
        }
}

#else

// Serial version
ThreadPool::ThreadPool(void):
    ThreadCount(1),
    TaskCount(0),
    Running(false),
    Threads(NULL) {

}

ThreadPool::~ThreadPool(void) {

}

void ThreadPool::pushTask(ThreadTask* task) {

    task->run();
}

ThreadTask* ThreadPool::popTask(void) {

    return NULL;
}

void ThreadPool::wait(void) {

}

void ThreadPool::worker(void) {

}

#endif

