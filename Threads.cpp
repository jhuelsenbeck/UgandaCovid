#include "Threads.hpp"

ThreadTask::ThreadTask(void) {

}

ThreadTask::~ThreadTask(void) {

}


#if 1
// Threaded version

ThreadPool::ThreadPool(void):
    ThreadCount(std::thread::hardware_concurrency()),
    TaskCount(0),
    Running(true),
    Threads(new std::thread[ThreadCount])
{
    for (int i = 0; i < ThreadCount; i++)
        Threads[i] = std::thread(&ThreadPool::Worker, this);
}

ThreadPool::~ThreadPool(void) {

    Wait();
    if (Threads)
        {
        Running = false;
        TaskCount = 1;
        CheckCondition.notify_all();
        for (auto* t = Threads; t < Threads + ThreadCount; ++t)
            t->join();
        delete[] Threads;
        }
}

void ThreadPool::PushTask(ThreadTask* task) {

    {
      std::lock_guard<std::mutex> lock(TaskMutex);
      ++TaskCount;
      Tasks.push(task);
    }
    std::unique_lock mlock(CheckMutex);
    CheckCondition.notify_one();
}
    
ThreadTask* ThreadPool::PopTask(void) {

    {
        std::unique_lock mlock(CheckMutex);
        CheckCondition.wait(mlock, [this]{return TaskCount > 0;});
    }

    std::lock_guard<std::mutex> tasklock(TaskMutex);
    if (Tasks.empty())
        return NULL;
    else 
        {
        auto task = Tasks.front();
        Tasks.pop();
        return task;
        }
}

void ThreadPool::Wait(void) {

    for (;;) 
        {
        // This WaitCondition is signaled when all tasks are completed
//        std::unique_lock lock(WaitMutex);
//        WaitCondition.wait(lock);

        // According to the documentation, this needs to be double-checked because of "spurious" wakeups

//        size_t count;
//        {
//            std::lock_guard<std::mutex> tasklock(TaskMutex);
//            count = TaskCount;
//        }
        if (TaskCount == 0)
            break;
        else
            std::this_thread::yield();
        }
}

void ThreadPool::Worker(void) {

    MathCache cache;
    while (Running)
        {
        ThreadTask* task = PopTask();
        if (task)
            {
            task->Run(cache);

            --TaskCount;
//            if (TaskCount == 0)

  //          size_t count;
//            {
//                std::lock_guard<std::mutex> tasklock(TaskMutex);
//                count = --TaskCount;
//            }
//            if (count == 0)
//                WaitCondition.notify_one();
            }
//        else
            std::this_thread::yield();
    }
}

#else

// Serial version

ThreadPool::ThreadPool():
    ThreadCount(1),
    TaskCount(0),
    Running(false),
    Threads(NULL)
{
}

ThreadPool::~ThreadPool(void) {
}

void ThreadPool::PushTask(ThreadTask* task) {

    task->Run();
}

ThreadTask* ThreadPool::PopTask(void) {
    return NULL;
}

void ThreadPool::Wait(void) {

}

void ThreadPool::Worker(void) {

}

#endif

