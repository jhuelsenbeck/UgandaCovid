#ifndef CondLikeJobMngr_hpp
#define CondLikeJobMngr_hpp

#include <mutex>
#include <vector>
class CondLikeJob;
class ThreadPool;
class Tree;



class CondLikeJobMngr {

    public:
                                    CondLikeJobMngr(void) = delete;
                                    CondLikeJobMngr(const CondLikeJobMngr& m) = delete;
                                    CondLikeJobMngr(Tree* t, ThreadPool* tp, int na);
                                   ~CondLikeJobMngr(void);
        void                        calculateConditionalLikelihoods(void);
        void                        calculateTime(int n);
        double                      getScaler(void);
        Tree*                       getTree(void) { return tree; }
        void                        print(void);
        void                        lock(void) { mtx.lock(); }
        void                        unlock(void) { mtx.unlock(); }
                
    private:
        CondLikeJob*                addJob(void);
        void                        jobSizeDistribution(void);
        void                        zeroResolvedDependencies(void);
        ThreadPool*                 threadPool;
        Tree*                       tree;
        std::vector<CondLikeJob*>   jobs;
        int                         jobId;
        int                         numStates;
        std::mutex                  mtx;
};

#endif

