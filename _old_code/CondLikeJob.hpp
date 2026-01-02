#ifndef CondLikeJob_hpp
#define CondLikeJob_hpp

#include <mutex>
#include <set>
#include <vector>
#include "Threads.hpp"
class CondLikeJobMngr;
class Node;
class Tree;



class CondLikeJob : public ThreadTask {

    public:
                                    CondLikeJob(void) = delete;
                                    CondLikeJob(CondLikeJobMngr* m, ThreadPool* tp, int na);
                                   ~CondLikeJob(void);
        void                        addDependency(CondLikeJob* j) { numDependencies++; }
        void                        addNode(Node* p) { nodes.push_back(p); }
        int                         getJobId(void) { return jobId; }
        std::vector<Node*>&         getJobNodes(void) { return nodes; }
        double                      getLnScaler(void) { return lnScaler; }
        int                         getNumDependencies(void) { return numDependencies; }
        int                         numNodesInJob(void) { return (int)nodes.size(); }
        void                        print(void);
        void                        resolveDependency(void);
        virtual void                run(MathCache& cache);
        void                        setJobId(int x) { jobId = x; }
        void                        setNumDepencies(int x) { numDependencies = x; }
        void                        setNumResolvedDependencies(int x) { numResolvedDependencies = x; }
    
    private:
        void                        conditionalLikelihood(void);
        double*                     clSum;
        double*                     clSumEnd;
        CondLikeJobMngr*            myManager;
        ThreadPool*                 threadPool;
        std::mutex                  mtx;
        double                      lnScaler;
        int                         jobId;
        std::vector<Node*>          nodes;
        int                         numDependencies;
        int                         numResolvedDependencies;
        int                         numStates;
};

#endif
