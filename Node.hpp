#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
class CondLikeJob;
class History;
class TransitionProbabilities;



class Node {

    public:
                                    Node(void);
                                   ~Node(void);
        void                        addDescendant(Node* p) { descendants.insert(p); }
        Node*                       getAncestor(void) { return ancestor; }
        int                         getAreaId(void) { return areaId; }
        std::string                 getAreaName(void) { return areaName; }
        bool                        getIsAreaFixed(void) { return isAreaFixed; }
        int                         getBrlen(void) { return brlen; }
        double                      getBrlenExact(void) { return brlenExact; }
        double*                     getConditionalLikelihood(void) { return cl; }
        double*                     getConditionalLikelihoodEnd(void) { return clEnd; }
        CondLikeJob*                getDependentJob(void) { return dependentJob; }
        std::set<Node*>&            getDescendants(void) { return descendants; }
        Node*                       getFirstDescendant(void);
        bool                        getGoodTime(void) { return goodTime; }
        History*                    getHistory(void) { return history; }
        int                         getIndex(void) { return index; }
        bool                        getIsTip(void) { return isTip; }
        CondLikeJob*                getJob(void) { return job; }
        std::string                 getName(void) { return name; }
        int                         getNumDescendants(void) { return (int)descendants.size(); }
        int                         getOffset(void) { return offset; }
        double                      getOldestDescendant(void);
        bool                        getScratchBool(void) { return scratchBool; }
        int                         getScratchInt(void) { return scratchInt; }
        double                      getTime(void) { return time; }
        TransitionProbabilities*    getTransitionProbability(void) { return tiProb; }
        void                        removeAllDescendants(void) { descendants.clear(); }
        void                        removeDescendant(Node* p) { descendants.erase(p); }
        void                        setAncestor(Node* p) { ancestor = p; }
        void                        setAreaId(int x) { areaId = x; }
        void                        setAreaName(std::string s) { areaName = s; }
        void                        setIsAreaFixed(bool tf) { isAreaFixed = tf; }
        void                        setBrlen(int x) { brlen = x; }
        void                        setBrlenExact(double x) { brlenExact = x; }
        void                        setConditionalLikelihood(double* p) { cl = p; }
        void                        setConditionalLikelihoodEnd(double* p) { clEnd = p; }
        void                        setDependentJob(CondLikeJob* j) { dependentJob = j; }
        void                        setGoodTime(bool tf) { goodTime = tf; }
        void                        setHistory(History* h) { history = h; }
        void                        setIndex(int x) { index = x; }
        void                        setIntervalIdx(int x) { intervalIdx = x; }
        void                        setIsTip(bool tf) { isTip = tf; }
        void                        setJob(CondLikeJob* j) { job = j; }
        void                        setName(std::string s);
        void                        setOffset(int x) { offset = x;}
        void                        setScratchBool(bool tf) { scratchBool = tf; }
        void                        setScratchInt(int x) { scratchInt = x; }
        void                        setTime(double x) { time = x; }
        void                        setTransitionProbability(TransitionProbabilities* p) { tiProb = p; }
        
    private:
        Node*                       ancestor;
        double*                     cl;
        double*                     clEnd;
        TransitionProbabilities*    tiProb;
        History*                    history;
        int                         index;
        int                         offset;
        int                         brlen;
        double                      brlenExact;
        double                      time;
        int                         areaId;
        std::string                 areaName;
        bool                        isAreaFixed;
        std::set<Node*>             descendants;
        std::string                 name;
        bool                        isTip;
        int                         scratchInt;
        bool                        scratchBool;
        bool                        goodTime;
        CondLikeJob*                job;
        CondLikeJob*                dependentJob;
        int                         intervalIdx;
};

#endif
