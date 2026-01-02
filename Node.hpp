#ifndef Node_hpp
#define Node_hpp

#include <string>
class CondLikeJob;
class History;
class TransitionProbabilities;


// -------------------------------------------------------------------
// LCRS (Left-Child, Right-Sibling) tree representation
// -------------------------------------------------------------------
// Instead of std::set<Node*> descendants (40 bytes overhead per child),
// we use two pointers per node:
//   - firstChild: leftmost child (nullptr if leaf)
//   - nextSibling: next sibling in parent's child list (nullptr if last)
//
// This is O(1) memory overhead regardless of child count, and enables
// simple iteration: for (Node* d = firstChild; d; d = d->nextSibling)
// -------------------------------------------------------------------

class Node {

    public:
                                    Node(void);
                                   ~Node(void);
        
        // --- LCRS tree structure ---
        void                        addDescendant(Node* p);
        void                        removeDescendant(Node* p);
        void                        removeAllDescendants(void);
        Node*                       getAncestor(void) { return ancestor; }
        Node*                       getFirstChild(void) { return firstChild; }
        Node*                       getNextSibling(void) { return nextSibling; }
        int                         getNumDescendants(void);
        
        // Legacy compatibility - returns first child
        Node*                       getFirstDescendant(void) { return firstChild; }
        
        // --- Node properties ---
        int                         getAreaId(void) { return areaId; }
        std::string                 getAreaName(void) { return areaName; }
        bool                        getIsAreaFixed(void) { return isAreaFixed; }
        int                         getBrlen(void) { return brlen; }
        double                      getBrlenExact(void) { return brlenExact; }
        double*                     getConditionalLikelihood(void) { return cl; }
        double*                     getConditionalLikelihoodEnd(void) { return clEnd; }
        CondLikeJob*                getDependentJob(void) { return dependentJob; }
        bool                        getGoodTime(void) { return goodTime; }
        History*                    getHistory(void) { return history; }
        int                         getIndex(void) { return index; }
        int                         getIntervalIdx(void) { return intervalIdx; }
        bool                        getIsTip(void) { return isTip; }
        CondLikeJob*                getJob(void) { return job; }
        std::string                 getName(void) { return name; }
        int                         getOffset(void) { return offset; }
        double                      getOldestDescendant(void);
        bool                        getScratchBool(void) { return scratchBool; }
        int                         getScratchInt(void) { return scratchInt; }
        double                      getTime(void) { return time; }
        TransitionProbabilities*    getTransitionProbability(void) { return tiProb; }
        
        void                        setAncestor(Node* p) { ancestor = p; }
        void                        setFirstChild(Node* p) { firstChild = p; }
        void                        setNextSibling(Node* p) { nextSibling = p; }
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
        // LCRS tree pointers (24 bytes total for tree structure)
        Node*                       ancestor;       // parent node
        Node*                       firstChild;     // leftmost child (nullptr if leaf)
        Node*                       nextSibling;    // next sibling (nullptr if last child)
        
        // Computational data pointers
        double*                     cl;
        double*                     clEnd;
        TransitionProbabilities*    tiProb;
        History*                    history;
        
        // Node identifiers
        int                         index;
        int                         offset;
        
        // Branch length (discretized and exact)
        int                         brlen;
        double                      brlenExact;
        
        // Time and interval
        double                      time;
        int                         intervalIdx;
        
        // Geographic state
        int                         areaId;
        std::string                 areaName;
        bool                        isAreaFixed;
        
        // Node metadata
        std::string                 name;
        bool                        isTip;
        
        // Scratch variables for algorithms
        int                         scratchInt;
        bool                        scratchBool;
        bool                        goodTime;
        
        // Job scheduling
        CondLikeJob*                job;
        CondLikeJob*                dependentJob;
};

#endif
