#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
class CondLikeJob;
class TransitionProbabilities;



class Node {

    public:
                                    Node(void);
                                   ~Node(void);
        void                        addDescendant(Node* p) { descendants.insert(p); }
        Node*                       getAncestor(void) { return ancestor; }
        int                         getAreaId(void) { return areaId; }
        int                         getBrlen(void) { return brlen; }
        double*                     getConditionalLikelihood(void) { return cl; }
        double*                     getConditionalLikelihoodEnd(void) { return clEnd; }
        CondLikeJob*                getDependentJob(void) { return dependentJob; }
        std::set<Node*>&            getDescendants(void) { return descendants; }
        int                         getIndex(void) { return index; }
        bool                        getIsTip(void) { return isTip; }
        CondLikeJob*                getJob(void) { return job; }
        std::string                 getName(void) { return name; }
        int                         getNumDescendants(void) { return (int)descendants.size(); }
        int                         getOffset(void) { return offset; }
        int                         getScratch(void) { return scratch; }
        TransitionProbabilities*    getTransitionProbability(void) { return tiProb; }
        void                        removeAllDescendants(void) { descendants.clear(); }
        void                        removeDescendant(Node* p) { descendants.erase(p); }
        void                        setAncestor(Node* p) { ancestor = p; }
        void                        setAreaId(int x) { areaId = x; }
        void                        setBrlen(int x) { brlen = x; }
        void                        setConditionalLikelihood(double* p) { cl = p; }
        void                        setConditionalLikelihoodEnd(double* p) { clEnd = p; }
        void                        setDependentJob(CondLikeJob* j) { dependentJob = j; }
        void                        setIndex(int x) { index = x; }
        void                        setIsTip(bool tf) { isTip = tf; }
        void                        setJob(CondLikeJob* j) { job = j; }
        void                        setName(std::string s);
        void                        setOffset(int x) { offset = x;}
        void                        setScratch(int x) { scratch = x; }
        void                        setTransitionProbability(TransitionProbabilities* p) { tiProb = p; }
        
    private:
        Node*                       ancestor;
        double*                     cl;
        double*                     clEnd;
        TransitionProbabilities*    tiProb;
        int                         index;
        int                         offset;
        int                         brlen;
        int                         areaId;
        std::set<Node*>             descendants;
        std::string                 name;
        bool                        isTip;
        int                         scratch;
        CondLikeJob*                job;
        CondLikeJob*                dependentJob;
};

#endif
