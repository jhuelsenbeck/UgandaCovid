#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
#include "jph.hpp"
class TransitionProbabilityPair;



class Node {

    public:
                                    Node(void);
                                   ~Node(void);
        void                        addDescendant(Node* p) { descendants.insert(p); }
        Node*                       getAncestor(void) { return ancestor; }
        int                         getAreaId(void) { return areaId; }
        int                         getBrlen(void) { return brlen; }
        real*                       getConditionalLikelihood(void) { return cl; }
        real*                       getConditionalLikelihoodEnd(void) { return clEnd; }
        std::set<Node*>&            getDescendants(void) { return descendants; }
        int                         getIndex(void) { return index; }
        bool                        getIsTip(void) { return isTip; }
        std::string                 getName(void) { return name; }
        int                         getNumDescendants(void) { return (int)descendants.size(); }
        int                         getOffset(void) { return offset; }
        TransitionProbabilityPair*  getTransitionProbabilityPair(void) { return tiPair; }
        void                        removeAllDescendants(void) { descendants.clear(); }
        void                        removeDescendant(Node* p) { descendants.erase(p); }
        void                        setAncestor(Node* p) { ancestor = p; }
        void                        setAreaId(int x) { areaId = x; }
        void                        setBrlen(int x) { brlen = x; }
        void                        setConditionalLikelihood(real* p) { cl = p; }
        void                        setConditionalLikelihoodEnd(real* p) { clEnd = p; }
        void                        setIndex(int x) { index = x; }
        void                        setIsTip(bool tf) { isTip = tf; }
        void                        setName(std::string s) { name = s; }
        void                        setOffset(int x) { offset = x;}
        void                        setTransitionProbabilityPair(TransitionProbabilityPair* p) { tiPair = p; }
        
    private:
        Node*                       ancestor;
        real*                       cl;
        real*                       clEnd;
        TransitionProbabilityPair*  tiPair;
        int                         index;
        int                         offset;
        int                         brlen;
        int                         areaId;
        std::set<Node*>             descendants;
        std::string                 name;
        bool                        isTip;
};

#endif
