#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
class TransitionProbabilityPair;



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
        void                        setConditionalLikelihood(double* p) { cl = p; }
        void                        setConditionalLikelihoodEnd(double* p) { clEnd = p; }
        void                        setIndex(int x) { index = x; }
        void                        setIsTip(bool tf) { isTip = tf; }
        void                        setName(std::string s) { name = s; }
        void                        setOffset(int x) { offset = x;}
        void                        setTransitionProbabilityPair(TransitionProbabilityPair* p) { tiPair = p; }
        
    private:
        Node*                       ancestor;
        double*                     cl;
        double*                     clEnd;
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
