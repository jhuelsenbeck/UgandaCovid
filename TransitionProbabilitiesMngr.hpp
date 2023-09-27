#ifndef TransitionProbabilitiesMngr_hpp
#define TransitionProbabilitiesMngr_hpp

#include <iostream>
#include <map>
class Model;
class ThreadPool;
class TransitionProbabilities;
class Tree;

struct TransitionProbabilityPair {

    TransitionProbabilities*    tipr[2];
};

typedef std::map<int,TransitionProbabilityPair> ti_map;



class TransitionProbabilitiesMngr {

    public:
                                    TransitionProbabilitiesMngr(void) = delete;
                                    TransitionProbabilitiesMngr(const TransitionProbabilitiesMngr& m) = delete;
                                    TransitionProbabilitiesMngr(Model* m, Tree* t, size_t d, ThreadPool* tp);
                                   ~TransitionProbabilitiesMngr(void);
        int                         getActiveTiProb(void) { return activeTiProb; }
        TransitionProbabilities*    getTiProb(int brlen);
        TransitionProbabilityPair*  getTiPair(int brlen);
        void                        switchActive(void);
        void                        updateTransitionProbabilities(real rate);
                
    private:
        Model*                      modelPtr;
        ThreadPool*                 threadPool;
        size_t                      dim;
        int                         activeTiProb;
        ti_map                      tiMap;
};

#endif
