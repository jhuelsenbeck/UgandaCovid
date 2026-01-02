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

typedef std::map<std::pair<int,int>,TransitionProbabilities*> ti_map;
//typedef std::map<int,TransitionProbabilities*> ti_map;



class TransitionProbabilitiesMngr {

    public:
                                    TransitionProbabilitiesMngr(void) = delete;
                                    TransitionProbabilitiesMngr(const TransitionProbabilitiesMngr& m) = delete;
                                    TransitionProbabilitiesMngr(Model* m, Tree* t, size_t d, ThreadPool* tp);
                                   ~TransitionProbabilitiesMngr(void);
        TransitionProbabilities*    getTiProb(int brlen);
        void                        printMap(void);
        void                        updateTransitionProbabilities(double rate);
                
    private:
        void                        checkTiProbs(void);
        Model*                      modelPtr;
        ThreadPool*                 threadPool;
        size_t                      dim;
        ti_map                      tiMap;
        std::pair<int,int>          key;
};

#endif
