#include "Model.hpp"
#include "Node.hpp"
#include "Threads.hpp"
#include "TransitionProbabilities.hpp"
#include "TransitionProbabilitiesMngr.hpp"
#include "Tree.hpp"

#define MIN_BRLEN 10e-5



TransitionProbabilitiesMngr::TransitionProbabilitiesMngr(Model* m, Tree* t, size_t d, ThreadPool* tp) {

    // pointer to the model object
    modelPtr = m;
    
    threadPool = tp;
    
    // set the dimensions of the rate matrices
    dim = d;
    
    // set the active transition probabilities to be space 0
    activeTiProb = 0;
    
    // allocate matrices for the unique branch lengths
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (int i=0, n=(int)dpSeq.size(); i<n; i++)
        {
        Node* p = dpSeq[i];
        int v = p->getBrlen();
        ti_map::iterator it = tiMap.find(v);
        if (it == tiMap.end())
            {
            TransitionProbabilities* ti1 = new TransitionProbabilities(dim);
            TransitionProbabilities* ti2 = new TransitionProbabilities(dim);
            TransitionProbabilityPair tiPair;
            tiPair.tipr[0] = ti1;
            tiPair.tipr[1] = ti2;
            ti1->setBrlen(v);
            ti2->setBrlen(v);
            tiMap.insert( std::make_pair(v,tiPair) );
            }
        }
}

TransitionProbabilitiesMngr::~TransitionProbabilitiesMngr(void) {

    for (auto p : tiMap)
        {
        TransitionProbabilities* ti1 = p.second.tipr[0];
        TransitionProbabilities* ti2 = p.second.tipr[1];
        delete ti1;
        delete ti2;
        }
    tiMap.clear();
}

TransitionProbabilityPair* TransitionProbabilitiesMngr::getTiPair(int brlen) {

    ti_map::iterator it = tiMap.find(brlen);
    if (it != tiMap.end())
        return &(it->second);
    return nullptr;
}

TransitionProbabilities* TransitionProbabilitiesMngr::getTiProb(int brlen) {

    ti_map::iterator it = tiMap.find(brlen);
    if (it != tiMap.end())
        return it->second.tipr[activeTiProb];
    return nullptr;
}

void TransitionProbabilitiesMngr::switchActive(void) {

    if (activeTiProb == 0)
        activeTiProb = 1;
    else
        activeTiProb = 0;
}

void TransitionProbabilitiesMngr::updateTransitionProbabilities(double rate) {

    RateMatrix* Q = modelPtr->getRateMatrix();

    auto tasks = new TransitionProbabilitiesTask[tiMap.size()];
    auto task = tasks;
    
    // update the main tree
    int id = 0;
    for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
        {
        double v = (double)it->first;
        if (it->first == 0)
            v = MIN_BRLEN;
        task->init(id, (int)dim, v, (DoubleMatrix*)Q, it->second.tipr[activeTiProb]);
        threadPool->PushTask(task);
        ++task;
        ++id;
        }
        
    threadPool->Wait();
    delete [] tasks;
}
