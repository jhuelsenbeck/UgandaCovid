#include "Model.hpp"
#include "Msg.hpp"
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
    
    // allocate matrices for the unique branch lengths
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (int i=0, n=(int)dpSeq.size(); i<n; i++)
        {
        Node* p = dpSeq[i];
        int v = p->getBrlen();
        ti_map::iterator it = tiMap.find(v);
        if (it == tiMap.end())
            {
            TransitionProbabilities* ti = new TransitionProbabilities(dim);
            ti->setBrlen(v);
            tiMap.insert( std::make_pair(v,ti) );
            }
        }
    std::cout << "     Tree has " << tiMap.size() << " unique branch lengths" << std::endl;
}

TransitionProbabilitiesMngr::~TransitionProbabilitiesMngr(void) {

    for (auto p : tiMap)
        {
        TransitionProbabilities* ti = p.second;
        delete ti;
        }
    tiMap.clear();
}

void TransitionProbabilitiesMngr::checkTiProbs(void) {

    for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
        {
        TransitionProbabilities* p = it->second;
        size_t n = p->dim();
        std::cout << "matrix " << it->first << std::endl;
        for (size_t i=0; i<n; i++)
            {
            double sum = 0.0;
            for (size_t j=0; j<n; j++)
                {
                double x = (*p)(i,j);
                sum += x;
                if (x < 0.0)
                    Msg::warning("Negative transition probability!!!!!");
                }
            std::cout << sum << std::endl;
            }
        }
}

TransitionProbabilities* TransitionProbabilitiesMngr::getTiProb(int brlen) {

    ti_map::iterator it = tiMap.find(brlen);
    if (it != tiMap.end())
        return it->second;
    return nullptr;
}

void TransitionProbabilitiesMngr::printMap(void) {

    for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
        {
        std::cout << it->first << " (" << it->second << ")" << std::endl;
        }
}

void TransitionProbabilitiesMngr::updateTransitionProbabilities(double rate) {

    RateMatrix* Q = modelPtr->getRateMatrix();

    auto tasks = new TransitionProbabilitiesTask[tiMap.size()];
    auto task = tasks;
    
    // update the main tree
    int id = 0;
    double r = modelPtr->getSubstitutionRate();
    for (ti_map::iterator it = tiMap.begin(); it != tiMap.end(); it++)
        {
        double v = (double)it->first * r;
        if (it->first == 0)
            v = MIN_BRLEN;
        task->init(id, (int)dim, v, (RealMatrix*)Q, it->second);
        threadPool->PushTask(task);
        ++task;
        ++id;
        }
        
    threadPool->Wait();
    delete [] tasks;
    
#   if 0
    checkTiProbs();
#   endif
}
