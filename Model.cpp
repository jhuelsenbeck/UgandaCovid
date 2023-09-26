#include <chrono>
#include <iomanip>
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "RateMatrix.hpp"
#include "TransitionProbabilities.hpp"
#include "TransitionProbabilitiesMngr.hpp"
#include "Tree.hpp"



Model::Model(Tree* tp, RateMatrix* m, ThreadPool* thp) {

    threadPool = thp;
    numStates = m->getNumRows();

    // set up the model parameters
    activeSubstitutionRate = 0;
    substitutionRate[0] = 0.1; // pick a good starting value?
    substitutionRate[1] = substitutionRate[0];
    activePi = 0;
    pi[0].resize(numStates);
    pi[1].resize(numStates);
    std::fill(pi[0].begin(), pi[0].end(), 1.0/numStates);
    pi[1] = pi[0];
    activeR = 0;
    r[0].resize(numStates*(numStates-1)/2);
    std::fill(r[0].begin(), r[0].end(), 1.0);
    r[1] = r[0];

    // set up the rate matrix and tree
    tree = tp;
    activeRateMatrix = 0;
    q[0] = m;
    q[1] = new RateMatrix(*m);
    q[0]->set(&pi[0][0], &r[0][0]);
    (*q[1]) = (*q[0]);
    //std::cout << *q[0] << std::endl;
    //std::cout << *q[1] << std::endl;
    
    // set up the transition probabilities
    tiMngr = new TransitionProbabilitiesMngr(this, tree, numStates, threadPool);
    tiMngr->updateTransitionProbabilities(0.1);
        
    // set up the conditional likelihoods (must be done after tiMngr is instantiated
    initializeConditionalLikelihoods();
}

Model::~Model(void) {

    delete tree;
    delete q[0];
    delete q[1];
    delete tiMngr;
}

void Model::initializeConditionalLikelihoods(void) {
        
    size_t numAllocatedDoubles = 0;
    size_t memSize = numStates * sizeof(double);
    
    std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    for (auto p=dpSeq.begin(); p != dpSeq.end(); p++)
        {
        // allocate and initialize the conditional likelihoods
        double* condLikes = new double[memSize];
        if (!condLikes)
            Msg::error("Failure to allocate conditional likelihood for node " + std::to_string((*p)->getIndex()));
            
        for (int i=0; i<memSize; i++)
            condLikes[i] = 0.0;
        if ((*p)->getIsTip() == true)
            {
            int areaId = (*p)->getAreaId();
            if (areaId != -1)
                {
                condLikes[areaId] = 1.0;
                }
            else
                {
                for (int i=0; i<numStates; i++)
                    condLikes[i] = 1.0;
                }
            }
        (*p)->setConditionalLikelihood(condLikes);
        (*p)->setConditionalLikelihoodEnd(condLikes+numStates);
        numAllocatedDoubles += numStates;
        
        // set up the transition probability pair for the node
        TransitionProbabilityPair* tip = tiMngr->getTiPair((*p)->getBrlen());
        if (tip == nullptr)
            Msg::error("Could not find transition probabilities for node " + std::to_string((*p)->getIndex()));
        (*p)->setTransitionProbabilityPair(tip);
        }
    std::cout << "Successfully allocated a total of " << numAllocatedDoubles << " doubles for the conditional likelihoods" << std::endl;

    double lnL = lnLikelihood();
    std::cout << "lnL = " << lnL << std::endl;
}

#if 1

// serial version of likelihood calculation
double Model::lnLikelihood(void) {

    auto begin = std::chrono::high_resolution_clock::now();

    // allocate a temporary vector for the sums
    double* lnSum = (double*)malloc(numStates*sizeof(double));
    double* lnSumEnd = lnSum + numStates;
    
    // calculate the conditional likelihoods
    int activeTi = tiMngr->getActiveTiProb();
    std::vector<Node*>& dpSeq = tree->getInteriorDownPassSequence();
    double lnFactor = 0.0;
    for (auto p=dpSeq.begin(); p != dpSeq.end(); p++)
        {
        for (double* v=lnSum; v != lnSumEnd; v++)
            (*v) = 0.0;
            
        std::set<Node*>& pDesc = (*p)->getDescendants();
        for (Node* d : pDesc)
            {
            double* P = d->getTransitionProbabilityPair()->tipr[activeTi]->begin();
            double* clBegin = d->getConditionalLikelihood();
            double* clEnd = d->getConditionalLikelihoodEnd();
            for (double* s=lnSum; s != lnSumEnd; s++)
                {
                double sum = 0.0;
                for (double* L=clBegin; L != clEnd; L++)
                    {
                    sum += (*P) * (*L);
                    P++;
                    }
                *s += log(sum);
                }
            }
            
        // rescale
        double lnLargestCl = lnSum[0];
        for (double* v=&lnSum[1]; v != lnSumEnd; v++)
            {
            if (*v > lnLargestCl)
                lnLargestCl = *v;
            }
        lnFactor += lnLargestCl;
        for (double* v=lnSum; v != lnSumEnd; v++)
            *v -= lnLargestCl;

        double* LP = (*p)->getConditionalLikelihood();
        for (double* v=lnSum; v != lnSumEnd; v++)
            {
            *LP = exp(*v);
            LP++;
            }

#       if 0
        for (int i=0; i<numStates; i++)
            std::cout << "cl[" << i << "] = " << std::fixed << std::setprecision(50) << LP[i] << std::endl;
        getchar();
#       endif
        }
        
    // calculate likelihood
    double* LR = tree->getRoot()->getConditionalLikelihood();
    double* f = &pi[activePi][0];
    double like = 0.0;
    for (int i=0; i<numStates; i++)
        {
        like += (*f) * (*LR);
        LR++;
        f++;
        }
        
    free(lnSum);

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "s" << std::endl;

    return log(like) + lnFactor;
}

#else

// threaded version
double Model::lnLikelihood(void) {

    auto begin = std::chrono::high_resolution_clock::now();

    // allocate a temporary vector for the sums
    double* lnSum = (double*)malloc(numStates*sizeof(double));
    double* lnSumEnd = lnSum + numStates;
    
    // calculate the conditional likelihoods
    int activeTi = tiMngr->getActiveTiProb();
    std::vector<Node*>& dpSeq = tree->getInteriorDownPassSequence();
    double lnFactor = 0.0;
    for (auto p=dpSeq.begin(); p != dpSeq.end(); p++)
        {
        for (double* v=lnSum; v != lnSumEnd; v++)
            (*v) = 0.0;
            
        std::set<Node*>& pDesc = (*p)->getDescendants();
        auto tasks = new ConditionalLikelihoodTask[pDesc.size()];
        auto task = tasks;
        for (Node* d : pDesc)
            {
            double* P = d->getTransitionProbabilityPair()->tipr[activeTi]->begin();
            double* clBegin = d->getConditionalLikelihood();
            double* bottomClBegin = d->getBottomConditionalLikelihood();

            task->init((int)numStates, clBegin, bottomClBegin, P);
            threadPool->PushTask(task);
            ++task;
            }
            
        threadPool->Wait();
        delete [] tasks;
        
        // calculate log sums
        for (Node* d : pDesc)
            {
            double* s = d->getBottomConditionalLikelihood();
            for (double* v = lnSum; v != lnSumEnd; v++, s++)
                *v += *s;
            }
            
        // rescale
        double lnLargestCl = lnSum[0];
        for (double* v=&lnSum[1]; v != lnSumEnd; v++)
            {
            if (*v > lnLargestCl)
                lnLargestCl = *v;
            }
        lnFactor += lnLargestCl;
        for (double* v=lnSum; v != lnSumEnd; v++)
            *v -= lnLargestCl;

        double* LP = (*p)->getConditionalLikelihood();
        for (double* v=lnSum; v != lnSumEnd; v++)
            {
            *LP = exp(*v);
            LP++;
            }

#       if 0
        for (int i=0; i<numStates; i++)
            std::cout << "cl[" << i << "] = " << std::fixed << std::setprecision(50) << LP[i] << std::endl;
        getchar();
#       endif
        }
        
    // calculate likelihood
    double* LR = tree->getRoot()->getConditionalLikelihood();
    double* f = &pi[activePi][0];
    double like = 0.0;
    for (int i=0; i<numStates; i++)
        {
        like += (*f) * (*LR);
        LR++;
        f++;
        }
        
    free(lnSum);

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "s" << std::endl;

    return log(like) + lnFactor;
}

#endif

double Model::lnPriorProbability(void) {

    return 0.0;
}
