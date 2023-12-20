#include <iostream>
#include <limits>
#include <mutex>
#include "CondLikeJob.hpp"
#include "CondLikeJobMngr.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "TransitionProbabilities.hpp"



CondLikeJob::CondLikeJob(CondLikeJobMngr* m, ThreadPool* tp, int na) : ThreadTask() {

    myManager = m;
    threadPool = tp;
    jobId = 0;
    numStates = na;
    numDependencies = 0;
    nodes = new std::vector<Node*>;
    clSum = (double*)malloc(numStates*sizeof(double));
    clSumEnd = clSum + numStates;
}

CondLikeJob::~CondLikeJob(void) {

    delete nodes;
    delete [] clSum;
}

void CondLikeJob::conditionalLikelihood(void) {

    lnScaler = 0.0;
    for (auto p=nodes->begin(); p != nodes->end(); p++)
        {
        for (double* s=clSum; s != clSumEnd; s++)
            (*s) = 0.0;
            
        std::set<Node*>& pDesc = (*p)->getDescendants();
        for (Node* d : pDesc)
            {
            double* p_ij      = d->getTransitionProbability()->begin();
            double* clD_begin = d->getConditionalLikelihood();
            double* clD_end   = d->getConditionalLikelihoodEnd();
            for (double* s=clSum; s != clSumEnd; s++)
                {
                double sum = 0.0;
                for (double* clD=clD_begin; clD != clD_end; clD++)
                    {
                    sum += (*p_ij) * (*clD);
                    p_ij++;
                    }
                *s += log(sum);
                }
            }

        // rescale and convert back into probabilities
        double largestLogProb = *clSum;
        for (double* s=clSum; s != clSumEnd; s++)
            {
            if (*s > largestLogProb)
                largestLogProb = *s;
            }
        double* clP = (*p)->getConditionalLikelihood();
        for (double* s=clSum; s != clSumEnd; s++)
            {
            *s -= largestLogProb;
            *clP = exp(*s);
            clP++;
            }
        lnScaler += largestLogProb;
                    
        // resolve dependence for p
        CondLikeJob* dependentJob = (*p)->getDependentJob();
        if (dependentJob != nullptr)
            dependentJob->resolveDependency();
        }
}

void CondLikeJob::print(void) {

    std::cout << "Job " << jobId << std::endl;
    std::cout << "   Number of nodes : " << nodes->size() << std::endl;
    std::cout << "   Nodes           : ";
    for (auto p=nodes->begin(); p != nodes->end(); p++)
        std::cout << (*p)->getIndex() << " ";
    std::cout << std::endl;
    std::cout << "   No. Dependencies: " << numDependencies << std::endl;
}

void CondLikeJob::resolveDependency(void) {

    mtx.lock();
    numResolvedDependencies++;
    if (numResolvedDependencies == numDependencies)
        threadPool->PushTask(this);
    mtx.unlock();
}

void CondLikeJob::Run(MathCache&) {

    // compute conditional likelihoods here
    conditionalLikelihood(); // actual work to be accomplished
}
