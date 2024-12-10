#include <chrono>
#include <iostream>
#include <vector>
#include "CondLikeJob.hpp"
#include "CondLikeJobMngr.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "Tree.hpp"



CondLikeJobMngr::CondLikeJobMngr(Tree* t, ThreadPool* tp, int na) : tree(t), threadPool(tp) {
    
    numStates = na;
    int numNodesInJob = 500;
    std::vector<Node*>& downPassSeq = t->getDownPassSequence();
    
    // divide the tree up into jobs with roughly the same number of nodes for each
    CondLikeJob* currentJob = addJob();
    currentJob->setJobId(0);
    int jobNum = 0;
    for (int i=0, n=(int)downPassSeq.size(); i<n; i++)
        {
        Node* p = downPassSeq[i];
        
        // sum the number of nodes in the subtree
        if (p->getIsTip() == true)
            {
            p->setScratchInt(1);
            }
        else
            {
            std::set<Node*>& pDesc = p->getDescendants();
            int sum = 0;
            for (Node* d : pDesc)
                sum += d->getScratchInt();
            p->setScratchInt(sum + 1);
            }
            
        if (p->getScratchInt() != 1 || p == t->getRoot())
            currentJob->addNode(p);
            
        if (p->getScratchInt() > numNodesInJob || p == t->getRoot())
            {
            jobNum++;
            p->setJob(currentJob);
            p->setScratchInt(1);
            if (p != t->getRoot())
                {
                currentJob = addJob();
                currentJob->setJobId(jobNum);
                }
            }
        }
        
    // set the jobs for all nodes
    for (int i=0, n=(int)downPassSeq.size(); i<n; i++)
        {
        Node* p = downPassSeq[i];
        p->setJob(nullptr);
        p->setDependentJob(nullptr);
        }
    for (std::vector<CondLikeJob*>::iterator it = jobs.begin(); it != jobs.end(); it++)
        {
        CondLikeJob* jobPtr = *it;
        std::vector<Node*>* jobNodes = (*it)->getJobNodes();
        for (auto p=jobNodes->begin(); p != jobNodes->end(); p++)
            {
            if ((*p)->getJob() != nullptr)
                Msg::error("Node is already assigned to a job");
            (*p)->setJob(jobPtr);
            }
        }
        
    // resolve the job dependencies
    for (std::vector<CondLikeJob*>::iterator it = jobs.begin(); it != jobs.end(); it++)
        {
        std::vector<Node*>* jobNodes = (*it)->getJobNodes();
        for (auto p=jobNodes->begin(); p != jobNodes->end(); p++)
            {
            CondLikeJob* pJob = (*p)->getJob();
            if (pJob != *it)
                Msg::error("P should be in this job!");
            std::set<Node*>& pDescendants = (*p)->getDescendants();
            for (Node* d : pDescendants)
                {
                CondLikeJob* dJob = d->getJob();
                if (pJob != dJob && dJob != nullptr)
                    {
                    d->setDependentJob(pJob);
                    pJob->addDependency(dJob);
                    }
                if (d->getIsTip() == false && dJob == nullptr)
                    Msg::error("Descendant job should already be assigned");
                }
            }
        }

    int numZeroDependencyJobs = 0;
    for (std::vector<CondLikeJob*>::iterator it = jobs.begin(); it != jobs.end(); it++)
        {
        if ((*it)->getNumDependencies() == 0)
            numZeroDependencyJobs++;
        }
    std::cout << "   * Initializing conditional likelihood jobs" << std::endl;
    std::cout << "     Divided tree into " << jobs.size() << " regions for parallelization" << std::endl;
    std::cout << "     " << numZeroDependencyJobs << " jobs/regions have no dependencies" << std::endl;
    //print();
}

CondLikeJobMngr::~CondLikeJobMngr(void) {

}

CondLikeJob* CondLikeJobMngr::addJob(void) {

    CondLikeJob* newJob = new CondLikeJob(this, threadPool, numStates);
    jobs.push_back(newJob);
    return newJob;
}

void CondLikeJobMngr::calculateConditionalLikelihoods(void) {

    zeroResolvedDependencies();
    jobId = 0;
    for (std::vector<CondLikeJob*>::iterator it = jobs.begin(); it != jobs.end(); it++)
        {
        if ((*it)->getNumDependencies() == 0)
            threadPool->pushTask(*it);
        }
    threadPool->wait();
}

void CondLikeJobMngr::calculateTime(int n) {

    auto begin = std::chrono::high_resolution_clock::now();
    for (int i=0; i<n; i++)
        {
        calculateConditionalLikelihoods();
        }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "     Test time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " milliseconds" << std::endl;
}

double CondLikeJobMngr::getScaler(void) {

    double lnScaler = 0.0;
    for (std::vector<CondLikeJob*>::iterator it = jobs.begin(); it != jobs.end(); it++)
        lnScaler += (*it)->getLnScaler();
    return lnScaler;
}

void CondLikeJobMngr::jobSizeDistribution(void) {

    std::vector<int> dist(20000, 0);
    int sum = 0;
    int maxSize = 0;
    for (CondLikeJob* j : jobs)
        {
        int n = j->numNodesInJob();
        sum += n;
        if (n > maxSize)
            maxSize = n;
        dist[n]++;
        }
    for (int i=0; i<=maxSize; i++)
        {
        if (dist[i] != 0)
            std::cout << i << " -- " << dist[i] << std::endl;
        else
            std::cout << i << " -- " << std::endl;
        }
    std::cout << "Num job nodes = " << sum << std::endl;
}

void CondLikeJobMngr::print(void) {

    for (std::vector<CondLikeJob*>::iterator it = jobs.begin(); it != jobs.end(); it++)
        {
        (*it)->print();
        }
}

void CondLikeJobMngr::zeroResolvedDependencies(void) {

    for (std::vector<CondLikeJob*>::iterator it = jobs.begin(); it != jobs.end(); it++)
        {
        (*it)->setNumResolvedDependencies(0);
        }
}
