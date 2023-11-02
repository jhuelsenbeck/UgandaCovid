#include <chrono>
#include <iomanip>
#include <arm_neon.h>
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "TransitionProbabilities.hpp"
#include "TransitionProbabilitiesMngr.hpp"
#include "Tree.hpp"
#undef INTEL_SSE // define/undef for using or not using SSE
#ifdef INTEL_SSE
#   include <xmmintrin.h>
#   include <emmintrin.h>
#   include <pmmintrin.h>
#endif

#define SHOW_LNL_TIME



Model::Model(Tree* tp, RateMatrix* m, ThreadPool* thp) {

    threadPool = thp;
    numStates = m->getNumRows();
    updateType = "";
    
    // check that we can use the SSE instruction set
#   ifdef INTEL_SSE
    if (numStates % 2 != 0)
        Msg::error("Trying to use the SSE instruction set for an odd number of states");
#   endif
    std::cout << "   * Model has " << numStates << " areas (" << numStates << " X " << numStates << ")" << std::endl;

    // set up the model parameters
    substitutionRate[0] = 0.01; // pick a good starting value?
    substitutionRate[1] = substitutionRate[0];
    pi[0].resize(numStates);
    pi[1].resize(numStates);
    std::fill(pi[0].begin(), pi[0].end(), 1.0/numStates);
    pi[1] = pi[0];
    r[0].resize(numStates*(numStates-1)/2);
    std::fill(r[0].begin(), r[0].end(), 1.0);
    double sum = 0.0;
    for (int i=0; i<r[0].size(); i++)
        sum += r[0][i];
    for (int i=0; i<r[0].size(); i++)
        r[0][i] /= sum;
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
    tiMngr->updateTransitionProbabilities(substitutionRate[0]);
        
    // set up the conditional likelihoods (must be done after tiMngr is instantiated
    initializeConditionalLikelihoods();
}

Model::~Model(void) {

    delete tree;
    delete q[0];
    delete q[1];
    delete tiMngr;
    delete [] condLikes;
}

void Model::accept(void) {

    // adjust for acceptance
    substitutionRate[0] = substitutionRate[1];
    pi[0] = pi[1];
    r[0] = r[1];
}

void Model::checkConditionalLikelihoods(void) {

    std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    for (int i=0, n=(int)dpSeq.size(); i<n; i++)
        {
        Node* p = dpSeq[i];
        double* m = p->getConditionalLikelihood();
        double* mEnd = p->getConditionalLikelihoodEnd();
        double sum = 0.0;
        for (double* d = m; d<mEnd; d++)
            sum += *d;
        if (sum < 0.0 || sum > numStates)
            Msg::warning("Problem in conditional likelihoods");
        if (p->getIsTip() == true && sum == 0.0)
            Msg::warning("Did not initialize tip correctly");
        }
}

void Model::initializeConditionalLikelihoods(void) {
        
    size_t numAllocatedDoubles = 0;
    size_t numNodes = tree->getNumNodes();
    size_t memSize = numNodes * numStates * sizeof(double);
    
    // allocate the conditional likelihoods
    condLikes = (double*)malloc(memSize);
    if (!condLikes)
        Msg::error("Failure to allocate conditional likelihood for tree");
    memset(condLikes, 0, memSize);
    
    double* m = condLikes;
    std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    for (auto p=dpSeq.begin(); p != dpSeq.end(); p++)
        {
        // initialize the conditional likelihoods
        if ((*p)->getIsTip() == true)
            {
            int areaId = (*p)->getAreaId();
            if (areaId != -1)
                {
                m[areaId] = 1.0;
                }
            else
                {
                for (int i=0; i<numStates; i++)
                    m[i] = 1.0;
                }
            }
        (*p)->setConditionalLikelihood(m);
        (*p)->setConditionalLikelihoodEnd(m+numStates);
        m += numStates;
        numAllocatedDoubles += numStates;
        
        // set up the transition probability for the node
        TransitionProbabilities* tip = tiMngr->getTiProb((*p)->getBrlen());
        if (tip == nullptr)
            {
            std::cout << "Branch Length = " << (*p)->getBrlen() << std::endl;
            tiMngr->printMap();
            Msg::error("Could not find transition probabilities for node " + std::to_string((*p)->getIndex()));
            }
        (*p)->setTransitionProbability(tip);
        }
    if (m != condLikes + (numNodes*numStates))
        Msg::warning("Something unexpected when setting up conditional likelihoods");
    std::cout << "Successfully allocated a total of " << numAllocatedDoubles << " doubles for the conditional likelihoods" << std::endl;

    double lnL = lnLikelihood();
    std::cout << "lnL = " << lnL << std::endl;
    
#   if 0
    checkConditionalLikelihoods();
#   endif
}

#ifdef INTEL_SSE

// vector version of likelihood calculation
double Model::lnLikelihood(void) {

#   if defined(SHOW_LNL_TIME)
    auto begin = std::chrono::high_resolution_clock::now();
#   endif

    // allocate a temporary vector for the sums
    double* clSum = (double*)malloc(numStates*sizeof(double));
    double* clSumEnd = clSum + numStates;
    
    // calculate the conditional likelihoods
    std::vector<Node*>& dpSeq = tree->getInteriorDownPassSequence();
    double lnFactor = 0.0;
    double lnNodeScalar = 0.0;
    for (auto p=dpSeq.begin(); p != dpSeq.end(); p++)
        {
        for (double* v=clSum; v != clSumEnd; v++)
            (*v) = 1.0;
        std::set<Node*>& pDesc = (*p)->getDescendants();
        for (Node* d : pDesc)
            {
            double* P = d->getTransitionProbability()->begin();
            double* clBegin = d->getConditionalLikelihood();
            double* clEnd = d->getConditionalLikelihoodEnd();
            double largestCl = 0.0;
            for (double* s=clSum; s != clSumEnd; s++)
                {
                __m128d sum01;
                double sum = 0.0;
                for (double* L=clBegin; L != clEnd; L+=2)
                    {
                    __m128d cl_01 = _mm_load_pd(L);

                    __m128d ti_01 = _mm_load_pd(P);

                    __m128d p01 = _mm_mul_pd(cl_01,ti_01);

                    sum01 = _mm_hadd_pd(sum01,p01);

                    P += 2;
                    }
                _mm_store_pd(&sum,sum01);
                *s *= sum;
                if (*s > largestCl)
                    largestCl = *s;
                
                
//                double sum = 0.0;
//                for (double* L=clBegin; L != clEnd; L++)
//                    {
//                    sum += (*P) * (*L);
//                    P++;
//                    }
//                *s *= sum;
//
//                if (*s > largestCl)
//                    largestCl = *s;
                }
               
            double factor = 1.0 / largestCl;
            lnNodeScalar += log(largestCl); // is this right?
            for (double* s=clSum; s != clSumEnd; s++)
                *s *= factor;
            }
        lnFactor += lnNodeScalar;

        double* LP = (*p)->getConditionalLikelihood();
        for (double* v=clSum; v != clSumEnd; v++)
            {
            *LP = *v;
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
    double* f = &pi[1][0];
    double like = 0.0;
    for (int i=0; i<numStates; i++)
        {
        like += (*f) * (*LR);
        LR++;
        f++;
        }
        
    free(clSum);

#   if defined(SHOW_LNL_TIME)
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "s" << std::endl;
#   endif

    return log(like) + lnFactor;
}


#else

// serial version of likelihood calculation
double Model::lnLikelihood(void) {

#   if defined(SHOW_LNL_TIME)
    auto begin = std::chrono::high_resolution_clock::now();
#   endif

    // allocate a temporary vector for the sums
    double* clSum = (double*)malloc(numStates*sizeof(double));
    double* clSumEnd = clSum + numStates;
    
    // calculate the conditional likelihoods
    std::vector<Node*>& dpSeq = tree->getInteriorDownPassSequence();
    double lnFactor = 0.0;
    double lnNodeScalar = 0.0;
    for (auto p=dpSeq.begin(); p != dpSeq.end(); p++)
        {
        for (double* v=clSum; v != clSumEnd; v++)
            (*v) = 1.0;
        std::set<Node*>& pDesc = (*p)->getDescendants();
        for (Node* d : pDesc)
            {
            double* P = d->getTransitionProbability()->begin();
            double* clBegin = d->getConditionalLikelihood();
            double* clEnd = d->getConditionalLikelihoodEnd();
            double largestCl = 0.0;
            for (double* s=clSum; s != clSumEnd; s++)
                {
                double sum = 0.0;
                for (double* L=clBegin; L != clEnd; L++)
                    {
                    sum += (*P) * (*L);
                    P++;
                    }
                *s *= sum;
                
                if (*s > largestCl)
                    largestCl = *s;
                }
               
            double factor = 1.0 / largestCl;
            lnNodeScalar += log(largestCl); // is this right?
            for (double* s=clSum; s != clSumEnd; s++)
                *s *= factor;
            }
        lnFactor += lnNodeScalar;

        double* LP = (*p)->getConditionalLikelihood();
        for (double* v=clSum; v != clSumEnd; v++)
            {
            *LP = *v;
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
    double* f = &pi[1][0];
    double like = 0.0;
    for (int i=0; i<numStates; i++)
        {
        like += (*f) * (*LR);
        LR++;
        f++;
        }
        
    free(clSum);

#   if defined(SHOW_LNL_TIME)
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "s" << std::endl;
#   endif

    return log(like) + lnFactor;
}

#endif

double Model::lnPriorProbability(void) {

    double lnPriorProb = 0.0;
    
    // exponential prior on the substitution rate
    double expParm = 100.0;
    lnPriorProb += Probability::Exponential::lnPdf(expParm, substitutionRate[1]);
    
    // flat Dirichlet distribution on pi
    lnPriorProb += Probability::Helper::lnGamma(numStates);

    // flat Dirichlet distribution on r
    lnPriorProb += Probability::Helper::lnGamma(r[0].size());
    
    return lnPriorProb;
}

void Model::reject(void) {

    // adjust for rejection
    substitutionRate[1] = substitutionRate[0];
    pi[1] = pi[0];
    r[1] = r[0];
    switchActiveRateMatrix();
}

void Model::switchActiveRateMatrix(void) {

    if (activeRateMatrix == 0)
        activeRateMatrix = 1;
    else
        activeRateMatrix = 0;
}

double Model::update(void) {

    RandomVariable& rng = RandomVariable::getInstance();
    double u = rng.uniformRv();
    
    double lnProposalProb = 0.0;
    if (u <= 0.33)
        {
        lnProposalProb = updateSubstitutionRate();
        }
    else if (u > 0.33 && u <= 0.67)
        {
        lnProposalProb = updatePi();
        updateRateMatrix();
        }
    else
        {
        lnProposalProb = updateR();
        updateRateMatrix();
        }
        
    // update the transition probabilities
    tiMngr->updateTransitionProbabilities(substitutionRate[1]);
        
    return lnProposalProb;
}

void Model::updateRateMatrix(void) {

    switchActiveRateMatrix();
    q[activeRateMatrix]->set(&pi[1][0], &r[1][0]);
}

double Model::updateSubstitutionRate(void) {

    updateType = "substitution rate";
    RandomVariable& rng = RandomVariable::getInstance();
    double oldValue = substitutionRate[0];
    double tuning = log(4.0);
    double factor = exp((rng.uniformRv()-0.5)*tuning);
    double newValue = oldValue * factor;
    substitutionRate[1] = newValue;
    return log(newValue) - log(oldValue);
}

double Model::updatePi(void) {

    updateType = "equilibrium frequencies";
    return updateSimplex(pi[0], pi[1], 10000.0, 10);
}

double Model::updateR(void) {

    updateType = "exchangeability rates";
    return updateSimplex(r[0], r[1], 10000.0, 10);
}

double Model::updateSimplex(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0) {

    RandomVariable& rng = RandomVariable::getInstance();
    
    std::vector<double> alphaForward(r[0].size());
    for (int i=0, n=(int)alphaForward.size(); i<n; i++)
        alphaForward[i] = oldVec[i] * alpha0;
        
    Probability::Dirichlet::rv(&rng, alphaForward, newVec);
    Probability::Helper::normalize(newVec, 10e-7);
    
    std::vector<double> alphaReverse(r[0].size());
    for (int i=0, n=(int)alphaReverse.size(); i<n; i++)
        alphaReverse[i] = newVec[i] * alpha0;

    // CHECK THIS!!!
    double lnForwardProb = Probability::Dirichlet::lnPdf(alphaReverse, oldVec) - Probability::Dirichlet::lnPdf(alphaForward, newVec);
    return lnForwardProb;
}

double Model::updateSimplex(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, int k) {

    RandomVariable& rng = RandomVariable::getInstance();
    
    // choose which elements to update
    int n = (int)oldVec.size();
    if (n == k)
        Msg::error("Cannot update that many simplex variables!!!");
        
    std::set<int> indices;
    while (indices.size() < k)
        {
        int idx = (int)(rng.uniformRv() * n);
        indices.insert(idx);
        }
    
    // fill in old values
    std::vector<double> oldValues(k+1);
    double sum = 0.0;
    int i = 0;
    for (int idx : indices)
        {
        double x = oldVec[idx];
        oldValues[i++] = x;
        sum += x;
        }
    oldValues[i] = 1.0 - sum;
    
    std::vector<double> alphaForward(k + 1);
    for (int i=0, n=(int)alphaForward.size(); i<n; i++)
        alphaForward[i] = oldValues[i] * alpha0;
        
    std::vector<double> newValues(k+1);
    Probability::Dirichlet::rv(&rng, alphaForward, newValues);
    Probability::Helper::normalize(newVec, 10e-7);
    
    std::vector<double> alphaReverse(k + 1);
    for (int i=0, n=(int)alphaReverse.size(); i<n; i++)
        alphaReverse[i] = newValues[i] * alpha0;
        
    // fill in the vector from the updated (newValues) vector
    double factor = newValues[k] / oldValues[k];
    for (int i=0, n=(int)newVec.size(); i<n; i++)
        newVec[i] *= factor;
    i = 0;
    for (int idx : indices)
        newVec[idx] = newValues[i++];

    // CHECK THIS!!!
    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldVec) - Probability::Dirichlet::lnPdf(alphaForward, newVec);
    lnProposalProb += (n - k - 1) * log(factor);
    return lnProposalProb;
}
