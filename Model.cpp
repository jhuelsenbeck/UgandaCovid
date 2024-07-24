#include <chrono>
#include <iomanip>
#include "CondLikeJobMngr.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "TransitionProbabilities.hpp"
#include "TransitionProbabilitiesMngr.hpp"
#include "Tree.hpp"

#undef SHOW_LNL_TIME



Model::Model(Tree* tp, RateMatrix* m, ThreadPool* thp, CondLikeJobMngr* mngr) {

    clManager = mngr;
    threadPool = thp;
    numStates = m->getNumRows();
    updateType = "";
    
    // check that we can use the SSE instruction set
    std::cout << "     Model has " << numStates << " areas (" << numStates << " X " << numStates << ")" << std::endl;

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

    //clMngr.calculateTime(10);
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
    std::cout << "     Successfully allocated a total of " << numAllocatedDoubles << " doubles for the conditional likelihoods" << std::endl;

    lnLikelihood();
    
#   if 1
    checkConditionalLikelihoods();
#   endif
}

#if 1

double Model::lnLikelihood(void) {

#   if defined(SHOW_LNL_TIME)
    auto begin = std::chrono::high_resolution_clock::now();
#   endif

    // calculate conditional likelihoods (threaded)
    clManager->calculateConditionalLikelihoods();
    double lnFactor = clManager->getScaler();
    
    // average likelihoods at root
    Node* r = tree->getRoot();
    double* cl = r->getConditionalLikelihood();
    double* f = &pi[1][0];
    double like = 0.0;
    for (int i=0; i<numStates; i++)
        {
        like += (*f) * (*cl);
        cl++;
        f++;
        }

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

        lnFactor += largestLogProb;
        double* clP = (*p)->getConditionalLikelihood();
        for (double* s=clSum; s != clSumEnd; s++)
            {
            *s -= largestLogProb;
            *clP = exp(*s);
            clP++;
            }
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

void Model::map(void) {

    std::cout << "Begin stochastic mapping" << std::endl;
    
    // allocate a vector for probabilities
    double* probs = new double[numStates];
    double* probsEnd = probs + numStates;
    
    // get a reference to the random number generator
    RandomVariable& rng = RandomVariable::getInstance();
    
    // calculate conditional likelihoods to the root (threaded)
    clManager->calculateConditionalLikelihoods();
    
    // pick a state at the root
    Node* r = tree->getRoot();
    double* cl = r->getConditionalLikelihood();
    double* f = &pi[1][0];
    double sum = 0.0;
    for (double* x=probs; x != probsEnd; x++)
        {
        (*x) = (*f) * (*cl);
        sum += (*x);
        cl++;
        f++;
        }
    double u = rng.uniformRv() * sum;
    sum = 0.0;
    int areaId = 0;
    for (double* x=probs; x != probsEnd; x++)
        {
        sum += (*x);
        if (u < sum)
            {
            r->setAreaId(areaId);
            break;
            }
        areaId++;
        }

    // move up the tree, picking states for each node
    std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    for (int n=(int)dpSeq.size()-1; n>=0; n--)
        {
        Node* p = dpSeq[n];
        if (p->getIsAreaFixed() == false && p != tree->getRoot())
            {
            int pAncState     = p->getAncestor()->getAreaId();
            double* p_ij      = p->getTransitionProbability()->begin();
            p_ij += pAncState * numStates;
            double* clP_begin = p->getConditionalLikelihood();
            double* clP_end   = p->getConditionalLikelihoodEnd();
            
            sum = 0.0;
            for (auto cl=clP_begin, x=probs; cl != clP_end; cl++, x++)
                {
                (*x) = (*cl) * (*p_ij);
                sum += (*x);
                p_ij++;
                }
            
            u = rng.uniformRv() * sum;
            sum = 0.0;
            int areaId = 0;
            for (double* x=probs; x != probsEnd; x++)
                {
                sum += (*x);
                if (u < sum)
                    {
                    p->setAreaId(areaId);
                    break;
                    }
                areaId++;
                }
                
            }
        }

    // simulate along each branch
    int numChanges = 0;
    RateMatrix* rateMatrix = q[activeRateMatrix];
    for (int n=0; n<dpSeq.size(); n++)
        {
        Node* p = dpSeq[n];
        if (p != tree->getRoot())
            {
            int begState = p->getAncestor()->getAreaId();
            int endState = p->getAreaId();
            TransitionProbabilities* p_ij = p->getTransitionProbability();
            int v_int = p_ij->getBrlen();
            double v = (double)v_int * substitutionRate[0];
            if (v_int == 0)
                v = 10e-5;
                
            int curState;
            int numRejections = 0;
            do
                {
                curState = begState;
                int count = 0;
                double pt = 0.0;
                while (pt < v)
                    {
                    double rate = -(*rateMatrix)(curState,curState);
                    if (count == 0 && begState != endState)
                        pt += -log(1.0 - rng.uniformRv()*(1.0 - exp(-rate*v))) / rate;
                    else
                        pt += -log(rng.uniformRv()) / rate;
                    if (pt < v)
                        {
                        double u = rng.uniformRv() * rate;
                        sum = 0.0;
                        for (int j=0; j<numStates; j++)
                            {
                            if (j != curState)
                                {
                                sum += (*rateMatrix)(curState, j);
                                if (u < sum)
                                    {
                                    curState = j;
                                    break;
                                    }
                                }
                            }
                        count++;
                        }
                    }
                    
                if (curState != endState)
                    numRejections++;
                else
                    numChanges += count;
                if (numRejections > 1000)
                    Msg::warning("Large number of rejections!");
                } while(curState != endState);
            
            std::cout << p->getIndex() << " -- " << begState << " -> " << endState << std::endl;
            if (p->getAreaId() == -1 || p->getAncestor()->getAreaId() == -1)
                exit(1);
            }
        }
    delete [] probs;
    
    std::cout << "Number of area changes = " << numChanges << std::endl;
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
    double tuning = 0.02;
    double factor = exp((rng.uniformRv()-0.5)*tuning);
    double newValue = oldValue * factor;
    substitutionRate[1] = newValue;
    return log(newValue) - log(oldValue);
}

double Model::updatePi(void) {

    updateType = "equilibrium frequencies";
    return updateSimplex(pi[0], pi[1], 5000.0, 1);
}

double Model::updateR(void) {

    updateType = "exchangeability rates";
    return updateSimplex(r[0], r[1], 5000.0, 1);
}

double Model::updateSimplex(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0) {

    RandomVariable& rng = RandomVariable::getInstance();
    
    std::vector<double> alphaForward(r[0].size());
    for (int i=0, n=(int)alphaForward.size(); i<n; i++)
        alphaForward[i] = oldVec[i] * alpha0 + 1.0;
        
    Probability::Dirichlet::rv(&rng, alphaForward, newVec);
    Probability::Helper::normalize(newVec, 10e-7);
    
    std::vector<double> alphaReverse(r[0].size());
    for (int i=0, n=(int)alphaReverse.size(); i<n; i++)
        alphaReverse[i] = newVec[i] * alpha0 + 1.0;

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
        alphaForward[i] = oldValues[i] * alpha0 + 1.0;
        
    std::vector<double> newValues(k+1);
    Probability::Dirichlet::rv(&rng, alphaForward, newValues);
    Probability::Helper::normalize(newVec, 10e-7);
    
    std::vector<double> alphaReverse(k + 1);
    for (int i=0, n=(int)alphaReverse.size(); i<n; i++)
        alphaReverse[i] = newValues[i] * alpha0 + 1.0;
        
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
