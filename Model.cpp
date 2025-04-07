#include <chrono>
#include <fstream>
#include <iomanip>
#include "BitVector.hpp"
#include "CondLikeJobMngr.hpp"
#include "History.hpp"
#include "MetaData.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "Samples.hpp"
#include "TransitionProbabilities.hpp"
#include "TransitionProbabilitiesMngr.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#define MIN_BRLEN 10e-5
#undef SHOW_LNL_TIME



Model::Model(Tree* tp, MetaData* md, ThreadPool* thp, CondLikeJobMngr* mngr) {

    // initialize the rate matrix
    RateMatrix* Q = new RateMatrix(md->getAreas());

    clManager = mngr;
    threadPool = thp;
    numStates = Q->getNumRows();
    updateType = "";
    
    // check that we can use the SSE instruction set
    std::cout << "     Model has " << numStates << " areas (" << numStates << " X " << numStates << ")" << std::endl;

    // set up the model parameters
    initializeParameters(tp, Q);
        
    // set up histories for each node
    initializeHistories();
    
    // set up the conditional likelihoods (must be done after tiMngr is instantiated
    initializeConditionalLikelihoods();
}

Model::~Model(void) {

    deleteHistories();
    delete tree;
    delete q[0];
    delete q[1];
    delete uniformizedRateMatrix;
    delete tiMngr;
    delete [] condLikes;
    for (int i=0; i<matrixPowers.size(); i++)
        delete matrixPowers[i];
}

void Model::accept(void) {

    // adjust for acceptance
    substitutionRate[0] = substitutionRate[1];
    pi[0] = pi[1];
    r[0] = r[1];
}

void Model::assignNodeStates(RandomVariable* rng) {

    // allocate a vector for probabilities
    double* probs = new double[numStates];
    double* probsEnd = probs + numStates;

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
    double u = rng->uniformRv() * sum;
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
            
            u = rng->uniformRv() * sum;
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

    delete [] probs;
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

void Model::deleteHistories(void) {

    std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    for (auto p=dpSeq.begin(); p != dpSeq.end(); p++)
        {
        History* h = (*p)->getHistory();
        delete h;
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

void Model::initializeHistories(void) {

    std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    for (auto p=dpSeq.begin(); p != dpSeq.end(); p++)
        {
        History* h = new History;
        (*p)->setHistory(h);
        }
}

void Model::initializeMatrixPowers(int num) {

    if (num < 2)
        num = 2;
    if (matrixPowers.size() < num)
        {
        for (size_t i=0, n=num-matrixPowers.size(); i<n; i++)
            matrixPowers.push_back(new RateMatrix(*uniformizedRateMatrix));
        }
    
    matrixPowers[0]->setIdentity();
    (*matrixPowers[1]) = (*uniformizedRateMatrix);
    for (size_t i=2; i<num; i++)
        matrixPowers[i]->multiply(*uniformizedRateMatrix, *matrixPowers[i-1]);
}

void Model::initializeParameters(Tree* tp, RateMatrix* m) {
    
    // read a file with initial parameter values
    bool foundStartingValues = false;
    UserSettings& settings = UserSettings::getUserSettings();
    std::vector<double> samplePi;
    std::vector<double> sampleR;
    double sampleRate = -1.0;
    if (settings.getInitialParameterValues() != "")
        {
        std::vector<Samples*> samples = readParameterFile(settings.getInitialParameterValues());
        for (Samples* s : samples)
            {
            if (s->getParmName() == "Rate")
                sampleRate = s->lastSample();
            else if (s->getParmName() == "R")
                sampleR.push_back(s->lastSample());
            else if (s->getParmName() == "Pi")
                samplePi.push_back(s->lastSample());
            }

        foundStartingValues = true;

        if (samplePi.size() != numStates)
            foundStartingValues = false;
        if (sampleRate < 0.0)
            foundStartingValues = false;
        if (sampleR.size() != numStates * (numStates-1) / 2)
            foundStartingValues = false;

        double sum = 0.0;
        for (int i=0; i<samplePi.size(); i++)
            sum += samplePi[i];
        for (int i=0; i<samplePi.size(); i++)
            samplePi[i] /= sum;
        
        sum = 0.0;
        for (int i=0; i<sampleR.size(); i++)
            sum += sampleR[i];
        for (int i=0; i<sampleR.size(); i++)
            sampleR[i] /= sum;
        }
    
    // set up subsitution rate
    if (foundStartingValues == true)
        substitutionRate[0] = sampleRate;
    else
        substitutionRate[0] = 0.01;
    substitutionRate[1] = substitutionRate[0];
    
    // set up stationary frequencies
    pi[0].resize(numStates);
    pi[1].resize(numStates);
    if (foundStartingValues == true)
        pi[0] = samplePi;
    else
        std::fill(pi[0].begin(), pi[0].end(), 1.0/numStates);
    pi[1] = pi[0];
    
    // set up substitution rate
    r[0].resize(numStates*(numStates-1)/2);
    if (foundStartingValues == true)
        r[0] = sampleR;
    else
        std::fill(r[0].begin(), r[0].end(), 1.0);
    double sum = 0.0;
    for (int i=0; i<r[0].size(); i++)
        sum += r[0][i];
    for (int i=0; i<r[0].size(); i++)
        r[0][i] /= sum;
    r[1] = r[0];
    
    // set colonization rate

    // set up the rate matrix and tree
    tree = tp;
    activeRateMatrix = 0;
    q[0] = m;
    q[1] = new RateMatrix(*m);
    q[0]->set(&pi[0][0], &r[0][0]);
    (*q[1]) = (*q[0]);
    //std::cout << *q[0] << std::endl;
    //std::cout << *q[1] << std::endl;
    uniformizedRateMatrix = new RateMatrix(*m);
    
    // set up the transition probabilities
    tiMngr = new TransitionProbabilitiesMngr(this, tree, numStates, threadPool);
    tiMngr->updateTransitionProbabilities(substitutionRate[0]);
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
        
    // get a reference to the random number generator
    RandomVariable& rng = RandomVariable::getInstance();
    
    //tiMngr->updateTransitionProbabilities(getSubstitutionRate());
    
    // calculate conditional likelihoods to the root (threaded)
    clManager->calculateConditionalLikelihoods();
    
    // assign states to each end of every branch
    assignNodeStates(&rng);

    // sample history for each branch
    //sampleHistoriesUsingRejectionSamplign(&rng);
    int numChanges = sampleHistoriesUsingUniformization(&rng);
    std::cout << "           Number of changes = " << numChanges << std::endl;
}

int Model::parsimonyScore(void) {

    std::vector<BitVector*> stateSets(tree->getNumNodes());
    for (int i=0; i<stateSets.size(); i++)
        stateSets[i] = new BitVector(numStates);
    std::vector<int> descendantStateCounts(numStates);

    std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    int n = 0;
    for (Node* p : dpSeq)
        {
        if (p->getIsTip() == true)
            {
            int areaId = p->getAreaId();
            if (areaId == -1)
                stateSets[p->getIndex()]->set();
            else
                stateSets[p->getIndex()]->set(p->getAreaId());
            }
        else
            {
            std::set<Node*>& pDesc = p->getDescendants();
            BitVector* pStateSet = stateSets[p->getIndex()];
            pStateSet->set();
            for (size_t i=0; i<numStates; i++)
                descendantStateCounts[i] = 0;
            for (Node* d : pDesc)
                {
                BitVector* dStateSet = stateSets[d->getIndex()];
                *pStateSet &= *dStateSet;
                int nSet = 0;
                for (size_t i=0; i<numStates; i++)
                    {
                    bool tf = (*dStateSet)[i];
                    if (tf == true)
                        {
                        nSet++;
                        descendantStateCounts[i]++;
                        }
                    }
                if (nSet != 0)
                    {
                    
                    }
//As for calculating parsimony scores on nonbinary trees… in PAUP* I haven’t done too much optimization because my tree searches are always on binary trees and calculation of tree lengths for user-specified trees is fast enough that I never worried much about it. Basically, I just use a simple application of the usual dynamic programming approach where the state set assigned to each node in a postorder traversal is the intersection of the state sets of the children, if it is nonempty. If this intersection is empty, the new state set is the union of the most frequent states in the state sets of the children. The tree length is increased by one each time the union is required. It’s essentially the algorithm of Hartigan 1973 (attached). I use some bit tricks to parallelize over sites (essentially what is described in the attached paper by Goloboff—he claimed his method was faster that PAUP’s and if that’s true it’s probably because he’s branching on special cases because otherwise the algorithms are essentially the same).

                }
            }
        }

    for (int i=0; i<stateSets.size(); i++)
        delete stateSets[i];
        
    return n;
}

void Model::printMatrixPowers(void) {

    for (int i=0; i<matrixPowers.size(); i++)
        {
        std::cout << "M^{" << i << "}" << std::endl;
        std::cout << (*matrixPowers[i]) << std::endl;
        }
}

std::vector<Samples*> Model::readParameterFile(std::string fn) {

    std::vector<Samples*> parmSamples;
    
    std::ifstream parmStrm( fn.c_str(), std::ios::in );
    if (!parmStrm)
        Msg::error("Cannot open file \"" + fn + "\"");
    
    std::string lineString = "";
    int lineNum = 0;
    while( getline(parmStrm, lineString).good() )
        {
        std::istringstream linestream(lineString);
        int ch;
        std::string word = "";
        int wordNum = 0;
        std::string cmdString = "";
        do
            {
            word = "";
            linestream >> word;
            wordNum++;
            if (lineNum == 0)
                {
                /* read the number of taxa/chars from the first line */
                parmSamples.push_back( new Samples(word) );
                }
            else
                {
                double x;
                std::istringstream buf(word);
                buf >> x;
                parmSamples[wordNum-1]->addSample(x);
                }
            } while ( (ch=linestream.get()) != EOF );
            
        lineNum++;
        }

    //for (Samples* s : parmSamples)
    //    s->print();
    
    parmStrm.close();
    
    return parmSamples;
}

void Model::reject(void) {

    // adjust for rejection
    substitutionRate[1] = substitutionRate[0];
    pi[1] = pi[0];
    r[1] = r[0];
    switchActiveRateMatrix();
}

int Model::sampleHistoriesUsingRejectionSamplign(RandomVariable* rng) {

    int numChanges = 0;
    RateMatrix* rateMatrix = q[activeRateMatrix];
    std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    for (int n=0; n<dpSeq.size(); n++)
        {
        Node* p = dpSeq[n];
        if (p != tree->getRoot())
            {
            History* h = p->getHistory();
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
                h->clearHistory();
                curState = begState;
                int count = 0;
                double pt = 0.0;
                while (pt < v)
                    {
                    double rate = -(*rateMatrix)(curState,curState);
                    if (count == 0 && begState != endState)
                        pt += -log(1.0 - rng->uniformRv()*(1.0 - exp(-rate*v))) / rate;
                    else
                        pt += -log(rng->uniformRv()) / rate;
                    if (pt < v)
                        {
                        double u = rng->uniformRv() * rate;
                        double sum = 0.0;
                        for (int j=0; j<numStates; j++)
                            {
                            if (j != curState)
                                {
                                sum += (*rateMatrix)(curState, j);
                                if (u < sum)
                                    {
                                    h->addChange(curState, j, pt);
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
                if (numRejections > 10000)
                    Msg::warning("Large number of rejections!");
                } while(curState != endState);
            
            if (p->getAreaId() == -1 || p->getAncestor()->getAreaId() == -1)
                Msg::error("Area assignment at at least one end of branch is ambiguous");
            }
        }
    return numChanges;
}

int Model::sampleHistoriesUsingUniformization(RandomVariable* rng) {

    double nFactorial[20];
    nFactorial[0] = 1.0;
    nFactorial[1] = 1.0;
    for (int i=2; i<20; i++)
        nFactorial[i] = nFactorial[i-1] * i;
    
    updateRateMatrix();
    RateMatrix* rateMatrix = q[activeRateMatrix];
    double mu = rateMatrix->uniformize(uniformizedRateMatrix);
    tiMngr->updateTransitionProbabilities(getSubstitutionRate());

    //std::cout << "R" << std::endl;
    //std::cout << *uniformizedRateMatrix << std::endl;
    initializeMatrixPowers(20);
    //printMatrixPowers();
    
    std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    std::vector<double> changeTimes;
    std::set<Change*> changesToRemove;
    double r = getSubstitutionRate();
    int numChanges = 0;
    for (int nde=0; nde<dpSeq.size(); nde++)
        {
        Node* p = dpSeq[nde];
        if (p != tree->getRoot())
            {
            // get information on branch
            TransitionProbabilities* P = p->getTransitionProbability();
            double v = p->getBrlen() * r;
            if (v < MIN_BRLEN)
                v = MIN_BRLEN;
            
            double rate = mu * v;
            History* h = p->getHistory();
            int a = p->getAncestor()->getAreaId();
            int b = p->getAreaId();
            
            // sample number of events
            double g = (*P)(a,b) * rng->uniformRv();
            int n = 0;
            double sum = 0.0;
            for (int i=0; i<20; i++)
                {
                double nProb = (*matrixPowers[i])(a,b) * pow(rate,(double)i) * exp(-rate) / nFactorial[i];
                sum += nProb;
                if (g < sum)
                    {
                    n = i;
                    break;
                    }
                }
                
            // sample the series of events
            std::vector<int> intermediateStates(n+1);
            intermediateStates[0] = a;
            for (int i=1; i<intermediateStates.size(); i++)
                {
                double sum = 0.0;
                for (int j=0; j<numStates; j++)
                    sum += (*matrixPowers[1])(intermediateStates[i-1],j) * (*matrixPowers[n-i])(j,b);
                double u = rng->uniformRv() * sum;
                sum = 0.0;
                for (int j=0; j<numStates; j++)
                    {
                    sum += (*matrixPowers[1])(intermediateStates[i-1],j) * (*matrixPowers[n-i])(j,b);
                    if (u < sum)
                        {
                        intermediateStates[i] = j;
                        break;
                        }
                    }
                }

            std::vector<double> intermediateTimes;
            for (int i=0; i<n; i++)
                intermediateTimes.push_back(rng->uniformRv()*v);
            intermediateTimes.push_back(0.0);
            sort(intermediateTimes.begin(), intermediateTimes.end());
                        
            h->clearHistory();
            for (int i=1; i<intermediateStates.size(); i++)
                {
                if (intermediateStates[i-1] != intermediateStates[i])
                    {
                    h->addChange(intermediateStates[i-1], intermediateStates[i], intermediateTimes[i]);
                    numChanges++;
                    }
                }
                    
            }
        }
    
    return numChanges;
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

    updateType = "colonization rate";
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
