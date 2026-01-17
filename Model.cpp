#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include "CondLikeJobMngr.hpp"
#include "FileManager.hpp"
#include "History.hpp"
#include "json.hpp"
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



Model::Model(RandomVariable* r, Tree* tp, MetaData* md, ThreadPool* thp, CondLikeJobMngr* mngr) : 
    clManager(mngr), metaData(md), rng(r), threadPool(thp) {

    // initialize the rate matrix
    RateMatrix* Q = new RateMatrix(md->getAreas());

    numStates = Q->getNumRows();
    updateType = "";
    
    ugandaIdx = -1;
    for (size_t i=0; i<numStates; i++)
        {
        if (md->getAreas()[i] == "Uganda")
            {
            ugandaIdx = i;
            break;
            }
        }
    if (ugandaIdx == -1)
        Msg::error("Could not find the Uganda index");
    
    // check that we can use the SSE instruction set
    std::cout << "     Model has " << numStates << " areas (" << numStates << " X " << numStates << ")" << std::endl;

    // set the model type
    UserSettings& settings = UserSettings::getUserSettings();
    switch (settings.getLikelihoodModel())
        {
        case LikelihoodModel::JC69:       modelType = jc69;       break;
        case LikelihoodModel::F81:        modelType = f81;        break;
        case LikelihoodModel::CUSTOM_F81: modelType = custom_f81; break;
        case LikelihoodModel::GTR:        modelType = gtr;        break;
        }
    variableUgandaRate = settings.isUgandaRateVariable();
    if (variableUgandaRate == true && modelType != custom_f81)
        Msg::error("You can only run a variable Uganda rate model with the custom_f81 model");

    // set up the model parameters
    initializeParameters(tp, Q);
        
    // set up histories for each node
    initializeHistories();
    
    // set up the conditional likelihoods (must be done after tiMngr is instantiated
    initializeConditionalLikelihoods();
    
    // precalculate factorials
    nFactorial[0] = 1.0;
    nFactorial[1] = 1.0;
    for (int i=2; i<MAX_NUM_CHANGES; i++)
        nFactorial[i] = nFactorial[i-1] * i;
}

Model::~Model(void) {

    deleteHistories();
    delete tree;
    delete q[0];
    delete q[1];
    delete uniformizedRateMatrix;
    delete tiMngr;
    delete [] condLikes;
    for (size_t i=0; i<matrixPowers.size(); i++)
        delete matrixPowers[i];
    delete [] intervalDwellTimes[0];
    delete [] intervalDwellTimes;
    for (int n=0; n<numIntervals; n++)
        {
        delete [] intervalTransitions[n][0];
        delete [] intervalTransitions[n];
        }
    delete [] intervalTransitions;
}

void Model::accept(void) {

    // adjust for acceptance
    substitutionRate[0] = substitutionRate[1];
    pi[0] = pi[1];
    r[0] = r[1];
    kappa[0] = kappa[1];
    kappaLockdown[0] = kappaLockdown[1];
    kappaNoLockdown[0] = kappaNoLockdown[1];
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

void Model::checkPoint(std::string fileName) {

    nlohmann::json j;
    j["substitutionRate"] = substitutionRate[0];
    j["kappa"] = kappa[0];
    j["pi"] = pi[0];
    j["r"] = r[0];

    std::ofstream out(fileName);
    out << j.dump();
}

void Model::computeDwellTimes(double t, double tAnc, double* t0, double* t1, double* t2) {

    // defensive: allow caller to pass null, but do something sane
    double dummy0 = 0.0, dummy1 = 0.0, dummy2 = 0.0;
    if (t0 == nullptr) t0 = &dummy0;
    if (t1 == nullptr) t1 = &dummy1;
    if (t2 == nullptr) t2 = &dummy2;

    *t0 = 0.0;
    *t1 = 0.0;
    *t2 = 0.0;

    // Make sure we have an increasing segment [a,b]
    double a = t;
    double b = tAnc;
    if (a > b)
        std::swap(a, b);

    // Degenerate branch
    if (b <= a)
        return;

    interval_map& intervals = metaData->getIntervalMap();

    for (const auto& kv : intervals) {

        const int    sInt = kv.first.first;
        const int    eInt = kv.first.second;
        const int    idx  = kv.second;

        // Skip malformed/empty intervals
        if (eInt <= sInt)
            continue;

        const double s = static_cast<double>(sInt);
        const double e = static_cast<double>(eInt);

        // overlap of [a,b] with [s,e] (treating intervals as half-open [s,e))
        const double left  = std::max(a, s);
        const double right = std::min(b, e);
        const double overlap = right - left;

        if (overlap <= 0.0)
            continue;

        if (idx == 0)
            *t0 += overlap;
        else if (idx == 1)
            *t1 += overlap;
        else if (idx == 2)
            *t2 += overlap;
        else {
            // If you ever extend to >3 intervals, you can handle it here.
            // For now, silently ignore unexpected indices.
        }
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
                for (size_t i=0; i<numStates; i++)
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
    
    interval_map& intervalInfo = metaData->getIntervalMap();
    numIntervals = (int)intervalInfo.size();
    intervalDwellTimes = new double*[numIntervals];
    intervalDwellTimes[0] = new double[numIntervals * numStates];
    for (int i=1; i<numIntervals; i++)
        intervalDwellTimes[i] = intervalDwellTimes[i-1] + numStates;
    for (int i=0; i<numIntervals; i++)
        for (size_t j=0; j<numStates; j++)
            intervalDwellTimes[i][j] = 0.0;
    
    intervalTransitions = new int**[numIntervals];
    for (int n=0; n<numIntervals; n++)
        {
        intervalTransitions[n] = new int*[numStates];
        intervalTransitions[n][0] = new int[numStates * numStates];
        for (size_t i=1; i<numStates; i++)
            intervalTransitions[n][i] = intervalTransitions[n][i-1] + numStates;
        for (size_t i=0; i<numStates; i++)
            for (size_t j=0; j<numStates; j++)
                intervalTransitions[n][i][j] = 0;
        }
}

void Model::initializeMatrixPowers(size_t num) {

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
    
    // set colonization rate
    substitutionRate[0] = 0.0001;
    substitutionRate[1] = substitutionRate[0];
    
    // set up stationary frequencies
    pi[0].resize(numStates);
    pi[1].resize(numStates);
    std::fill(pi[0].begin(), pi[0].end(), 1.0/numStates);
    pi[1] = pi[0];
    alphaPi.resize(numStates);
    for (size_t i=0; i<numStates; i++)
        alphaPi[i] = 1.0;
    
    // set up exchangeability parameters
    r[0].resize(numStates*(numStates-1)/2);
    std::fill(r[0].begin(), r[0].end(), 1.0);
    double sum = 0.0;
    for (size_t i=0; i<r[0].size(); i++)
        sum += r[0][i];
    for (size_t i=0; i<r[0].size(); i++)
        r[0][i] /= sum;
    r[1] = r[0];
    alphaR.resize(r[0].size());
    for (size_t i=0; i<r[0].size(); i++)
        alphaR[i] = 1.0;
    
    // set up Uganda bias factor
    kappa[0] = 1.0;
    kappa[1] = kappa[0];
    
    // see if the values should come from the checkpoint file
    if (UserSettings::getUserSettings().readCheckpointFile() == true)
        loadCheckpont();
        
    // set up the rate matrix and tree
    tree = tp;
    activeRateMatrix = 0;
    q[0] = m;
    q[1] = new RateMatrix(*m);
    q[0]->set(modelType, &pi[0][0], &r[0][0], kappa[0]);
    (*q[1]) = (*q[0]);
    //std::cout << *q[0] << std::endl;
    //std::cout << *q[1] << std::endl;
    uniformizedRateMatrix = new RateMatrix(*m);
    
    // set up the transition probabilities
    tiMngr = new TransitionProbabilitiesMngr(this, tree, numStates, threadPool, ugandaIdx);
    tiMngr->updateTransitionProbabilities();
}

bool Model::isValidSimplex(const std::vector<double>& x, double eps) {

    double s = 0.0;
    for (size_t i=0; i<x.size(); i++)
        {
        if (!(x[i] > 0.0))
            return false;
        s += x[i];
        }
    return std::fabs(s - 1.0) < eps;
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
    for (size_t i=0; i<numStates; i++)
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
        // LCRS iteration over children
        for (Node* d = (*p)->getFirstChild(); d != nullptr; d = d->getNextSibling())
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
    double expParm = 1000.0;
    lnPriorProb += Probability::Exponential::lnPdf(expParm, substitutionRate[1]);
    
    // flat Dirichlet distribution on pi
    if (modelType != jc69)
        lnPriorProb += Probability::Helper::lnGamma(numStates);

    // flat Dirichlet distribution on r
    if (modelType == gtr)
        lnPriorProb += Probability::Helper::lnGamma(r[0].size());
    
    // kappa
    if (modelType == custom_f81 && variableUgandaRate == false)
        lnPriorProb += -2.0 * log(1.0 + kappa[0]);
    else if (modelType == custom_f81 && variableUgandaRate == true)
        lnPriorProb += (-2.0 * log(1.0 + kappaLockdown[0])) + (-2.0 * log(1.0 + kappaNoLockdown[0]));
   
    return lnPriorProb;
}

void Model::loadCheckpont(void) {

    std::string parmOutFile = UserSettings::getUserSettings().getOutputFile();
    std::string checkPointFile = FileManager::getParentPath(parmOutFile);
    checkPointFile += "/check_point.json";

    // check if file exists
    if (!FileManager::fileExists(checkPointFile))
        {
        std::cout << "Checkpoint file not found: " << checkPointFile << std::endl;
        return;
        }

    std::ifstream in(checkPointFile);
    if (!in)
        Msg::error("Failed to open checkpoint file: " + checkPointFile);

    nlohmann::json j;
    try
        {
        in >> j;
        }
    catch (const nlohmann::json::parse_error& e)
        {
        Msg::error("Failed to parse checkpoint file: " + std::string(e.what()));
        }

    // extract values with key checking
    if (!j.contains("substitutionRate") || !j.contains("kappa") || !j.contains("pi") || !j.contains("r"))
        Msg::error("Checkpoint file missing required fields");
    
    double jSubstitionRate = j["substitutionRate"].get<double>();
    double jKappa = j["kappa"].get<double>();
    std::vector<double> jPi = j["pi"].get<std::vector<double>>();
    std::vector<double> jR = j["r"].get<std::vector<double>>();
    
    substitutionRate[0] = jSubstitionRate;
    substitutionRate[1] = substitutionRate[0];
    
    kappa[0] = jKappa;
    kappa[1] = kappa[0];
    
    if (jPi.size() == pi[0].size())
        {
        pi[0] = jPi;
        pi[1] = pi[0];
        }

    if (jR.size() == r[0].size())
        {
        r[0] = jR;
        r[1] = r[0];
        }
}

void Model::map(void) {
        
    // calculate conditional likelihoods to the root (threaded)
    clManager->calculateConditionalLikelihoods();
    
    // assign states to each end of every branch
    assignNodeStates(rng);

    // sample history for each branch
    //sampleHistoriesUsingRejectionSamplign(&rng);
    if (modelType != gtr)
        updateRateMatrix();
    int numChanges = sampleHistoriesUsingUniformization(rng);
    std::cout << "           Number of changes = " << numChanges << std::endl;
}

void Model::printMatrixPowers(void) {

    for (size_t i=0; i<matrixPowers.size(); i++)
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
    kappa[1] = kappa[0];
    kappaLockdown[1] = kappaLockdown[0];
    kappaNoLockdown[1] = kappaNoLockdown[0];
    switchActiveRateMatrix();
}

int Model::sampleHistoriesUsingRejectionSamplign(RandomVariable* rng) {

    int numChanges = 0;
    RateMatrix* rateMatrix = q[activeRateMatrix];
    std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    for (size_t n=0; n<dpSeq.size(); n++)
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
                        for (int j=0; j<(int)numStates; j++)
                            {
                            if (j != curState)
                                {
                                sum += (*rateMatrix)(curState, j);
                                if (u < sum)
                                    {
                                    h->addChange(curState, j, pt, 0); // note that you should figure out interval information
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

    // zero out summary information for mapping
    for (int n=0; n<numIntervals; n++)
        {
        for (size_t i=0; i<numStates; i++)
            intervalDwellTimes[n][i] = 0.0;
        
        for (size_t i=0; i<numStates; i++)
            for (size_t j=0; j<numStates; j++)
                intervalTransitions[n][i][j] = 0;
        }

    // subsitution rate (scales tree length)
    double r = getSubstitutionRate();

    // update the rate matrix
    updateRateMatrix();
    
    // uniformization
    RateMatrix* rateMatrix = q[activeRateMatrix];
    double mu = rateMatrix->uniformize(uniformizedRateMatrix);
    initializeMatrixPowers(MAX_NUM_CHANGES);

    // update the transition probabilities
    tiMngr->updateTransitionProbabilities();

    // variables to help with mappings
    std::vector<double> changeTimes;
    std::set<Change*> changesToRemove;
    int numChanges = 0;
    std::vector<double> probNumberEvents(MAX_NUM_CHANGES+1);
    
    // loop over all branches of the tree
    std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    for (int nde=0, numNodes=(int)dpSeq.size(); nde<numNodes; nde++)
        {
        Node* p = dpSeq[nde];
        if (p != tree->getRoot())
            {
            // begin branch map
            
            // get information on branch
            //TransitionProbabilities* P = p->getTransitionProbability();
            double vExact = p->getBrlenExact();
            double v = (double)(p->getBrlen()) * r; // v: integer branch length X rate
            if (v < MIN_BRLEN)
                v = MIN_BRLEN;
            double rate = mu * v;
            double branchDuration = p->getTime() - p->getAncestor()->getTime();

            // get information on states at ends of branch
            History* h = p->getHistory();
            int a = p->getAncestor()->getAreaId();
            int b = p->getAreaId();
            
            // calculate the probability of n events on the branch
            double sum = 0.0;
            for (int i=0; i<MAX_NUM_CHANGES; i++)
                {
                double nProb = (*matrixPowers[i])(a,b) * pow(rate,(double)i) * exp(-rate) / nFactorial[i];
                sum += nProb;
                probNumberEvents[i] = sum;
                }
            
            // sample number of events
            double g = sum * rng->uniformRv();
            int n = 0;
            bool foundN = false;
            for (int i=0; i<MAX_NUM_CHANGES; i++)
                {
                if (g < probNumberEvents[i])
                    {
                    n = i;
                    foundN = true;
                    break;
                    }
                }
            if (foundN == false)
                Msg::error("More than "+ std::to_string(MAX_NUM_CHANGES) + " changes along the branch");
                
            // sample the series of events
            std::vector<size_t> intermediateStates(n+1);
            intermediateStates[0] = a;
            for (size_t i=1; i<intermediateStates.size(); i++)
                {
                double sum = 0.0;
                for (size_t j=0; j<numStates; j++)
                    sum += (*matrixPowers[1])(intermediateStates[i-1],j) * (*matrixPowers[n-i])(j,b);
                double u = rng->uniformRv() * sum;
                sum = 0.0;
                for (size_t j=0; j<numStates; j++)
                    {
                    sum += (*matrixPowers[1])(intermediateStates[i-1],j) * (*matrixPowers[n-i])(j,b);
                    if (u < sum)
                        {
                        intermediateStates[i] = j;
                        break;
                        }
                    }
                }

            // sample the times
            std::vector<double> intermediateTimes;
            for (int i=0; i<n; i++)
                intermediateTimes.push_back(rng->uniformRv()*branchDuration);
            intermediateTimes.push_back(0.0);
            sort(intermediateTimes.begin(), intermediateTimes.end());
                        
            // sample changes
            h->clearHistory();
            Change* curChange = h->addChange(a, a, 0.0, p->getAncestor()->getIntervalIdx());
            curChange->anc = NULL;
            for (size_t i=1; i<intermediateStates.size(); i++)
                {
                if (intermediateStates[i-1] != intermediateStates[i])
                    {
                    double pos = p->getAncestor()->getTime() + intermediateTimes[i];
                    //double pos = p->getAncestor()->getTime() + (intermediateTimes[i] / v) * branchDuration;
                    int iid = metaData->getIntervalId(pos);
                    //std::cout << "pos=" << pos << " at=" << p->getAncestor()->getTime() << " dt=" << p->getTime() << " iid=" << iid << std::endl;
                    Change* newChange = h->addChange(intermediateStates[i-1], intermediateStates[i], intermediateTimes[i], iid);
                    newChange->anc = curChange;
                    curChange = newChange;
                    numChanges++;
                    
                    if (curChange->anc != NULL)
                        {
                        if (curChange->intervalId != curChange->anc->intervalId)
                            {
                            }
                        else
                            {
                            intervalDwellTimes[curChange->intervalId][curChange->begState] += curChange->time - curChange->anc->time;
                            }
                        }
                    }
                }
            Change* lastChange = h->addChange(curChange->endState, curChange->endState, vExact, p->getIntervalIdx());
            lastChange->anc = curChange;

            // summarize history
            std::set<Change*>& branchHistory = h->getChanges();
            for (Change* ch : branchHistory)
                {
                if (ch->begState != ch->endState)
                    {
                    if (ch->begState != ch->endState)
                        intervalTransitions[ch->intervalId][ch->begState][ch->endState]++;
                    
                    // insert code here for dwell times
                    incrementDwellTimes(h, p, ch);
                    }
                }
            incrementDwellTimes(h, p, lastChange);

            // end branch map
            }
        }
    
    for (size_t j=0; j<numStates; j++)
        {
        std::cout << std::fixed << std::setprecision(1);
        std::cout << j << " -- ";
        for (int i=0; i<numIntervals; i++)
            std::cout << std::setw(10) << intervalDwellTimes[i][j] << " ";
        std::cout << std::endl;
        }
    
    return numChanges;
}

void Model::incrementDwellTimes(History* h, Node* p, Change* change) {
  
    Change* changeAnc = change->anc;
    if (changeAnc == nullptr)
        Msg::error("Null ancestor for change in incrementDwellTimes");
    
    double pAncTime = p->getAncestor()->getTime();
    double changeTime = change->time + pAncTime;
    double changeAncTime = changeAnc->time + pAncTime;
    
    if (change->intervalId == changeAnc->intervalId)
        {
        intervalDwellTimes[change->intervalId][change->begState] += change->time - changeAnc->time;
        }
    else if (change->intervalId == 1 && changeAnc->intervalId == 0)
        {
        int boundary = 0;
        interval_map& intervalInfo = metaData->getIntervalMap();
        for (auto [key,val] : intervalInfo)
            {
            if (val == 1)
                boundary = key.first;
            }
        if (changeTime - boundary < 0.0 || boundary - changeAncTime < 0.0)
            {
            std::cout << "p->time           = " << p->getTime() << std::endl;
            std::cout << "p->anc->time      = " << p->getAncestor()->getTime() << std::endl;
            std::cout << "p->interval       = " << p->getIntervalIdx() << std::endl;
            std::cout << "p->anc->interval  = " << p->getAncestor()->getIntervalIdx() << std::endl;
            std::cout << "boundary          = " << boundary << std::endl;
            std::cout << "change->time      = " << changeTime << std::endl;
            std::cout << "change->anc->time = " << changeAncTime << std::endl;
            Msg::error("Negative times in 0");
            }
        intervalDwellTimes[1][change->begState] += changeTime - boundary;
        intervalDwellTimes[0][change->begState] += boundary - changeAncTime;
        }
    else if (change->intervalId == 0 && changeAnc->intervalId == 1)
        {
        Msg::error("01");
        }
    else if (change->intervalId == 2 && changeAnc->intervalId == 1)
        {
        int boundary = 0;
        interval_map& intervalInfo = metaData->getIntervalMap();
        for (auto [key,val] : intervalInfo)
            {
            if (val == 2)
                boundary = key.first;
            }
        if (changeTime - boundary < 0.0 || boundary - changeAncTime < 0.0)
            {
            Msg::error("Negative times in 1");
            }
        intervalDwellTimes[2][change->begState] += changeTime - boundary;
        intervalDwellTimes[1][change->begState] += boundary - changeAncTime;
        }
    else if (change->intervalId == 1 && changeAnc->intervalId == 2)
        {
        Msg::error("12");
        }
    else if (change->intervalId == 2 && changeAnc->intervalId == 0)
        {
        int boundary0 = 0, boundary1 = 0;
        interval_map& intervalInfo = metaData->getIntervalMap();
        for (auto [key,val] : intervalInfo)
            {
            if (val == 1)
                boundary0 = key.first;
            else if (val == 2)
                boundary1 = key.first;
            }
        if (changeTime - boundary1 < 0.0 || boundary1 - boundary0 < 0.0 || boundary0 - changeAncTime < 0.0)
            {
            std::cout << "boundary0         = " << boundary0 << std::endl;
            std::cout << "boundary1         = " << boundary1 << std::endl;
            std::cout << "change->time      = " << changeTime << std::endl;
            std::cout << "change->anc->time = " << changeAncTime << std::endl;
            Msg::error("Negative times in 2");
            }
        intervalDwellTimes[2][change->begState] += changeTime - boundary1;
        intervalDwellTimes[1][change->begState] += boundary1 - boundary0;
        intervalDwellTimes[0][change->begState] += boundary0 - changeAncTime;
        }
    else if (change->intervalId == 0 && changeAnc->intervalId == 2)
        {
        Msg::error("02");
        }
    else
        {
        std::cout << "Branch: " << p->getBrlen() << " " <<  p->getAncestor()->getAreaId() << "->" << p->getAreaId() << std::endl;
        std::cout << " Branch: " << p->getIntervalIdx() << " " << p->getAncestor()->getIntervalIdx() << std::endl;
        std::set<Change*>& hChanges = h->getChanges();
        for (Change* c : hChanges)
            std::cout << c->time << " " << c->begState << "->" << c->endState << " " << c->intervalId << std::endl;

        Msg::error(std::to_string(change->intervalId) + " " + std::to_string(change->anc->intervalId));
        }

}

void Model::switchActiveRateMatrix(void) {

    if (activeRateMatrix == 0)
        activeRateMatrix = 1;
    else
        activeRateMatrix = 0;
}

double Model::update(void) {

    double lnProposalProb = 0.0;

    if (modelType == jc69)
        {
        // JC69
        lnProposalProb = updateSubstitutionRate();
        }
    else if (modelType == f81)
        {
        // F81
        double u = rng->uniformRv();
        if (u < 0.25)
            {
            lnProposalProb = updateSubstitutionRate();
            }
        else 
            {
            lnProposalProb = updatePi();
            //updateRateMatrix(); // unnecessary b/c ti probabilities are analytical
            }
        }
    else if (modelType == custom_f81)
        {
        // Uganda biased F81
        double u = rng->uniformRv();
        if (u < 0.20)
            {
            lnProposalProb = updateSubstitutionRate();
            }
        else if (u < 0.80)
            {
            lnProposalProb = updatePi();
            //updateRateMatrix(); // unnecessary b/c ti probabilities are analytical
            }
        else 
            {
            if (variableUgandaRate == false)
                lnProposalProb = updateKappa();
            else 
                {
                if (rng->uniformRv() < 0.5)
                    lnProposalProb = updateKappaLockdown();
                else 
                    lnProposalProb = updateKappaNoLockdown();
                }
            //updateRateMatrix(); // unnecessary b/c ti probabilities are analytical
            }
        }
    else 
        {
        // GTR model
        double u = rng->uniformRv();
        if (u < 0.20)
            lnProposalProb = updateSubstitutionRate();
        else if (u < 0.40)
            {
            lnProposalProb = updatePi();
            updateRateMatrix();
            }
        else 
            {
            lnProposalProb = updateR();
            updateRateMatrix();
            }
        }
    
    // update the transition probabilities
    tiMngr->updateTransitionProbabilities();
        
    return lnProposalProb;
}

void Model::updateRateMatrix(void) {

    switchActiveRateMatrix();
    q[activeRateMatrix]->set(modelType, &pi[1][0], &r[1][0], kappa[1]);
}

double Model::updateSubstitutionRate(void) {

    updateType = "colonization rate";
    double oldValue = substitutionRate[0];
    double tuning = 0.005;
    double factor = exp((rng->uniformRv()-0.5)*tuning);
    double newValue = oldValue * factor;
    substitutionRate[1] = newValue;
    return log(newValue) - log(oldValue);
}

double Model::updateKappa(void) {

    updateType = "Uganda migration bias";
    double oldValue = kappa[0];
    double tuning = 0.005;
    double factor = exp((rng->uniformRv()-0.5)*tuning);
    double newValue = oldValue * factor;
    kappa[1] = newValue;
    return log(newValue) - log(oldValue);
}

double Model::updateKappaLockdown(void) {

    updateType = "Uganda migration bias during lockdown";
    double oldValue = kappaLockdown[0];
    double tuning = 0.005;
    double factor = exp((rng->uniformRv()-0.5)*tuning);
    double newValue = oldValue * factor;
    kappaLockdown[1] = newValue;
    return log(newValue) - log(oldValue);
}

double Model::updateKappaNoLockdown(void) {

    updateType = "Uganda migration bias before and after lockdown";
    double oldValue = kappaNoLockdown[0];
    double tuning = 0.005;
    double factor = exp((rng->uniformRv()-0.5)*tuning);
    double newValue = oldValue * factor;
    kappaNoLockdown[1] = newValue;
    return log(newValue) - log(oldValue);
}

static inline double randStdNormal_BoxMuller(RandomVariable* rng) {

    // one-shot Box–Muller using rng->uniformRv()
    double u1 = rng->uniformRv();
    double u2 = rng->uniformRv();

    if (u1 <= 0.0)
        u1 = std::numeric_limits<double>::min();

    const double twoPi = 6.283185307179586476925286766559;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(twoPi * u2);
}

static inline double sumLogVec(const std::vector<double>& x) {

    double s = 0.0;
    for (size_t i=0; i<x.size(); i++)
        {
        if (x[i] <= 0.0)
            return -std::numeric_limits<double>::infinity();
        s += std::log(x[i]);
        }
    return s;
}

static inline void alrForward(const std::vector<double>& x, std::vector<double>& y) {

    const size_t K = x.size();
    y.resize(K - 1);

    double xK = x[K-1];
    for (size_t i=0; i<K-1; i++)
        y[i] = std::log(x[i] / xK);
}

static inline void alrInverse(const std::vector<double>& y, std::vector<double>& x) {

    const size_t d = y.size();
    const size_t K = d + 1;

    // stable inverse: subtract max
    double m = y[0];
    for (size_t i=1; i<d; i++)
        {
        if (y[i] > m)
            m = y[i];
        }

    std::vector<double> e(d, 0.0);
    double sumE = 0.0;
    for (size_t i=0; i<d; i++)
        {
        e[i] = std::exp(y[i] - m);
        sumE += e[i];
        }

    double expNegM = std::exp(-m);
    double denom = expNegM + sumE;

    x.assign(K, 0.0);
    for (size_t i=0; i<d; i++)
        x[i] = e[i] / denom;
    x[K-1] = expNegM / denom;
}

double Model::updateSimplexALRMVN(std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal) {

    const size_t K = oldVec.size();
    if (oldVec.size() != newVec.size())
        Msg::error("updateSimplexALRMVN: unequal vector lengths");
    if (K < 2)
        Msg::error("updateSimplexALRMVN: K < 2");
    if (!(sigma > 0.0))
        Msg::error("updateSimplexALRMVN: sigma <= 0");
    if (!isValidSimplex(oldVec, 1e-10))
        Msg::error("updateSimplexALRMVN: oldVec not valid simplex");

    // forward ALR transform
    std::vector<double> y;
    alrForward(oldVec, y);

    // MVN random walk in ALR space (isotropic sigma^2 I)
    for (size_t i=0; i<y.size(); i++)
        y[i] += sigma * randStdNormal_BoxMuller(rng);

    // inverse transform
    alrInverse(y, newVec);

    // For this kernel, do NOT "normalize-with-floor" (that changes the mapping).
    // Instead, hard-reject if any component violates minVal.
    if (minVal > 0.0)
        {
        for (size_t i=0; i<newVec.size(); i++)
            {
            if (newVec[i] < minVal)
                {
                newVec = oldVec;
                return -std::numeric_limits<double>::infinity();
                }
            }
        }

    if (!isValidSimplex(newVec, 1e-10))
        Msg::error("updateSimplexALRMVN: newVec not valid simplex");

    // Hastings/Jacobian correction:
    // proposal symmetric in ALR space => only Jacobian ratio remains
    return sumLogVec(newVec) - sumLogVec(oldVec);
}

double  Model::updateSimplexALRMVN(const std::vector<double>& oldVec, std::vector<double>& newVec, double sigma, double minVal, size_t blockSize) {

    const size_t K = oldVec.size();
    if (oldVec.size() != newVec.size())
        Msg::error("updateSimplexALRMVN_blocked: unequal vector lengths");
    if (K < 2)
        Msg::error("updateSimplexALRMVN_blocked: K < 2");
    if (!(sigma > 0.0))
        Msg::error("updateSimplexALRMVN_blocked: sigma <= 0");
    if (blockSize < 2)
        Msg::error("updateSimplexALRMVN_blocked: blockSize < 2");
    if (blockSize > K)
        Msg::error("updateSimplexALRMVN_blocked: blockSize > K");
    if (!isValidSimplex(oldVec, 1e-10))
        Msg::error("updateSimplexALRMVN_blocked: oldVec not valid simplex");

    // default: start with no change
    //newVec = oldVec; // not needed because both vectors start every MCMC cycle identical

    // choose a random block of indices 
    // We want a set of distinct indices of size blockSize (Fisher-Yates partial shuffle)
    std::vector<size_t> idx(K);
    for (size_t i=0; i<K; i++)
        idx[i] = i;

    for (size_t i=0; i<blockSize; i++)
        {
        size_t j = i + (size_t)std::floor( rng->uniformRv() * (double)(K - i) );
        if (j >= K)
            j = K - 1;
        std::swap(idx[i], idx[j]);
        }

    idx.resize(blockSize);

    // extract block mass and block proportions
    double blockMass = 0.0;
    for (size_t t=0; t<blockSize; t++)
        blockMass += oldVec[idx[t]];

    if (!(blockMass > 0.0))
        {
        // should never happen for a valid simplex, but be safe
        return -1e300;
        }

    // proportions p live on an m-simplex (m = blockSize)
    std::vector<double> p_old(blockSize, 0.0);
    for (size_t t=0; t<blockSize; t++)
        p_old[t] = oldVec[idx[t]] / blockMass;

    // if you enforce a global minVal, then inside the block the implied minimum on p is minVal / blockMass.
    double pMin = 0.0;
    if (minVal > 0.0)
        pMin = minVal / blockMass;

    // propose new block proportions with your existing ALR-MVN code
    std::vector<double> p_new(blockSize, 0.0);

    /* This returns the Hastings correction for the move in p-space.
       Because x_block = blockMass * p, the constant blockMass cancels in Jacobians and in proposal densities.
       So the Hastings correction for the overall x-move is IDENTICAL to the one for p. */
    double logH = updateSimplexALRMVN(p_old, p_new, sigma, pMin);

    // if updateALRMVN hard-rejected (returns LOG_ZERO), propagate rejection.
    if (logH <= -1e250)
        {
        newVec = oldVec;
        return logH;
        }

    // write back block and validate global minVal if requested
    for (size_t t=0; t<blockSize; t++)
        newVec[idx[t]] = blockMass * p_new[t];

    if (minVal > 0.0)
        {
        for (size_t i=0; i<newVec.size(); i++)
            {
            if (newVec[i] < minVal)
                {
                newVec = oldVec;
                return -1e300;
                }
            }
        }

    if (!isValidSimplex(newVec, 1e-10))
        Msg::error("updateSimplexALRMVN_blocked: newVec not valid simplex");

    return logH;
}

double Model::updateSimplexPrior(const std::vector<double>& oldVec, std::vector<double>& newVec, std::vector<double>& alpha, double minVal) {

    Probability::Dirichlet::rv(rng, alpha, newVec);
    return Probability::Dirichlet::lnPdf(alpha, oldVec) - Probability::Dirichlet::lnPdf(alpha, newVec);
}

static inline bool isValidSimplex(const std::vector<double>& x, double eps) {

    double s = 0.0;
    for (size_t i=0; i<x.size(); i++)
        {
        if (!(x[i] > 0.0))
            return false;
        s += x[i];
        }
    return std::fabs(s - 1.0) < eps;
}

double Model::updatePi(void) {

    double minVal = 0.00001;

    // Mixture of simplex kernels:
    //  - ALR MVN: good for informative posteriors (tight, correlated)
    //  - Dirichlet subset move (k=1): safe local exploration
    //  - Mass transfer: cheap polishing move
    double u = rng->uniformRv();
    if (u < 0.34)
        {
        // tune sigma to get reasonable acceptance (often ~0.05–0.20)
        updateType = "equilibrium frequencies: ALR MVN";
       return updateSimplexALRMVN(pi[0], pi[1], 0.001, minVal);
        }
    else if (u < 0.67)
        {
        updateType = "equilibrium frequencies: Dirichlet";
        return updateSimplex(pi[0], pi[1], 20000.0, 1, minVal);
        }
    else
        {
        updateType = "equilibrium frequencies: mass transfer";
        return updateSimplexTransfer(pi[0], pi[1], 300.0, minVal);
        }
}

double Model::updateR(void) {

    double minVal = 0.000001;

    // Mixture of simplex kernels:
    //  - ALR MVN: good for informative posteriors (tight, correlated)
    //  - Dirichlet subset move (k=1): safe local exploration
    //  - Mass transfer: cheap polishing move
    double u = rng->uniformRv();
    if (u < 0.34)
        {
        // tune sigma to get reasonable acceptance (often ~0.05–0.20)
        if (r[0].size() > 100)
            {
            updateType = "exchangeability rates: blocked ALR MVN";
            return updateSimplexALRMVN(r[0], r[1], 0.001, minVal, 100);
            }
        else 
            {
            updateType = "exchangeability rates: ALR MVN";
            return updateSimplexALRMVN(r[0], r[1], 0.001, minVal);
            }
        }
    else if (u < 0.67)
        {
        updateType = "exchangeability rates: Dirichlet";
        return updateSimplex(r[0], r[1], 10000.0, 1, minVal);
        }
    else
        {
        updateType = "exchangeability rates: mass transfer";
        return updateSimplexTransfer(r[0], r[1], 300.0, minVal);
        }
}

double Model::updateSimplexTransfer(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal) {

    const size_t K = oldVec.size();
    if (oldVec.size() != newVec.size())
        Msg::error("proposeSimplex_MassTransfer: unequal vector lengths");
    if (K < 2)
        Msg::error("proposeSimplex_MassTransfer: K < 2");
    if (!(alpha0 > 0.0))
        Msg::error("proposeSimplex_MassTransfer: alpha0 <= 0");
    if (!isValidSimplex(oldVec, 1e-10))
        Msg::error("proposeSimplex_MassTransfer: x_curr not valid simplex");

    // pick a random pair (i,j), i<j
    size_t i = (size_t)std::floor(rng->uniformRv() * (double)K);
    size_t j = (size_t)std::floor(rng->uniformRv() * (double)(K-1));
    if (j >= i)
        j += 1;
    if (j < i)
        std::swap(i, j);

    double xi = oldVec[i];
    double xj = oldVec[j];
    double s  = xi + xj;

    if (!(s > 0.0))
        Msg::error("proposeSimplex_MassTransfer: s <= 0 (should not happen)");

    double u = xi / s;

    // Beta parameters centered at u
    // To avoid a,b becoming too tiny (can cause numerical ugliness), clamp u away from 0/1 a hair.
    const double eps = 1e-15;
    double uClamped = std::min(1.0 - eps, std::max(eps, u));

    double a = alpha0 * uClamped;
    double b = alpha0 * (1.0 - uClamped);

    double uProp = Probability::Beta::rv(rng, a, b);

    // construct proposal
    newVec[i] = s * uProp;
    newVec[j] = s * (1.0 - uProp);

    // (tiny numerical safety) renormalize in case of roundoff
    Probability::Helper::normalize(newVec, minVal);

    // reverse proposal is Beta centered at u_prop (from the proposed state)
    double uRev = newVec[i] / (newVec[i] + newVec[j]); // should equal u_prop up to rounding
    double uRevClamped = std::min(1.0 - eps, std::max(eps, uRev));

    double aRev = alpha0 * uRevClamped;
    double bRev = alpha0 * (1.0 - uRevClamped);

    // Hastings ratio (reverse - forward)
    // Forward density: u_prop ~ Beta(a,b)
    // Reverse density: u ~ Beta(a_rev, b_rev)
    double log_q_fwd = Probability::Beta::lnPdf(a, b, uProp);
    double log_q_rev = Probability::Beta::lnPdf(aRev, bRev, u);
    double lnHastings = log_q_rev - log_q_fwd;

    if (!isValidSimplex(newVec, 1e-10))
        Msg::error("proposeSimplex_MassTransfer: x_prop invalid");
    return lnHastings;
}

double Model::updateSimplex(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, double minVal) {
    
    std::vector<double> alphaForward(oldVec.size());
    for (size_t i=0, n=alphaForward.size(); i<n; i++)
        alphaForward[i] = oldVec[i] * alpha0 + 1.0;
        
    Probability::Dirichlet::rv(rng, alphaForward, newVec);
    Probability::Helper::normalize(newVec, minVal);
    
    std::vector<double> alphaReverse(newVec.size());
    for (size_t i=0, n=alphaReverse.size(); i<n; i++)
        alphaReverse[i] = newVec[i] * alpha0 + 1.0;

    double lnForwardProb = Probability::Dirichlet::lnPdf(alphaReverse, oldVec) - Probability::Dirichlet::lnPdf(alphaForward, newVec);
    return lnForwardProb;
}

double Model::updateSimplex(std::vector<double>& oldVec, std::vector<double>& newVec, double alpha0, size_t k, double minVal) {
    
    // choose which elements to update
    size_t n = oldVec.size();
    if (n == k)
        Msg::error("Cannot update that many simplex variables!!!");
        
    std::set<size_t> indices;
    while (indices.size() < k)
        {
        size_t idx = (size_t)(rng->uniformRv() * n);
        indices.insert(idx);
        }
    
    // fill in old values for the selected indices plus the "remainder" element
    std::vector<double> oldValues(k+1);
    double sum = 0.0;
    int i = 0;
    for (size_t idx : indices)
        {
        double x = oldVec[idx];
        oldValues[i++] = x;
        sum += x;
        }
    oldValues[i] = 1.0 - sum;  // remainder element
    
    // construct the forward Dirichlet parameters
    std::vector<double> alphaForward(k + 1);
    for (size_t j=0, m=alphaForward.size(); j<m; j++)
        alphaForward[j] = oldValues[j] * alpha0 + 1.0;
        
    // draw new values from Dirichlet
    std::vector<double> newValues(k+1);
    Probability::Dirichlet::rv(rng, alphaForward, newValues);
    
    // FIX 1: normalize newValues (the k+1 element vector), NOT newVec (the full n-element vector)
    Probability::Helper::normalize(newValues, minVal);
    
    // construct the reverse Dirichlet parameters
    std::vector<double> alphaReverse(k + 1);
    for (size_t j=0, m=alphaReverse.size(); j<m; j++)
        alphaReverse[j] = newValues[j] * alpha0 + 1.0;
        
    // fill in the full newVec from the updated (newValues) vector
    // The non-selected elements are scaled by the ratio of remainders
    double factor = newValues[k] / oldValues[k];
    
    // FIX 2: Initialize newVec from oldVec explicitly, then scale
    for (size_t j=0, m=newVec.size(); j<m; j++)
        newVec[j] = oldVec[j] * factor;
    
    // overwrite the selected indices with their new values
    i = 0;
    for (size_t idx : indices)
        newVec[idx] = newValues[i++];

    // FIX 3: Compute Hastings ratio using oldValues and newValues (both size k+1),
    // NOT oldVec and newVec (which are size n and would cause dimension mismatch)
    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldValues) - Probability::Dirichlet::lnPdf(alphaForward, newValues);
    lnProposalProb += (n - k - 1) * log(factor);
    return lnProposalProb;
}
