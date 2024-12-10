#include "History.hpp"
#include "HistorySummary.hpp"
#include "Node.hpp"
#include "Tree.hpp"



HistorySummary::HistorySummary(int ns) {

    numStates = ns;
    numSamples = 0;
    transitionCount = new IntMatrix(numStates, numStates);
}

HistorySummary::~HistorySummary(void) {

    delete transitionCount;
}

void HistorySummary::addSample(Tree& t) {

    numSamples++;
    
    int n = 0;
    std::vector<Node*>& dpSeq = t.getDownPassSequence();
    for (Node* p : dpSeq)
        {
        if (p != t.getRoot())
            {
            History* h = p->getHistory();
            std::set<Change*>& branchChanges = h->getChanges();
            for (Change* c : branchChanges)
                {
                (*transitionCount)(c->begState, c->endState)++;
                n++;
                }
            }
        }
    
    numChanges.push_back(n);
}

void HistorySummary::summary(void) {

    int n = 0;
    for (int x : numChanges)
        n += x;

    std::cout << "Number of samples = " << numSamples << std::endl;
    std::cout << "Average number of changes = " << (double)n / numSamples << std::endl;
    
    std::multimap<double, std::pair<int,int>> chgs;
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                double key = (double)(*transitionCount)(i,j) / numSamples;
                chgs.insert( std::make_pair(key, std::make_pair(i,j)) );
                }
            }
        }
        
    for (std::multimap<double, std::pair<int,int>>::reverse_iterator it = chgs.rbegin(); it != chgs.rend(); it++)
        std::cout << it->second.first << " -> " << it->second.second << " -- " << it->first << std::endl;
}
