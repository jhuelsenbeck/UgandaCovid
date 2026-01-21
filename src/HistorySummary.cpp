#include "History.hpp"
#include "HistorySummary.hpp"
#include "Node.hpp"
#include "Tree.hpp"



HistorySummary::HistorySummary(int ns) : numSamples(ns) {

    numSamples = 0;
}

HistorySummary::~HistorySummary(void) {

    for (size_t i=0; i<transitionCount.size(); i++)
        delete transitionCount[i];
}

void HistorySummary::addSample(int** counts) {

    IntMatrix* cnt = new IntMatrix(numStates, numStates);
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            (*cnt)(i,j) = counts[i][j];
            }
        }
    transitionCount.push_back(cnt);
}

void HistorySummary::summary(void) {
        
    // calculate averages for number of changes
    int ugandaId = 106;
    double aveInto = 0.0;
    double aveOutof = 0.0;
    double aveNumChanges = 0.0;
    
    for (size_t n=0; n<transitionCount.size(); n++)
        {
        int numIntoUganda = 0;
        int numOutofUganda = 0;
        int nc = 0;
        for (int i=0; i<numStates; i++)
            {
            if (i != ugandaId)
                {
                numIntoUganda += (*transitionCount[n])(i,ugandaId);
                numOutofUganda += (*transitionCount[n])(ugandaId,i);
                }
            for (int j=0; j<numStates; j++)
                {
                if (i != j)
                    nc += (*transitionCount[n])(i,j);
                }
            }
        aveInto += numIntoUganda;
        aveOutof += numOutofUganda;
        aveNumChanges += nc;

        if (n == transitionCount.size() - 1)
            {
            std::cout << "   Number changes = " << nc << std::endl;
            std::cout << "   Number into Uganda = " << numIntoUganda << std::endl;
            std::cout << "   Number out of Uganda = " << numOutofUganda << std::endl;
            }
        }

    std::cout << "   Average number changes = " << aveNumChanges / transitionCount.size() << std::endl;
    std::cout << "   Average number into Uganda    = " << aveInto / transitionCount.size() << std::endl;
    std::cout << "   Average number out of Uganda  = " << aveOutof / transitionCount.size() << std::endl;
}
