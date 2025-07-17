#include "History.hpp"
#include "HistorySummary.hpp"
#include "Node.hpp"
#include "Tree.hpp"



HistorySummary::HistorySummary(int ns, int ti) {

    numStates = ns;
    timeInterval = ti;
    numSamples = 0;
}

HistorySummary::~HistorySummary(void) {

    for (int i=0; i<transitionCount.size(); i++)
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
    int ugundaId = 106;
    double aveInto = 0.0;
    double aveOutof = 0.0;
    double aveNumChanges = 0.0;
    
    for (int n=0; n<transitionCount.size(); n++)
        {
        int numIntoUgunda = 0;
        int numOutofUgunda = 0;
        int nc = 0;
        for (int i=0; i<numStates; i++)
            {
            if (i != ugundaId)
                {
                numIntoUgunda += (*transitionCount[n])(i,ugundaId);
                numOutofUgunda += (*transitionCount[n])(ugundaId,i);
                }
            for (int j=0; j<numStates; j++)
                {
                if (i != j)
                    nc += (*transitionCount[n])(i,j);
                }
            }
        aveInto += numIntoUgunda;
        aveOutof += numOutofUgunda;
        aveNumChanges += nc;
        }
    
    std::cout << "   Average number number changes = " << aveNumChanges / transitionCount.size() << std::endl;
    std::cout << "   Average number into Ugunda    = " << aveInto / transitionCount.size() << std::endl;
    std::cout << "   Average number out of Ugunda  = " << aveOutof / transitionCount.size() << std::endl;
}
