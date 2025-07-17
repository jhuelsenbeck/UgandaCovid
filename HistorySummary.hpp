#ifndef HistorySummary_hpp
#define HistorySummary_hpp

#include <vector>
#include "Container.hpp"
class Tree;



class HistorySummary {

    public:
                                HistorySummary(void) = delete;
                                HistorySummary(int ns, int ti);
                               ~HistorySummary(void);
        void                    addSample(int** counts);
        void                    summary(void);
    
    private:
        int                     numSamples;
        int                     numStates;
        int                     timeInterval;
        std::vector<int>        numChanges;
        std::vector<double>     treeLength;
        std::vector<IntMatrix*> transitionCount;
};

#endif
