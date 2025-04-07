#ifndef HistorySummary_hpp
#define HistorySummary_hpp

#include <vector>
#include "Container.hpp"
class Tree;



class HistorySummary {

    public:
                            HistorySummary(void) = delete;
                            HistorySummary(int ns);
                           ~HistorySummary(void);
        void                addSample(Tree& t);
        void                summary(void);
    
    private:
        int                 numSamples;
        int                 numStates;
        std::vector<int>    numChanges;
        IntMatrix*          transitionCount;
};

#endif
