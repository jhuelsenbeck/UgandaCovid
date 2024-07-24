#ifndef History_hpp
#define History_hpp

struct Change {

    int begState;
    int endState;
    double time;
};

#include <set>
#include <vector>



class History {

    public:
                               ~History(void);
        void                    addChange(int a, int b, double t);
        void                    clearHistory(void);
        std::set<Change*>&      getChanges(void) { return changes; }
        int                     getNumChanges(void) { return numChanges; }
        void                    removeChanges(std::set<Change*>& goners);
    
    private:
        int                     numChanges;
        std::set<Change*>       changes;
        std::vector<Change*>    changesPool;
};

#endif
