#ifndef History_hpp
#define History_hpp

struct Change {

    int begState;
    int endState;
    double time;
    int intervalId;
    Change* anc;
};

#include <set>
#include <vector>



class History {

    public:
                                History(void);  // Initialize numChanges
                               ~History(void);
        Change*                 addChange(int a, int b, double t, int iid);
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
