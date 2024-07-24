#include "History.hpp"



History::~History(void) {

    for (int i=0; i<changesPool.size(); i++)
        delete changesPool[i];
}

void History::addChange(int a, int b, double t) {

    numChanges++;
    if (changesPool.size() < numChanges)
        {
        for (int i=0; i<numChanges-changesPool.size(); i++)
            changesPool.push_back(new Change);
        }
    
    Change* chg = changesPool[changesPool.size()-1];
    changesPool.pop_back();
    chg->begState = a;
    chg->endState = b;
    chg->time = t;
    changes.insert(chg);
}

void History::clearHistory(void) {
    
    numChanges = 0;
    for (Change* c : changes)
        changesPool.push_back(c);
    changes.clear();
}

void History::removeChanges(std::set<Change*>& goners) {

    for (Change* c : goners)
        {
        changesPool.push_back(c);
        changes.erase(c);
        }
    numChanges -= goners.size();
}
