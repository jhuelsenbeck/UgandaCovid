#include <regex>
#include "Node.hpp"



Node::Node(void) {

    ancestor     = nullptr;
    history      = nullptr;
    index        = 0;
    offset       = 0;
    brlen        = 0;
    brlenExact   = 0.0;
    time         = 0.0;
    name         = "";
    areaId       = 0;
    areaName     = "";
    isAreaFixed  = false;
    isTip        = false;
    cl           = nullptr;
    clEnd        = nullptr;
    tiProb       = nullptr;
    scratchInt   = 0;
    scratchBool  = false;
    job          = nullptr;
    dependentJob = nullptr;
    goodTime     = false;
    intervalIdx  = 0;
}

Node::~Node(void) {

}

Node* Node::getFirstDescendant(void) {
    
    return *descendants.begin();
}

double Node::getOldestDescendant(void) {

    double oldestTime = -1.0;
    for (Node* d : descendants)
        {
        if (d->getGoodTime() == true)
            {
            if (oldestTime < 0.0)
                oldestTime = d->getTime();
            else if (d->getTime() < oldestTime)
                oldestTime = d->getTime();
            }
        }
    return oldestTime;
}

void Node::setName(std::string s) {

    //std::string newS = std::regex_replace(s, std::regex("_"), " ");
    name = s;
}
