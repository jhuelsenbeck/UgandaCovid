#include <regex>
#include "Node.hpp"



Node::Node(void) : areaName(""), name("") {

    ancestor     = nullptr;
    firstChild   = nullptr;
    nextSibling  = nullptr;
    history      = nullptr;
    index        = 0;
    offset       = 0;
    brlen        = 0;
    brlenExact   = 0.0;
    time         = 0.0;
    areaId       = 0;
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

// Add a child to this node (append to end of sibling list)
void Node::addDescendant(Node* p) {
    
    if (firstChild == nullptr) {
        // No children yet - this becomes the first child
        firstChild = p;
    } else {
        // Find the last sibling and append
        Node* lastChild = firstChild;
        while (lastChild->nextSibling != nullptr)
            lastChild = lastChild->nextSibling;
        lastChild->nextSibling = p;
    }
    p->nextSibling = nullptr;  // Ensure new child has no sibling
}

// Remove a specific child from this node
void Node::removeDescendant(Node* p) {
    
    if (firstChild == nullptr)
        return;
    
    if (firstChild == p) {
        // Removing the first child
        firstChild = p->nextSibling;
        p->nextSibling = nullptr;
        return;
    }
    
    // Search through siblings to find predecessor
    Node* prev = firstChild;
    while (prev->nextSibling != nullptr && prev->nextSibling != p)
        prev = prev->nextSibling;
    
    if (prev->nextSibling == p) {
        prev->nextSibling = p->nextSibling;
        p->nextSibling = nullptr;
    }
}

// Remove all children from this node
void Node::removeAllDescendants(void) {
    
    // Clear sibling links for all children
    Node* child = firstChild;
    while (child != nullptr) {
        Node* next = child->nextSibling;
        child->nextSibling = nullptr;
        child = next;
    }
    firstChild = nullptr;
}

// Count the number of children (O(n) but rarely called in hot path)
int Node::getNumDescendants(void) {
    
    int count = 0;
    for (Node* d = firstChild; d != nullptr; d = d->nextSibling)
        count++;
    return count;
}

double Node::getOldestDescendant(void) {

    double oldestTime = -1.0;
    for (Node* d = firstChild; d != nullptr; d = d->nextSibling)
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
