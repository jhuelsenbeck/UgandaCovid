#include "Node.hpp"



Node::Node(void) {

    ancestor    = nullptr;
    index       = 0;
    offset      = 0;
    brlen       = 0.0;
    name        = "";
    areaId      = 0;
    isTip       = false;
    cl          = nullptr;
    clEnd       = nullptr;
    tiProb      = nullptr;
}

Node::~Node(void) {

}
