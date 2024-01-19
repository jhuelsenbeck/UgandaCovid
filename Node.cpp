#include <regex>
#include "Node.hpp"



Node::Node(void) {

    ancestor     = nullptr;
    index        = 0;
    offset       = 0;
    brlen        = 0.0;
    name         = "";
    areaId       = 0;
    isTip        = false;
    cl           = nullptr;
    clEnd        = nullptr;
    tiProb       = nullptr;
    scratch      = 0;
    job          = nullptr;
    dependentJob = nullptr;
}

Node::~Node(void) {

}

void Node::setName(std::string s) {

    //std::string newS = std::regex_replace(s, std::regex("_"), " ");
    name = s;
}
