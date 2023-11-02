#include "Msg.hpp"
#include "Node.hpp"
#include "Tree.hpp"
#include <cmath>
#include <fstream>
#include <iostream>


Tree::Tree(const Tree& t) {

    clone(t);
}

Tree::Tree(std::string fileName) {

    numNodes = 0;
    numTips = 0;

    // open file
    std::fstream treeFile(fileName, std::ios::in);
    if (treeFile.is_open() == false)
        Msg::error("Could not open file \"" + fileName + "\"");
    
    // read the file's contents
    std::string newickStr = "";
    std::string str = "";
    while (getline(treeFile, str))
        {
        newickStr += str;
        str = "";
        }

    // close the file
    treeFile.close();

    // read the Newick string information
    Node* p = nullptr;
    bool readingBranchLength = false;
    for (int i=0, n=(int)newickStr.length(); i<n; i++)
        {
        char c = newickStr[i];
        
        if (c == '(')
            {
            Node* newNode = addNode();
            if (p == nullptr)
                {
                root = newNode;
                }
            else
                {
                p->addDescendant(newNode);
                newNode->setAncestor(p);
                }
            p = newNode;
            readingBranchLength = false;
            }
        else if (c == ')' || c == ',')
            {
            if (p->getAncestor() == nullptr)
                Msg::error("Can't move down a node");
            p = p->getAncestor();
            readingBranchLength = false;
            }
        else if (c == ':')
            {
            readingBranchLength = true;
            }
        else if (c == ';')
            {
            if (p != root)
                Msg::error("Expected to end at the root of the tree");
            }
        else
            {
            std::string token = "";
            readAhead(token, newickStr, i);
            if (readingBranchLength == false)
                {
                // add a taxon
                Node* newNode = addNode();
                newNode->setName(token);
                newNode->setIsTip(true);
                p->addDescendant(newNode);
                newNode->setAncestor(p);
                p = newNode;
                numTips++;
                }
            else
                {
                // add a branch length
                double x = atof(token.c_str());
                int xI = round(x);
                p->setBrlen(xI);
                readingBranchLength = false;
                }
            }
        }
        
    initializeDownPassSequence();
    
    std::cout << "Successfully read tree file \"" << fileName << "\"" << std::endl;
    std::cout << "Tree has " << numTips << " tips and " << numNodes << " nodes" << std::endl;
#   if 0
    std::vector<int> descDist;
    for (auto p=downPassSequence.begin(); p != downPassSequence.end(); p++)
        {
        if ((*p)->getIsTip() == false)
            {
            int nd = (*p)->getNumDescendants();
            if (nd+1 > descDist.size())
                descDist.resize(nd+10);
            descDist[nd]++;
            }
        }
    for (int i=0; i<descDist.size(); i++)
        std::cout << i << " -- " << descDist[i] << std::endl;
#   endif
}

Tree::~Tree(void) {

    // free memory here
    deleteNodes();
}

Tree& Tree::operator=(const Tree& rhs) {

    if (this != &rhs)
        clone(rhs);
    return *this;
}

Node* Tree::addNode(void) {

    Node* newNode = new Node;
    if (newNode == nullptr)
        Msg::error("Failed to allocate node");
    newNode->setOffset((int)nodes.size());
    nodes.push_back(newNode);
    newNode->setIndex(numNodes);
    numNodes++;
    return newNode;
}

void Tree::clone(const Tree& t) {

    // set up nodes for this tree
    if (this->numNodes != t.numNodes)
        {
        deleteNodes();
        for (int i=0; i<t.numNodes; i++)
            addNode();
        }
    this->numNodes = t.numNodes;
    this->numAreas = t.numAreas;
    
    // specify the root
    this->root = nodes[t.root->getOffset()];
    
    // deep copy of node information
    for (int i=0; i<this->numNodes; i++)
        {
        Node* p = nodes[i];
        Node* q = t.nodes[i];
        
        p->setIndex(q->getIndex());
        p->setIsTip(q->getIsTip());
        p->setBrlen(q->getBrlen());
        p->setName(q->getName());
        p->setAreaId(q->getAreaId());
        
        if (q->getAncestor() != nullptr)
            p->setAncestor( nodes[q->getAncestor()->getOffset()] );
        else
            p->setAncestor(nullptr);
        p->removeAllDescendants();
        std::set<Node*>& qDesc = q->getDescendants();
        for (Node* r : qDesc)
            p->addDescendant( nodes[r->getOffset()] );

        if (p->getConditionalLikelihood() == nullptr)
            {
            double* x = new double[this->numAreas];
            p->setConditionalLikelihood(x);
            memcpy(x, q->getConditionalLikelihood(), numAreas*sizeof(double));
            }
        }
        
    this->downPassSequence.resize(t.downPassSequence.size());
    for (int i=0,n=(int)t.downPassSequence.size(); i<n; i++)
        this->downPassSequence[i] = nodes[t.downPassSequence[i]->getOffset()];
}

void Tree::deleteNodes(void) {

    for (int i=0, n=(int)nodes.size(); i<n; i++)
        delete nodes[i];
    nodes.clear();
}

void Tree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    interiorDownPassSequence.clear();
    passDown(root);
    endDownPassSequence = downPassSequence[0] + downPassSequence.size();
}

void Tree::passDown(Node* p) {

    if (p != nullptr)
        {
        std::set<Node*>& pDesc = p->getDescendants();
        for (Node* d : pDesc)
            passDown(d);
        downPassSequence.push_back(p);
        if (p->getIsTip() == false)
            interiorDownPassSequence.push_back(p);
        }
}

void Tree::print(void) {

    showNode(root, 0);
}

void Tree::readAhead(std::string& token, std::string& newickStr, int& i) {

    int j = i;
    bool goodChar = true;
    do
        {
        char c = newickStr[j];
        if (c != '(' && c != ')' && c != ',' && c != ';' && c != ':')
            {
            token += std::string(1,c);
            goodChar = true;
            j++;
            }
        else
            {
            goodChar = false;
            j--;
            }
        } while (goodChar == true);
    i = j;
}

void Tree::showNode(Node* p, int indent) {

    if (p != nullptr)
        {
        std::set<Node*>& pDesc = p->getDescendants();
        
        for (int i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex() << " ";
        std::cout << std::endl;
        
        for (Node* d : pDesc)
            showNode(d, indent + 1);
        }
}
