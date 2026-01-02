#include "Msg.hpp"
#include "Node.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>


Tree::Tree(const Tree& t) {

    clone(t);
}

Tree::Tree(Node* rootPattern) {

    std::map<Node*,Node*> nodeMap;
    root = copy(rootPattern, nodeMap);
    initializeDownPassSequence();
    print();
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
                p->setBrlenExact(x);
                readingBranchLength = false;
                }
            }
        }
        
    initializeDownPassSequence();
        
    std::cout << "     Successfully read tree file \"" << fileName << "\"" << std::endl;
    std::cout << "     Tree has " << numTips << " tips and " << numNodes << " nodes" << std::endl;
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
        // LCRS iteration over children
        for (Node* r = q->getFirstChild(); r != nullptr; r = r->getNextSibling())
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

Node* Tree::copy(Node* pPattern, std::map<Node*,Node*>& nodeMap) {

    if (pPattern == nullptr)
        return nullptr;
    
    std::map<Node*,Node*>::iterator it = nodeMap.find(pPattern);
    if (it != nodeMap.end())
        return it->second;

    Node* newP = addNode();
    nodeMap.insert( std::make_pair(pPattern, newP) );
    newP->setIndex(pPattern->getIndex());
    newP->setBrlen(pPattern->getBrlen());
    newP->setAreaId(pPattern->getAreaId());
    newP->setName(pPattern->getName());
    newP->setIsTip(pPattern->getIsTip());
    newP->setAncestor( copy(pPattern->getAncestor(), nodeMap) );

    // LCRS iteration over children
    for (Node* d = pPattern->getFirstChild(); d != nullptr; d = d->getNextSibling())
        newP->addDescendant( copy(d, nodeMap) );

    return newP;
}

void Tree::deleteNodes(void) {

    for (int i=0, n=(int)nodes.size(); i<n; i++)
        delete nodes[i];
    nodes.clear();
}

Node* Tree::findTaxonNamed(std::string tName) {
    
    for (Node* p : downPassSequence)
        {
        if (p->getIsTip() == true)
            {
            if (p->getName() == tName)
                return p;
            }
        }
    return nullptr;
}

std::string Tree::getNewickString(void) {

    std::stringstream strm;
    writeTree(root, strm);
    strm << ";";
    return strm.str();
}

void Tree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    interiorDownPassSequence.clear();
    passDown(root);
}

void Tree::passDown(Node* p) {

    if (p != nullptr)
        {
        // LCRS iteration over children
        for (Node* d = p->getFirstChild(); d != nullptr; d = d->getNextSibling())
            passDown(d);
        downPassSequence.push_back(p);
        if (p->getIsTip() == false)
            interiorDownPassSequence.push_back(p);
        }
}

void Tree::print(void) {

    showNode(root, 0);
}

void Tree::print(Node* subtree) {

    showNode(subtree, 0);
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

void Tree::removeNodes(std::vector<Node*>& nodesToRemove) {
    
    // mark tips to be removed
    int numNodesMarkedForRemoval = 0;
    for (Node* p : downPassSequence)
        p->setScratchBool(false);
    for (Node* p : nodesToRemove)
        {
        p->setScratchBool(true);
        numNodesMarkedForRemoval++;
        if (p->getIsTip() == false)
            Msg::error("Attempting to remove an interior node!");
        }
    
    // mark interior nodes to be removed
    int numTipsBefore = 0;
    for (Node* p : downPassSequence)
        {
        if (p->getIsTip() == false)
            {
            // LCRS iteration over children
            int numDescendantsMarked = 0;
            int numDescendantsUnmarked = 0;
            for (Node* d = p->getFirstChild(); d != nullptr; d = d->getNextSibling())
                {
                if (d->getScratchBool() == true)
                    numDescendantsMarked++;
                else
                    numDescendantsUnmarked++;
                }
            if (numDescendantsUnmarked == 0 || numDescendantsUnmarked == 1)
                {
                p->setScratchBool(true);
                numNodesMarkedForRemoval++;
                }
            }
        else
            {
            numTipsBefore++;
            }
        }
    
    // check marking pattern
    for (Node* p : downPassSequence)
        {
        if (p->getScratchBool() == true)
            {
            // LCRS iteration over children
            int numUnmarkedDescendants = 0;
            for (Node* d = p->getFirstChild(); d != nullptr; d = d->getNextSibling())
                {
                if (d->getScratchBool() == false)
                    numUnmarkedDescendants++;
                }
            if (p->getIsTip() == true && numUnmarkedDescendants != 0)
                Msg::error("Should not have any descendants for a tip node");
            else if (p->getIsTip() == false && numUnmarkedDescendants > 1)
                Msg::error("Should have only 0 or 1 descendant unmarked for removal");
            }
        }
    
    // remove marked nodes
    std::vector<double> treeLengthBefore;
    int totalRemoved = 0;
    int numRemovedInPass = 0;
    std::set<Node*> removedNodes;
    do
        {
        treeLengthBefore.push_back(0.0);
        size_t idx = treeLengthBefore.size()-1;
        numRemovedInPass = 0;
        for (Node* p : downPassSequence)
            {
            if (p->getScratchBool() == true)
                {
                Node* pAnc = p->getAncestor();
                if (pAnc == nullptr)
                    Msg::error("Ancestor is null in removeNodes");
                if (p->getNumDescendants() == 0)
                    {
                    pAnc->removeDescendant(p);
                    p->setAncestor(nullptr);
                    p->setScratchBool(false);
                    removedNodes.insert(p);
                    numRemovedInPass++;
                    }
                else if (p->getNumDescendants() == 1)
                    {
                    Node* d = p->getFirstDescendant();
                    if (d == nullptr)
                        Msg::error("First descendant of p is null");
                    treeLengthBefore[idx] += p->getBrlenExact();
                    d->setBrlenExact( p->getBrlenExact() + d->getBrlenExact() );
                    pAnc->removeDescendant(p);
                    pAnc->addDescendant(d);
                    d->setAncestor(pAnc);
                    p->removeAllDescendants();
                    p->setAncestor(nullptr);
                    p->setScratchBool(false);
                    p->setBrlenExact(0.0);
                    removedNodes.insert(p);
                    numRemovedInPass++;
                    }
                else
                    {
                    treeLengthBefore[idx] += p->getBrlenExact();
                    }
                }
            else
                {
                if (p->getAncestor() != nullptr)
                    treeLengthBefore[idx] += p->getBrlenExact();
                }
            }
        totalRemoved += numRemovedInPass;
        initializeDownPassSequence();
        } while (numRemovedInPass > 0);
    
    std::cout << "     Removed " << totalRemoved << " nodes out of a total of " << numNodesMarkedForRemoval << " originally marked for removal" << std::endl;
    
    initializeDownPassSequence();
        
    // check that we don't have any weird things in the tree
    int numTipsAfter = 0;
    int numWeirdInterior = 0;
    double treeLengthAfter = 0.0;
    for (Node* p : downPassSequence)
        {
        if (p->getIsTip() == true)
            {
            numTipsAfter++;
            }
        else
            {
            if (p->getNumDescendants() < 2)
                numWeirdInterior++;
            }
        if (p->getAncestor() != nullptr)
            treeLengthAfter += p->getBrlenExact();
        }
    
    // set the number of tips
    numTips = numTipsAfter;
    
    // set the number of nodes
    nodes = downPassSequence;
    numNodes = (int)nodes.size();
    for (int i=0, n=(int)nodes.size(); i<n; i++)
        nodes[i]->setOffset(i);
    
    // free the removed nodes
    for (Node* delNode : removedNodes)
        delete delNode;
    
    if (numWeirdInterior > 0)
        Msg::error("We have " + std::to_string( numWeirdInterior) + " weird interior nodes");
    
    std::cout << "     Number of tips before removal = " << numTipsBefore << std::endl;
    std::cout << "     Number of tips after removal  = " << numTipsAfter << std::endl;
    std::cout << "     Tree length before removal    = ";
    for (int i=0; i<treeLengthBefore.size(); i++)
        std::cout << treeLengthBefore[i] << " ";
    std::cout << std::endl;
    std::cout << "     Tree length after removal     = " << treeLengthAfter << std::endl;
    std::cout << "     Tree length difference        = ";
    for (int i=0; i<treeLengthBefore.size(); i++)
        std::cout << treeLengthBefore[i] - treeLengthAfter << " ";
    std::cout << std::endl;
    std::cout << "     " << numTipsAfter << " + " << nodesToRemove.size() << " = " << numTipsAfter + nodesToRemove.size() << std::endl;
}


void Tree::showNode(Node* p, int indent) {

    if (p != nullptr)
        {
        for (int i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex();
        if (p->getAncestor() != nullptr)
            std::cout << " a_" << p->getAncestor()->getIndex() << " ( ";
        // LCRS iteration over children
        for (Node* d = p->getFirstChild(); d != nullptr; d = d->getNextSibling())
            std::cout << d->getIndex() << " ";
        std::cout << ") ";
        std::cout << p->getBrlen() << " ";
        std::cout << std::endl;
        
        // Recurse into children
        for (Node* d = p->getFirstChild(); d != nullptr; d = d->getNextSibling())
            showNode(d, indent + 3);
        }
}

void Tree::writeTree(Node* p, std::stringstream& strm) {

    if (p == nullptr)
        return;
    
    if (p->getIsTip() == false)
        strm << "(";
    else
        strm << p->getName() << ":" << p->getBrlen();
    
    // LCRS iteration over children
    bool first = true;
    for (Node* d = p->getFirstChild(); d != nullptr; d = d->getNextSibling())
        {
        if (!first)
            strm << ",";
        first = false;
        writeTree(d, strm);
        }
        
    if (p->getIsTip() == false)
        strm << ")" << ":" << p->getBrlen();
}
