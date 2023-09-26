#ifndef Tree_hpp
#define Tree_hpp

#include <string>
class Node;



class Tree {

    public:
                            Tree(void) = delete;
                            Tree(const Tree& t);
                            Tree(std::string fileName);
                           ~Tree(void);
        Tree&               operator=(const Tree& rhs);
        std::vector<Node*>& getDownPassSequence(void) { return downPassSequence; }
        std::vector<Node*>& getInteriorDownPassSequence(void) { return interiorDownPassSequence; }
        int                 getNumnodes(void) { return numNodes; }
        Node*               getRoot(void) { return root; }
        void                print(void);
        void                setNumAreas(int x) { numAreas = x; }
    
    private:
        Node*               addNode(void);
        void                clone(const Tree& t);
        void                deleteNodes(void);
        void                initializeDownPassSequence(void);
        void                passDown(Node* p);
        void                readAhead(std::string& token, std::string& newickStr, int& i);
        void                showNode(Node* p, int indent);
        Node*               root;
        int                 numNodes;
        int                 numTips;
        int                 numAreas;
        Node*               endDownPassSequence;
        std::vector<Node*>  nodes;
        std::vector<Node*>  downPassSequence;
        std::vector<Node*>  interiorDownPassSequence;
};

#endif
