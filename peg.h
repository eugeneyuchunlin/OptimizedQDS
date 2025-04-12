#ifndef __PEG_H__
#define __PEG_H__

#include <vector>
using namespace std;


enum NodeType{
    SYMBOL, CHECK
};

class Node{
private:
    vector<Node *> neighbors; 
    const int idx;
public:

    const NodeType type;
    Node(NodeType type, int idx) : type(type), idx(idx) {}
    int degree() const { return neighbors.size(); }
    vector<Node *> getNeighbors() {return neighbors; }
    void addNeighbor(Node* node) { 
        // printf("(%s:%d, %s:%d)\n", type==SYMBOL?"S" : "C", idx, type==SYMBOL?"S" : "C", node->idx);     
        neighbors.push_back(node); 
    }
    const int index() { return idx; }
};


class PEG{
protected:
    int _m, _n;
    vector<int> _degree;
    vector<Node *> _symbol_nodes, _check_nodes;

    virtual void connect(Node * symbol, Node *check);
    virtual Node * pickup(Node * symbol, vector<Node *>);
    int expansion(Node *node, vector<Node *> checks);

    virtual void initializeNodes();
public:
    void algorithm();
    PEG(int n, int m, vector<int> degree);
    PEG(int n, int m, int degree);
};

#endif