#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
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
        printf("(%s:%d, %s:%d)\n", type==SYMBOL?"S" : "C", idx, type==SYMBOL?"S" : "C", node->idx);     
        neighbors.push_back(node); 
    }
    const int index() { return idx; }
};

typedef struct{
    vector<vector<Node *> > level_nodes;
    unordered_set<Node *> tree_nodes;
}Subgraph;

template<typename T>
void vectorAppend(vector<T> &src, vector<T> other){
    for(unsigned int i = 0; i < other.size(); ++i){
        src.push_back(other[i]);
    }
}

int expansion(Node *node, vector<Node *> checks){
    bool stop_criteria = false;

    vector<Node *> current_level = {node};
    vector<Node *> next_level;

    unordered_set<Node *> Nset; 
    unordered_set<Node *> NsetComplement(checks.begin(), checks.end());

    unordered_set<Node *> prevNsetComplement;

    int lv = 0;
    while(!stop_criteria){
        ++lv;
        for(unsigned int i = 0; i < current_level.size(); ++i){
            next_level.insert(next_level.end(), 
                current_level[i]->getNeighbors().begin(), 
                current_level[i]->getNeighbors().end()
            );
        }
        
        int prev_N_set_cardinality = Nset.size();
        int prev_N_set_comp_cardinality = NsetComplement.size();
        prevNsetComplement = NsetComplement;
        if(lv % 2 == 1){ // checks
            for(unsigned int i = 0; i < next_level.size(); ++i){
                auto it = NsetComplement.find(next_level[i]);
                if(it != NsetComplement.end()){
                    Nset.insert(*it);
                    NsetComplement.erase(it);
                }
            }
            stop_criteria |= (prev_N_set_cardinality == Nset.size() && Nset.size() < checks.size()); // first stop criteria: N_{s_j}^l stop increasing and |N_{s_j}^l| < m
            stop_criteria |= (prev_N_set_comp_cardinality != 0 && NsetComplement.size() == 0); // second stop criteria: \bar{N_{s_j}^l != 0 && \bar{N_{s_j}^{l+1} != 0}
        }
    }
    // choose one node in prevNsetComplement
    vector<Node *> complement(prevNsetComplement.begin(), prevNsetComplement.end());
    
    sort(complement.begin(), complement.end(), [](Node *n1, Node *n2){
        return n1->degree() < n2->degree();
    });

    // cout << "candidates:" << endl;
    // for(auto node: complement){
    //     cout << "\tchecks" << node->index() << endl;
    // }

    node->addNeighbor(complement[0]);
    complement[0]->addNeighbor(node);
    return 0;
}

void progressiveEdgeGrowthAlgorithm(vector<Node *> symbol_nodes, vector<Node *> check_nodes, vector<int> degrees){

    auto comparator = [](const Node* a, const Node* b) { return a->degree() < b->degree();};
    for(int j = 0; j < symbol_nodes.size(); ++j){
        for(int k = 0; k < degrees[j]; ++k){
            if (k == 0){
                // pick a check node with the lowest degree
                Node* picked_node = *min_element(check_nodes.begin(), check_nodes.end(), comparator); 

                symbol_nodes[j]->addNeighbor(picked_node);
                picked_node->addNeighbor(symbol_nodes[j]);
            }else{
                expansion(symbol_nodes[j], check_nodes);
            }
        }
    }
}

int main(){
    int n = 8, m = 4;
    vector<int> degrees = {2, 2, 2, 2, 2, 2, 2, 2};



    vector<Node *> symbol_nodes;
    vector<Node *> check_nodes;

    for(int i = 0; i < n; ++i) symbol_nodes.push_back(new Node(SYMBOL, i));
    for(int i = 0; i < m; ++i) check_nodes.push_back(new Node(CHECK, i));
    progressiveEdgeGrowthAlgorithm(symbol_nodes, check_nodes, degrees);


    vector<vector<int> > matrix;
    for(int i = 0; i < m; ++i){
        matrix.push_back(vector<int>());
        for(int j = 0; j < n; ++j){
            matrix[i].push_back(0);
        }
    }

    for(unsigned int i = 0; i < symbol_nodes.size(); ++i){
        printf("Symbol[%d]:\n", symbol_nodes[i]->index());
        for(unsigned int j = 0; j < symbol_nodes[i]->degree(); ++j){
            int index = symbol_nodes[i]->getNeighbors()[j]->index();
            printf("\tCheck[%d]\n", index);
            matrix[index][symbol_nodes[i]->index()] = 1;
        }
    }

    for(int i = 0; i < m; ++i){
        for(int j = 0; j < n; ++j){
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}