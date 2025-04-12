#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include "peg.h"


using namespace std;

template<typename T>
void vectorExtend(vector<T> & src, vector<T> other){
    for(auto obj : other){
        src.push_back(obj);
    }
}


int PEG::expansion(Node *node, vector<Node *> checks){
    bool stop_criteria = false;

    vector<Node *> current_level;
    vector<Node *> next_level = {node};

    unordered_set<Node *> Nset; 
    unordered_set<Node *> NsetComplement(checks.begin(), checks.end());

    unordered_set<Node *> prevNsetComplement;

    int lv = 0;
    while(!stop_criteria){
        current_level = next_level;
        next_level.clear();

        ++lv;
        for(unsigned int i = 0; i < current_level.size(); ++i){
            vectorExtend(next_level, current_level[i]->getNeighbors());
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
    connect(node, pickup(node, complement));
    // node->addNeighbor(pickup(node, complement));
    // complement[0]->addNeighbor(node);
    return 0;
}

Node * PEG::pickup(Node * symbol, vector<Node *> nodes){
    sort(nodes.begin(), nodes.end(), [](Node *n1, Node *n2){
        if(n1->degree() < n2->degree()){
            return true;
        }else if(n1->degree() == n2->degree()){
            return n1->index() < n2->index();
        }else{
            return false;
        }
    });
    return nodes[0];
}

void PEG::connect(Node * symbol, Node * check){
    symbol->addNeighbor(check);
    check->addNeighbor(symbol);
}

void PEG::algorithm(){

    auto comparator = [](const Node* a, const Node* b) { return a->degree() < b->degree();};
    for(int j = 0; j < _symbol_nodes.size(); ++j){
        // cout <<"symbol: " << _symbol_nodes[j]->index() << endl;
        for(int k = 0; k < _degree[j]; ++k){
            if (k == 0){
                // pick a check node with the lowest degree
                Node* picked_node = PEG::pickup(_symbol_nodes[j], _check_nodes);
                connect(_symbol_nodes[j], picked_node);
            }else{
                expansion(_symbol_nodes[j], _check_nodes);
            }
        }
    }
}

void PEG::initializeNodes(){
    for(int i = 0; i < _m; ++i) _symbol_nodes.push_back(new Node(SYMBOL, i));
    for(int i = 0; i < _n; ++i) _check_nodes.push_back(new Node(CHECK, i));
}

PEG::PEG(int n, int m, vector<int> degree): _m(m), _n(n), _degree(degree){
    initializeNodes();
}

PEG::PEG(int n, int m, int degree): _m(m), _n(n), _degree(m, degree){
    initializeNodes();
}

// int main(){
//     int n = 8, m = 4;
//     vector<int> degrees = {2, 2, 2, 2, 2, 2, 2, 2};



//     vector<Node *> symbol_nodes;
//     vector<Node *> check_nodes;

//     for(int i = 0; i < n; ++i) symbol_nodes.push_back(new Node(SYMBOL, i));
//     for(int i = 0; i < m; ++i) check_nodes.push_back(new Node(CHECK, i));
//     progressiveEdgeGrowthAlgorithm(symbol_nodes, check_nodes, degrees);


//     vector<vector<int> > matrix;
//     for(int i = 0; i < m; ++i){
//         matrix.push_back(vector<int>());
//         for(int j = 0; j < n; ++j){
//             matrix[i].push_back(0);
//         }
//     }

//     for(unsigned int i = 0; i < symbol_nodes.size(); ++i){
//         printf("Symbol[%d]:\n", symbol_nodes[i]->index());
//         for(unsigned int j = 0; j < symbol_nodes[i]->degree(); ++j){
//             int index = symbol_nodes[i]->getNeighbors()[j]->index();
//             printf("\tCheck[%d]\n", index);
//             matrix[index][symbol_nodes[i]->index()] = 1;
//         }
//     }

//     for(int i = 0; i < m; ++i){
//         for(int j = 0; j < n; ++j){
//             cout << matrix[i][j] << " ";
//         }
//         cout << endl;
//     }
// }