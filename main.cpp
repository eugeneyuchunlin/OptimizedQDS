#include <iostream>
#include <vector>
#include <string>

#include "generator.h"
#include "expression.h"
#include "matrix.h"
#include "algorithm.h"

using namespace std;

// Inputs: 
//  Stabilizer code [[n, k]] and its generator set <S1, S2, ..., Sn-k> 
//  r

void testMatrix(vector<vector<int> > input){
    Matrix m(input);
    Matrix rref_mat = rref(m);
    Matrix null_space = nullSpace(rref_mat);
    cout << "null space : \n"<< null_space.print() << endl;
    Matrix null_space2 = nullSpace(rref(null_space.transpose()));
    
    cout << "null space 2 : \n"<< (null_space2).print() << endl;
}


int main(){
    // testMatrix({{1, 1, 0, 1, 0, 0, 1},{1, 0, 1, 0, 1, 0, 1}, {0, 1, 1, 0, 0, 1, 1}});
    // testMatrix({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    // testMatrix({{0, 0, 0}, {0, 0, 0}, {0, 0, 0}});
    // testMatrix({{1, 1, 1}});
    // testMatrix({{1, 0}, {0, 1}, {0, 0}});
    // testMatrix({{1, 0, 1}, {0, 1, 1}});
    // testMatrix({{1, 1, 0}, {0, 1, 1}, {1, 0, 1}});
    // testMatrix({{1, 0, 1}, {0, 1, 1}, {1, 1, 0}});
    // testMatrix({{1}, {0}, {1}});
    // testMatrix({{0}, {0}, {0}});

    // Matrix m({{1, 1, 0, 1, 0, 0, 1},{1, 0, 1, 0, 1, 0, 1}, {0, 1, 1, 0, 0, 1, 1}});
    // Matrix m({{1, 0, 0, 0, 1, 1, 0}, {0, 1, 0, 0, 1, 0, 1}, {0, 0, 1, 0, 0, 1, 1}, {0, 0, 0, 1, 1, 1, 1}});
    Matrix m({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    double sparsity_val = sparsity(m);
    Matrix rref_mat = rref(m);

    cout << "RREF: \n" << rref_mat.print() << endl;

    double new_sparsity_val = sparsity(rref_mat);
    cout << "Sparsity: " << sparsity_val << endl;
    cout << "New Sparsity: " << new_sparsity_val << endl;

    Matrix codewords_mat(codewords(m));
    cout << "Codewords: \n" << codewords_mat.print() << endl;

    int minimum_d = minimumDistance(codewords(m));
    cout << "minimum distance: " << minimum_d << endl;

    cout << "matrix: \n" << m.print() << endl;
    cout << "couting depth: " << countingDepth(m) << endl;
    cout << "correction depth: " << correctionDepth(m) << endl; 
    // cout << (sizeof(int) <<3) << endl;

    return 0;
}
