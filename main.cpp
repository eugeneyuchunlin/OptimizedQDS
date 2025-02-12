#include <iostream>
#include <vector>
#include <string>

#include "generator.h"
#include "expression.h"
#include "matrix.h"
#include "algorithm.h"
#include "genetic.h"

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
    // Matrix m({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    // double sparsity_val = sparsity(m);
    // Matrix rref_mat = rref(m);

    // cout << "RREF: \n" << rref_mat.print() << endl;

    // double new_sparsity_val = sparsity(rref_mat);
    // cout << "Sparsity: " << sparsity_val << endl;
    // cout << "New Sparsity: " << new_sparsity_val << endl;

    // Matrix codewords_mat(codewords(m));
    // cout << "Codewords: \n" << codewords_mat.print() << endl;

    // int minimum_d = minimumDistance(codewords(m));
    // cout << "minimum distance: " << minimum_d << endl;

    // cout << "matrix: \n" << m.print() << endl;
    // cout << "couting depth: " << countingDepth(m) << endl;
    // cout << "correction depth: " << correctionDepth(m) << endl; 
    // cout << (sizeof(int) <<3) << endl;

    // srand(time(NULL));

    GeneticAlgorithm ga(150, 7, 13, 0.8, 0.2, 0.2, 0.8);
    ga.run(2000);
    // Matrix pmat({{1, 0, 0, 0, 0, 1, 1}, {0, 0, 1, 0, 0, 0, 0}, {1, 1, 0, 0, 1, 0, 0}, {0, 1, 0, 1, 0, 0, 1}});
    // Matrix pmat2 = matrixOptimization(pmat);
    // int dist1 = minimumDistance(pmat);
    // int dist2 = minimumDistance(pmat2);
    // cout << "minimum distance1: " << dist1 << endl;
    // cout << "minimum distance2: " << dist2 << endl;

    // vector<vector<int> > codeword1 = codewords(nullSpace(rref(pmat)));
    // vector<vector<int> > codeword2 = codewords(nullSpace(rref(pmat2)));

    // Matrix cwd_mat1(codeword1);
    // Matrix cwd_mat2(codeword2);

    // cout << cwd_mat1.print() << endl;
    // cout << cwd_mat2.print() << endl;
    

    // printf("countingDepth(p1, p2) = (%d, %d)\n", countingDepth(pmat), countingDepth(pmat2));
    // printf("CorrectionDepth(p1, p2) = (%d, %d)\n", correctionDepth(pmat), correctionDepth(pmat2));

    // cout << "pmat 1 : \n" << rref(pmat).print() << endl;
    // cout << "pmat 2 : \n" << rref(pmat2).print() << endl;

    return 0;
}
