#include <iostream>
#include <vector>
#include <string>

#include "generator.h"
#include "expression.h"
#include "matrix.h"
#include "algorithm.h"
#include "genetic.h"
#include "csv.h"

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


int main(int argc, const char * argv[]){
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
    if(argc < 3){
        printf("Usage: [executable file] [v] [h]");
        exit(-1);
    }

    vector<map<string, string> > iteration_data;
    GeneticAlgorithm ga(150, atoi(argv[1]), atoi(argv[2]), 0.8, 0.2, 0.2, 0.8);
    Chromosome best_result = ga.run(500, iteration_data);

    csv_t csv;
    for(unsigned int i = 0; i < iteration_data.size(); ++i){
        csv.addData(iteration_data[i]);
    }

    string filename = "output_" + string(argv[1]) + "_" + string(argv[2]) + ".csv";
    csv.write(filename, "w");

    best_result.updateMatrix();
    Matrix matrix(best_result.matrixForm());

    cout << matrix.print() << endl;
    cout << best_result.optimized_matrix.print() << endl;
    return 0;
}
