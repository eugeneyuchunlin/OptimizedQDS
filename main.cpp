#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <csignal>

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

void outputMatrix(const std::string& filename, Matrix& matrix) {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    for (unsigned int i = 0; i < matrix.size_v(); ++i) {
        for (unsigned int j = 0; j < matrix.size_h(); ++j) {
            file << matrix.element(i, j);
            if (j < matrix.size_h() - 1) {
                file << ","; 
            }
        }
        file << "\n";
    }

    file.close();
}

int seed;
void signal_handler(int signal_num) 
{ 
    cout << "The interrupt signal is (" << signal_num 
         << "). \n"; 
    cout << "seed: " << seed << endl;
    exit(signal_num); 
} 


int main(int argc, const char * argv[]){
    signal(SIGINT, signal_handler);

    seed = time(NULL);
    srand(seed);
    if(argc < 10){
        printf("Usage: [executable file] [stabilizers] [v] [h] [alpha] [beta] [gamma] [delta] [zeta] [iterations]");
        exit(-1);
    }

    csv_t qcode_csv;
    qcode_csv.read(argv[1], "r", false);
    vector<vector<string> > data = qcode_csv.getData();
    int quantum_v = data.size(), quantum_h = 0;
    if(quantum_v){
        quantum_h = data[0].size();
    }
    Matrix quantum_code(quantum_v, quantum_h);
    for(int i = 0; i < data.size(); ++i){
        for(int j = 0; j < data[i].size(); ++j){
            if(data[i][j] != "I"){
                quantum_code.setElement(i, j, 1);
            }
        }
    }

    // cout << quantum_code.print() << endl;
    int v = atoi(argv[2]), h = atoi(argv[3]);
    Parameters params = {
        .alpha = atof(argv[4]),
        .beta = atof(argv[5]),
        .gamma = atof(argv[6]),
        .delta = atof(argv[7]),
        .zeta = atof(argv[8])
    };
    int iterations = atoi(argv[9]);


    vector<map<string, string> > iteration_data;
// rerun:
    GeneticAlgorithm ga(quantum_code, 200, v, h, 0.8, 0.2, 0.2, 0.8, params, {"./experiments/4_10/3/pmatrix_4_10.csv"});
    Chromosome best_result = ga.run(iterations, iteration_data);

    csv_t csv;
    for(unsigned int i = 0; i < iteration_data.size(); ++i){
        csv.addData(iteration_data[i]);
    }
    string prefix = "./output/";
    string filename = "output_" + to_string(v) + "_" + to_string(h) + ".csv";
    csv.write(prefix + filename, "w");

    best_result.updateMatrix();
    Matrix matrix(best_result.matrixForm());

     
    string pmatrix_csv_filename = "pmatrix_" + to_string(v) + "_" + to_string(h) + ".csv";
    outputMatrix(prefix + pmatrix_csv_filename, best_result.optimized_matrix);


    string gmatrix_csv_filename = "gmatrix_" + to_string(v) + "_" + to_string(h) + ".csv";
    Matrix gmatrix = nullSpace(best_result.optimized_matrix);
    outputMatrix(prefix + gmatrix_csv_filename, gmatrix);


    cout << best_result.optimized_matrix.print() << endl;
    cout << gmatrix.print() << endl;

    cout << "counting depth: " << countingDepth(best_result.optimized_matrix) << endl;
    cout << "correction depth: " << correctionDepth(best_result.optimized_matrix) << endl;
    cout << "encoding cost: " << encodingDepth(gmatrix, rowWeights(quantum_code)) << endl;
    cout << "minimum distance: " << minimumDistance(best_result.optimized_matrix)  << endl;
    int gr = girth(best_result.optimized_matrix);
    cout << "girth: " << girth(best_result.optimized_matrix) << endl;

    // if(gr < 7){
    //     goto rerun;
    // }

  return 0;
}