#include <iostream>
#include <vector>
#include <string>

#include "generator.h"
#include "expression.h"
#include "matrix.h"

using namespace std;

// Inputs: 
//  Stabilizer code [[n, k]] and its generator set <S1, S2, ..., Sn-k> 
//  r

void testMatrix(vector<vector<int> > input){
    Matrix m(input);
    Matrix rref_mat = rref(m);
    Matrix null_space = nullSpace(rref_mat);
    cout << "null space : \n"<< null_space.print() << endl;
}


int main(){
    testMatrix({{1, 1, 0, 1, 0, 0, 1},{1, 0, 1, 0, 1, 0, 1}, {0, 1, 1, 0, 0, 1, 1}});
    testMatrix({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    testMatrix({{0, 0, 0}, {0, 0, 0}, {0, 0, 0}});
    testMatrix({{1, 1, 1}});
    testMatrix({{1, 0}, {0, 1}, {0, 0}});
    testMatrix({{1, 0, 1}, {0, 1, 1}});
    testMatrix({{1, 1, 0}, {0, 1, 1}, {1, 0, 1}});
    testMatrix({{1, 0, 1}, {0, 1, 1}, {1, 1, 0}});
    testMatrix({{1}, {0}, {1}});
    testMatrix({{0}, {0}, {0}});
}
