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



int main(){
    Matrix m(2, 3);
    m.setElement(0, 0, 1);
    m.setElement(0, 2, 1);
    m.setElement(1, 1, 1);

    Matrix transposed_mat = m.transpose();

    cout << (m*transposed_mat).print() << endl;
}
