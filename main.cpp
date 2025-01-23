#include <iostream>
#include <vector>
#include <string>

#include "generator.h"
#include "expression.h"

using namespace std;

// Inputs: 
//  Stabilizer code [[n, k]] and its generator set <S1, S2, ..., Sn-k> 
//  r



int main(){
    Generator g1(1, 1, 2, 4);
    Generator g2(1, 1, 2, 4);

    Expression exp(2, 4);
    exp = g1 + g2;

    cout << exp.print() << endl;
}
