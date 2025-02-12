#ifndef __COMMON_H__
#define __COMMON_H__

#include <random>
#include <cstdlib>
#include <vector>

using namespace std;

inline int random(int min, int max, int avoid=-1){
    int num = rand() % (max - min) + min;
    while(num == avoid) num = rand() % (max - min) + min;
    return num;
}


inline double randomFloat(){
    return (double)rand() / RAND_MAX;
}


inline int weight(const vector<int> codeword){
    int result = 0;
    for(unsigned int i = 0; i < codeword.size(); ++i){
        result += codeword[i];
    }
    return result;
}

#endif