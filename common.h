#ifndef __COMMON_H__
#define __COMMON_H__

#include <random>
#include <cstdlib>
#include <vector>

using namespace std;

int random(int min, int max, int avoid=-1){
    int num = rand() % (max - min) + min;
    while(num == avoid) num = rand() % (max - min) + min;
    return num;
}


double randomFloat(){
    return (double)rand() / RAND_MAX;
}

#endif