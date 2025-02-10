#ifndef __GENETIC_ALGORITHM_H__
#define __GENETIC_ALGORITHM_H__

#include "matrix.h"
#include <vector>
#include <string>

using namespace std;

typedef enum{
    ALL_ZERO, RANDOM
}Initializer;

class ChromosomeBase{
protected:
    int size;
    int *genes;
    double fitness_val;
    void randomizeGenes();
    virtual void copyFrom(const ChromosomeBase &other);
public:
    ChromosomeBase():size(0){genes = NULL;}
    ChromosomeBase(int size, Initializer init=RANDOM);
    ChromosomeBase(const ChromosomeBase &);
    virtual ~ChromosomeBase();
    ChromosomeBase & operator=(const ChromosomeBase &);


    static void crossover(const ChromosomeBase &, const ChromosomeBase &, ChromosomeBase &, ChromosomeBase&);
    static void mutation(const ChromosomeBase &, ChromosomeBase &);
    string print();
    inline double fitnessValue() const { return fitness_val;}
    virtual double computeFitnessValue()=0;
};



#endif