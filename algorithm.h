#ifndef __ALGORITHM_H__
#define __ALGORITHM_H__

#include <iostream>
#include <vector>

#include "matrix.h"
#include "genetic.h"

double sparsity(const Matrix &mat);
vector<vector<int> > codewords(const Matrix &mat);
int minimumDistance(const vector<vector<int> > &codewords);
int minimumDistance(const Matrix &parity_check_mat);

class Chromosome : public ChromosomeBase{
protected:
    int ** mat;
    int v, h;
    double a, b, r, d;

    Matrix matrix;
public:
    Chromosome():ChromosomeBase(), mat(nullptr), v(0), h(0){}
    Chromosome(int v, int h, Initializer init=RANDOM);
    Chromosome(const Chromosome &);
    virtual ~Chromosome();
    virtual double computeFitnessValue();

    inline Matrix matrixForm(){return matrix;}
};

class GeneticAlgorithm{
private:
    int v, h;
    double crossover_rate;
    double mutation_rate;
    double elite_rate;
    double roulete_rate;

    int chromosome_size;
    int pop_size;

    vector<Chromosome *> available_offspring;
    vector<Chromosome *> population;


    vector<Chromosome *> crossover();
    vector<Chromosome *> mutation();

    vector<Chromosome *> eliteSelection();
    vector<Chromosome *> rouleteSelection();

    vector<Chromosome *> offspringRecycle(int num);
public:
    GeneticAlgorithm();
    GeneticAlgorithm(
        int initial_population_size,
        int v, int h,
        double c_rate, 
        double m_rate, 
        double e_rate, 
        double r_rate
    );
    ~GeneticAlgorithm(){
        for (auto *c: population){
            delete c;
        }
        population.clear();
    }

    Chromosome run(int iterations);
};


#endif