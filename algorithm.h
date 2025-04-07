#ifndef __ALGORITHM_H__
#define __ALGORITHM_H__

#include <iostream>
#include <vector>
#include <map>

#include "matrix.h"
#include "genetic.h"

double sparsity(const Matrix &mat);
vector<vector<int> > codewords(const Matrix &mat);
int minimumDistance(const vector<vector<int> > &codewords);
int minimumDistance(const Matrix &parity_check_mat);
int encodingDepth(const Matrix &G, const vector<int> &weights);
vector<int> rowWeights(const Matrix &stabilizer);

typedef struct{
    double alpha;
    double beta;
    double gamma;
    double delta;
    double zeta;
}Parameters;

class Chromosome : public ChromosomeBase{
private:
    static vector<int> twoDimensionCycleRand(int lb, int ub);
protected:
    int ** mat;
    int v, h;

    Matrix matrix;
    vector<int> *stabilizer_weights;
public:
    Chromosome():ChromosomeBase(), mat(nullptr), v(0), h(0), params(nullptr){}
    Chromosome(vector<int> & stabilizer_weights, Parameters & params, int v, int h, Initializer init=RANDOM);
    Chromosome(string filename, vector<int> & stabilizer_weights, Parameters & params, int v, int h, Initializer init=RANDOM);
    Chromosome(const Chromosome &);

    virtual ~Chromosome();
    virtual double computeFitnessValue();

    void updateMatrix();
    inline Matrix matrixForm(){return matrix;}

    static void crossover(const Chromosome &, const Chromosome &, Chromosome &, Chromosome &);
    static void mutation(const Chromosome &, Chromosome &);

    Chromosome & operator=(const Chromosome & other);
    Matrix optimized_matrix;
    Matrix generator_matrix;

    Parameters * params;
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
    
    vector<int> stabilizer_weights; 
    Parameters params;
public:
    GeneticAlgorithm();
    GeneticAlgorithm(
        Matrix stabilizers,
        int initial_population_size,
        int v, int h,
        double c_rate, 
        double m_rate, 
        double e_rate, 
        double r_rate,
        Parameters params,
        vector<string> injection_population
    );
    ~GeneticAlgorithm(){
        for (auto *c: population){
            delete c;
        }
        population.clear();
    }

    void populationInjection(vector<Chromosome *> pop);
    Chromosome run(int iterations, vector<map<string, string> > & iteration_data);
};


#endif