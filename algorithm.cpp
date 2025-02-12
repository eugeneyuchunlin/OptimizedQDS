#include "algorithm.h"
#include "matrix.h"
#include "common.h"

#include <iostream>
#include <vector>

using namespace std;

double sparsity(const Matrix &mat){
    int count = 0;
    for(int i = 0; i < mat.size_v(); ++i){
        for(int j = 0; j < mat.size_h(); ++j){
            if(mat.element(i, j) != 0){
                count++;
            }
        }
    }
    return (double)count / (mat.size_v() * mat.size_h());
}

vector<vector<int> > codewords(const Matrix &mat){
    vector<vector<int> > codewords;
    int v = mat.size_v();
    int h = mat.size_h();
    int num_codewords = 1 << h;

    vector<int> num_representation;
    // cout << mat.print() << endl;
    for(int i = 0; i < v; ++i){
        int num = 0;
        for(int j = 0; j < h; ++j){
            num += mat.mat[i][j] << (h-j-1);
        }  
        // cout << num << endl;
        num_representation.push_back(num);
    }

    for(int i = 0; i < num_codewords; ++i){
        vector<int> codeword;
        codeword.assign(v, 0);

        for(int k = 0; k < v; ++k) {
            codeword[k] = num_representation[k] & i;
        }
        for(int k = 0; k < v; ++k) {
            int num = 0;
            for(int j = 0; j < h; ++j){
                num += (codeword[k] >> j) & 1;
            }
            codeword[k] = num & 1;
        }

        // cout << i << ": ";
        // for(int k = 0; k < v; ++k) cout << codeword[k] << " ";
        // cout << endl;

        codewords.push_back(codeword);
    }

    return codewords;
}


int minimumDistance(const vector<vector<int> > &codewords){
    int mini_d = codewords[0].size();
    for(unsigned int i = 0; i < codewords.size(); ++i){
        int w = weight(codewords[i]);
        if(w < mini_d && w != 0){
            mini_d = w;
        }
    }
    return mini_d;
}

int minimumDistance(const Matrix &parity_check_matrix){
    Matrix generator_mat = nullSpace(rref(parity_check_matrix));
    vector<vector<int> > cwds = codewords(generator_mat);
    Matrix cwds_mat(cwds);
    // cout << generator_mat.print() << endl;
    // cout << cwds_mat.print() << endl;
    return minimumDistance(cwds);
}

inline int columnDegree(const vector<vector<int> > &mat, int col){
    int deg = 0;
    for(unsigned int i = 0; i < mat.size(); ++i){
        deg += mat[i][col];
    }
    return deg;
}

vector<int> degree(const Matrix &mat){
    int h = mat.size_h();

    vector<int> degrees;
    degrees.assign(h, 0);

    for(int i = 0; i < h; ++i){
        degrees[i] = columnDegree(mat.mat, i);
    }

    return degrees;
}


// this implies that we want to have 
int countingDepth(const Matrix &mat){
    
    int h = mat.size_h();

    int depth = 0;

    for(int i = 0; i < h; ++i){
        int deg = columnDegree(mat.mat, i);
        
        // cout << "deg " << i << " " << deg << " ";

        int single_depth = 0;
        for(int j = 1; j <= deg; ++j){
            // floor(log2(j)) + 1
            single_depth += ((sizeof(j)<<3) - __builtin_clz(j));
        }
        // cout << "single depth " << single_depth << endl;
        depth += single_depth;
    }

    return depth;
}

int correctionDepth(const Matrix &mat){

    int v = mat.size_v();
    int h = mat.size_h();

    vector<int> degrees = degree(mat);
    vector<int> stat_degree;
    stat_degree.assign(v+1, 0);
    for(int i = 0; i < h; ++i){
        ++stat_degree[degrees[i]];
    }

    // for(int i = 0; i < v+1; ++i){
    //     cout << i << ": " << stat_degree[i] << endl;
    // }

    int depth = 0;
    for(int i = 1; i < v+1; ++i){
        for(int j = i; j < v+1; ++j){
            depth += stat_degree[j];
        }
    }

    return depth;
}

Chromosome::Chromosome(const Chromosome & other): ChromosomeBase(other), h(other.h), a(other.a), b(other.b), r(other.r), d(other.d){
    mat = new int*[other.v];
    v = other.v;
    for(int i = 0; i < v; ++i) mat[i] = genes + i*h;
}

Chromosome::Chromosome(int v, int h, Initializer init):
    ChromosomeBase(v*h, init), 
    mat(NULL), v(v), h(h), matrix(v, h)
{
    mat = new int*[v];
    for(int i = 0; i < v; ++i){
        mat[i] = genes + h*i;
    } 

    a = 0.2;
    b = 0.3;
    r = 0.2;
    d = 0.3;
}



Chromosome::~Chromosome(){
    delete[] mat;
    mat = nullptr;
}


void Chromosome::updateMatrix(){
    for(int i = 0 ; i < v; ++i){
        for(int j = 0; j < h; ++j){
            matrix.setElement(i, j, mat[i][j]);
        }
    }
}

double Chromosome::computeFitnessValue(){
    // convert genes to a matrix
    updateMatrix(); 

    a = 0.2;
    b = 1;
    r = 5;
    d = 40;
    fitness_val = 0;
    // fit += a*
    Matrix optimizedMat = matrixOptimization(matrix);
    // Matrix optimizedMat = matrix;
    fitness_val += b*countingDepth(optimizedMat);
    fitness_val += r*correctionDepth(optimizedMat);
    fitness_val += d*(1/log10(minimumDistance(optimizedMat) + 1e-10));

    return fitness_val;
}

vector<int> Chromosome::twoDimensionCycleRand(int lb, int ub){
    int rnd1 = random(lb, ub);
    int rnd2 = random(lb, ub, rnd1);

    rnd2 += ub*(rnd2 < rnd1);
    vector<int> nums;
    for(int i = rnd1; i < rnd2; ++i) nums.push_back(i % ub);

    return nums;
}

void Chromosome::crossover(const Chromosome & p1, const Chromosome & p2, Chromosome & o1, Chromosome & o2){
    o1 = p1;
    o2 = p2;

    int v = p1.v, h = p1.h;

    vector<int> rows = twoDimensionCycleRand(0, v);
    vector<int> cols = twoDimensionCycleRand(0, h);

    // cout <<"new crossover is called" << endl;
    for(unsigned int i = 0; i < rows.size(); ++i){
        for(unsigned int j = 0; j < cols.size(); ++j){
            o2.mat[rows[i]][cols[j]] = p1.mat[rows[i]][cols[j]];
            o1.mat[rows[i]][cols[j]] = p2.mat[rows[i]][cols[j]];
        }
    }
}

void Chromosome::mutation(const Chromosome &p, Chromosome & o){
    int v = p.v, h = p.h;
    o = p;

    // cout <<"new mutation is called" << endl;
    vector<int> rows = twoDimensionCycleRand(0, v);
    vector<int> cols = twoDimensionCycleRand(0, h);

    for(unsigned int i = 0; i < rows.size(); ++i){
        for(unsigned int j = 0; j < cols.size(); ++j){
            o.mat[rows[i]][cols[j]] = p.mat[rows[i]][cols[j]] ^1; // flip the bit
        }
    }
}


Chromosome & Chromosome::operator=(const Chromosome &other){
    // cout << "Chromosome operator= is called" << endl;
    if (this == &other) return *this;

    ChromosomeBase::operator=(other);

    if(other.v != this->v){
        if(this->mat != nullptr){
            delete[] this->mat;
        } 
        this->mat = new int*[this->v];
        for(int i = 0; i < v; ++i){
            this->mat[i] = genes + h*i;
        }
        this->v = other.v;
        this->h = other.h;
    }

    return *this;
}


GeneticAlgorithm::GeneticAlgorithm(){
    pop_size = 0;
    chromosome_size = 0;
    crossover_rate = mutation_rate = 0;
    elite_rate = roulete_rate = 0;
}

GeneticAlgorithm::GeneticAlgorithm(
    int initial_population_size, 
    int v, int h,
    double c_rate,
    double m_rate, 
    double e_rate, 
    double r_rate):
    v(v), h(h),
    crossover_rate(c_rate),
    mutation_rate(m_rate),
    elite_rate(e_rate),
    roulete_rate(r_rate),
    chromosome_size(v*h),
    pop_size(initial_population_size)
{

    
    
    for(int i = 0 ; i < initial_population_size; ++i){
        population.push_back(new Chromosome(v, h));
    }
}

vector<Chromosome *> GeneticAlgorithm::offspringRecycle(int num){


    vector<Chromosome *> ready_to_use;
    int j = available_offspring.size() - 1;
    for(unsigned int i = 0; i < num && i < available_offspring.size(); ++i, --j){
        ready_to_use.push_back(available_offspring[j]);
        available_offspring.pop_back();
    }

    for(unsigned int i = ready_to_use.size(); i < num; ++i){
        ready_to_use.push_back(new Chromosome(v, h, RANDOM)); 
    }

    return ready_to_use;
}

vector<Chromosome *> GeneticAlgorithm::crossover(){
    int num_offspring = (int)(pop_size * crossover_rate);
    int times = num_offspring >> 1;

    vector<Chromosome *> ready_to_use_offspring = offspringRecycle(num_offspring);

    vector<Chromosome *> new_offspring;
    int current_pop_size = population.size();
    for(int i = 0, j = 0; i < times; ++i, j+=2){
        int n_p1 = random(0, current_pop_size);
        int n_p2 = random(0, current_pop_size, n_p1);

        Chromosome * o1 = ready_to_use_offspring[j];
        Chromosome * o2 = ready_to_use_offspring[j+1];
        Chromosome::crossover(*population[n_p1], *population[n_p2], *o1, *o2);
        new_offspring.push_back(o1);
        new_offspring.push_back(o2);
    } 

    return new_offspring;
}

vector<Chromosome *> GeneticAlgorithm::mutation(){
    int num_offspring = (int)(pop_size * mutation_rate);

    vector<Chromosome *> ready_to_use = offspringRecycle(num_offspring);

    vector<Chromosome *> new_offspring;
    int current_pop_size = population.size();
    for(int i = 0; i < num_offspring; ++i){
        int p = random(0, current_pop_size);

        Chromosome * o = ready_to_use[i];
        Chromosome::mutation(*population[p], *o);
        new_offspring.push_back(o); 
    }

    return new_offspring;
}

vector<Chromosome *> GeneticAlgorithm::eliteSelection(){
    vector<Chromosome *> elites; 
    int num_elites = pop_size * elite_rate;

    for(unsigned int i = 0, j = population.size() - 1; i < num_elites; ++i, --j){
        elites.push_back(population[j]);
        population.pop_back();
    }

    return elites;
}

vector<Chromosome *> GeneticAlgorithm::rouleteSelection(){
    typedef struct{
        Chromosome * c;
        double val;
    }Fitnessbundle;

    vector<Fitnessbundle> bundles;

    double normalization = 0, cumulative = 0;
    for(unsigned int i = 0; i < population.size(); ++i){
        normalization += 1/population[i]->fitnessValue();
        bundles.push_back({
            .c = population[i],
            .val = 1/population[i]->fitnessValue()
        });
    }

    for(unsigned int i = 0; i < bundles.size(); ++i){
        bundles[i].val /= normalization;
    }

    for(unsigned int i = 0; i < bundles.size(); ++i){
        cumulative = bundles[i].val += cumulative;
        // cout << i << ": " << bundles[i].val << endl;
    }

    vector<int> duplicacy;

    duplicacy.assign(bundles.size(), 0);

    vector<double> rnds;
    int num_selection = pop_size * roulete_rate;
    int failed_counts = num_selection;
    // TODO: Should optimize the following segment
    for(int i = 0; i < failed_counts; ++i) rnds.push_back(randomFloat());
    sort(rnds.begin(), rnds.end());

    int j = 0;
    for(int i = 0; i < failed_counts; ++i){
        while(j < bundles.size() && rnds[i] >= bundles[j].val) ++j;
        if(j == bundles.size()) j = bundles.size()-1;

        if(duplicacy[j] == 0) {
            duplicacy[j] = 1;
            failed_counts -=1;
        }
    }
    // cout << "failed counts: " << failed_counts << endl;
    rnds.clear();
    
    // from duplicacy[j] == 0, pick failed_counts out
    while(failed_counts > 0){
        int rnd;
        do{
            rnd = random(0, duplicacy.size());
        }while(duplicacy[rnd] == 1);
        duplicacy[rnd] = 1;
        --failed_counts;
    }

    vector<Chromosome *> next_gen;
    for(int i = 0; i < duplicacy.size(); ++i){
        if(duplicacy[i]) next_gen.push_back(bundles[i].c);
        else available_offspring.push_back(bundles[i].c);
    }
    return next_gen;
}


Chromosome GeneticAlgorithm::run(int iterations){

    // for(int i = 0; i < pop_size; ++i){
    //     double f_val = population[i]->computeFitnessValue();
    //     cout << i << ":" << population[i]->fitnessValue() << endl;
    // }


    for(int i = 0; i < iterations; ++i){
        // crossover
        vector<Chromosome *> crsv_ofsprg = crossover();
        vector<Chromosome *> mut_ofsprg = mutation();

        for(int j  = 0 ; j < crsv_ofsprg.size(); ++j) population.push_back(crsv_ofsprg[j]);
        for(int j  = 0 ; j < mut_ofsprg.size(); ++j) population.push_back(mut_ofsprg[j]);
        for(int j = 0; j < population.size(); ++j) population[j]->computeFitnessValue();

        sort(population.begin(), population.end(), [](const Chromosome *a, const Chromosome *b){
            return a->fitnessValue() > b->fitnessValue();
        });

        // for(int j = 0; j < population.size(); ++j){
        //     cout << j << ": " << population[j]->fitnessValue() << endl;
        // }

        // selection
        vector<Chromosome *> elites = eliteSelection();
        // for(int j = 0; j < population.size(); ++j){
        //     cout << j << ": " << population[j]->fitnessValue() << endl;
        // }

        // cout <<"============" << endl;
        // for(int j = 0; j < elites.size(); ++j){
        //     cout << j << ": " << elites[j]->fitnessValue() << endl;
        // }
        
        vector<Chromosome *> roulette_result = rouleteSelection();
        population.clear();
        for(int j = 0; j < elites.size(); ++j) population.push_back(elites[j]);
        for(int j = 0; j < roulette_result.size(); ++j) population.push_back(roulette_result[j]);

        // cout << "population size: "<< population.size() << endl;
        // cout << "available ofs: " << available_offspring.size()  << endl;
        cout << i << ": " <<  elites[0]->fitnessValue() << endl;
    }
    Matrix original = population[0]->matrixForm();
    Matrix matrix = matrixOptimization(original);
    cout << "original:\n" << original.print() << endl;
    cout << matrix.print() << endl;
    cout << rref(matrix).print() << endl;
    int distance = minimumDistance(matrix);
    cout << "minimum distance : " << distance << endl;
    cout << "counting depth: " << countingDepth(matrix)<< endl;
    cout << "correction depth: " << correctionDepth(matrix) << endl;

    Matrix generator_matrix = nullSpace(rref(matrix));
    cout << generator_matrix.print() << endl;
    return *population[0];
}


