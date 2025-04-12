#include "algorithm.h"
#include "matrix.h"
#include "common.h"

#include <utility>
#include <fstream>
#include <sstream>
#include <map>
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

int encodingDepth(const Matrix &G, const vector<int> &weights){
    int h = G.size_h(), v = G.size_v();
    if(h < weights.size()){
        // size incorrect
        throw std::runtime_error("generator matrix size does not align with the size of weights");  
    }

    int message_size = min((int)weights.size(), h);
    int total_weight = 0;
    for(int i = 0; i < v; ++i){
        int wt = 0;
        for(int j = 0; j < message_size; ++j){
            wt += G.element(i, j) * weights[j];
        }
        total_weight += wt;
    }
    return total_weight;
}

vector<int> rowWeights(const Matrix &stabilizers){
    vector<int> weights;
    int vs = stabilizers.size_v();
    for(int i = 0; i < vs; ++i){
        int wt = weight(stabilizers.getRow(i));
        weights.push_back(wt);
    }
    return weights;
}

Chromosome::Chromosome(const Chromosome & other): ChromosomeBase(other), h(other.h){
    mat = new int*[other.v];
    v = other.v;
    for(int i = 0; i < v; ++i) mat[i] = genes + i*h;
    matrix = Matrix(v, h);
    params = other.params;
}

Chromosome::Chromosome(vector<int> & stabilizer_weights, Parameters & params, int v, int h, Initializer init):
    ChromosomeBase(v*h, init), 
    mat(NULL), v(v), h(h), matrix(v, h), params(&params), stabilizer_weights(&stabilizer_weights)
{
    mat = new int*[v];
    for(int i = 0; i < v; ++i){
        mat[i] = genes + h*i;
    } 
}

Chromosome::Chromosome(Matrix default_mat, vector<int> & stabilizer_weights, Parameters & params, int v, int h, Initializer init):
    ChromosomeBase(v*h, init), 
    mat(NULL), v(v), h(h), matrix(default_mat), params(&params), stabilizer_weights(&stabilizer_weights)
{
    mat = new int*[v];
    for(int i = 0; i < v; ++i){
        mat[i] = genes + h*i;
        for(int j = 0; j < h; ++j){
            mat[i][j] = matrix.element(i, j);
        }
    } 

    cout << "default matrix: " << matrix.print() << endl;
}

Chromosome::Chromosome(string filename, vector<int> & stabilizer_weights, Parameters & params, int v, int h, Initializer init):
    ChromosomeBase(v*h, init), 
    mat(NULL), v(v), h(h), matrix(v, h), params(&params), stabilizer_weights(&stabilizer_weights)
{
    ifstream file(filename);
    vector<vector<int>> matrix;
    string line;

    while (getline(file, line)) {
        vector<int> row;
        stringstream ss(line);
        string value;

        while (getline(ss, value, ',')) {
            row.push_back(stoi(value));
        }
        matrix.push_back(row);
    }

    file.close();


    mat = new int*[v];
    for(int i = 0; i < v; ++i){
        mat[i] = genes + h*i;
        for(int j = 0; j < h; ++j){
            mat[i][j] = matrix.at(i).at(j);
        }
    } 

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
    optimized_matrix = matrixOptimization(matrix);
    generator_matrix = nullSpace(optimized_matrix);
}

double Chromosome::computeFitnessValue(){
    // convert genes to a matrix
    updateMatrix(); 

    fitness_val = 0;

    int gr = girth(optimized_matrix);
    if(gr < 0){
        gr = 2;
    }
    fitness_val += params->alpha*encodingDepth(generator_matrix, *stabilizer_weights);
    fitness_val += params->beta*countingDepth(optimized_matrix);
    fitness_val += params->gamma*correctionDepth(optimized_matrix);
    int min_distance = minimumDistance(matrix);
    fitness_val += params->delta *(1/floor((min_distance-1)/2) + 1/log10(pow(min_distance, 8) + 1e-10));
    fitness_val += params->zeta / pow(log10(gr), 8);
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
    params = other.params;

    return *this;
}


GeneticAlgorithm::GeneticAlgorithm(){
    pop_size = 0;
    chromosome_size = 0;
    crossover_rate = mutation_rate = 0;
    elite_rate = roulete_rate = 0;
}

GeneticAlgorithm::GeneticAlgorithm(
    Matrix stabilizers,
    int initial_population_size, 
    int v, int h,
    double c_rate,
    double m_rate, 
    double e_rate, 
    double r_rate,
    Parameters params,
    vector<Matrix> injection_files
    ):
    v(v), h(h),
    crossover_rate(c_rate),
    mutation_rate(m_rate),
    elite_rate(e_rate),
    roulete_rate(r_rate),
    chromosome_size(v*h),
    pop_size(initial_population_size),
    params(params)
{

    stabilizer_weights = rowWeights(stabilizers); 

    int number_of_files = injection_files.size();
    int number_of_random_population = initial_population_size - number_of_files;

    for(int i = 0; i < number_of_files; ++i){
        population.push_back(new Chromosome(injection_files[i], stabilizer_weights, params, v, h));
    }

    for(int i = 0 ; i < number_of_random_population; ++i){
        population.push_back(new Chromosome(stabilizer_weights, params, v, h));
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
        ready_to_use.push_back(new Chromosome(this->stabilizer_weights, params, v, h, RANDOM)); 
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
        ChromosomeBase::mutation(*population[p], *o);
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

void GeneticAlgorithm::populationInjection(vector<Chromosome *> pop){
    int size = min(pop.size(), population.size());
    for(int i = 0; i < size; ++i){
        *population[i] = *pop[i]; 
    }
}


Chromosome GeneticAlgorithm::run(int iterations, vector<map<string, string> > & iteration_data){

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

        cout << "iteration " << i << endl;
        map<string, string> data = {
            {"iteration", to_string(i)},
            {"fitness value", to_string(elites[0]->fitnessValue())},
            {"counting depth", to_string(countingDepth(elites[0]->optimized_matrix))},
            {"correction depth", to_string(correctionDepth(elites[0]->optimized_matrix))},
            {"encoding cost", to_string(encodingDepth(elites[0]->generator_matrix, stabilizer_weights))},
            {"minimum distance", to_string(minimumDistance(elites[0]->optimized_matrix))},
            {"girth", to_string(girth(elites[0]->optimized_matrix))}
        };
        cout << "fitness value " << elites[0]->fitnessValue() << " girth : " << girth(elites[0]->optimized_matrix) << endl;
        iteration_data.push_back(data);
    }
    return *population[0];
}

HeuristicPEG::HeuristicPEG(int n, int m, vector<int> degree): PEG(n, m, degree), _matrix(n, m){}
HeuristicPEG::HeuristicPEG(int v, int h, int degree): PEG(v, h, degree), _matrix(v, h){}

void HeuristicPEG::connect(Node * symbol, Node * check){
    PEG::connect(symbol, check);
    _matrix.setElement(check->index(), symbol->index(), 1);
}


Node * HeuristicPEG::pickup(Node * node, vector<Node *> nodes){
    // choose one node with the lowest cost
    vector<pair<Node *, int> > costs;
    Matrix tmp;
    tmp = _matrix;
    for(int i = 0; i < nodes.size(); ++i){
        tmp.setElement(nodes[i]->index(), node->index(), 1);
        int cost = countingDepth(tmp) + correctionDepth(tmp);
        costs.push_back({nodes[i], cost});

        tmp.setElement(nodes[i]->index(), node->index(), 0); // flip back
    }

    sort(costs.begin(), costs.end(), [](pair<Node *, int> p1, pair<Node *, int> p2){
        return p1.second < p2.second;
    });

    return costs[0].first;
}