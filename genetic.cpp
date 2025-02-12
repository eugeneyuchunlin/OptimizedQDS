#include "genetic.h"


#include <iostream>
using namespace std;

ChromosomeBase::ChromosomeBase(int size, Initializer init):size(size){
    // genes.assign(size, 0);
    genes = new int[size];
    if(genes == NULL){
        cout << "Memory allocation failed. New memory for int" << endl; 
    }else{
        memset(genes, 0, sizeof(int)*size);
        if(init == RANDOM) randomizeGenes();
    }
}

void ChromosomeBase::randomizeGenes(){
    for(int i = 0; i < size; ++i){
        genes[i] = rand() & 1; // only set it to 0 or 1
    }
}

string ChromosomeBase::print(){
    string s;
    for(int i = 0; i < size; ++i){
        s += to_string(genes[i]) + " ";
    }
    return s;
}

ChromosomeBase::~ChromosomeBase(){
    if(genes){
        delete[] genes;
        genes = NULL;
    }
}

void ChromosomeBase::copyFrom(const ChromosomeBase & other){
    this->genes = new int[other.size];
    this->size = other.size;
    memcpy(this->genes, other.genes, sizeof(int)*this->size);
}

ChromosomeBase::ChromosomeBase(const ChromosomeBase &other){
    copyFrom(other); 
}

ChromosomeBase & ChromosomeBase::operator=(const ChromosomeBase & other){
    if(this->size != other.size){
        if(this->genes != nullptr){
            delete[] this->genes;
        }
        this->genes = new int[other.size];
        this->size = other.size;
    }
    memcpy(this->genes, other.genes, sizeof(int)*this->size);
    return *this;
}

void ChromosomeBase::crossover(const ChromosomeBase & c1, const ChromosomeBase & c2, ChromosomeBase & o1, ChromosomeBase & o2){
    if(c1.size != c2.size){
        cout << "ChromosomeBase size doesn't align" << endl; 
    }

    int size = c1.size;

    int cut1, swp_size;
    cut1 = rand() % size;
    swp_size = rand() % (size - cut1);


    // cout << "cut : " << cut1 << endl;
    // cout << "swp_size: " << swp_size << endl;

    o1 = c1;
    o2 = c2;

    memcpy(o1.genes + cut1, c2.genes + cut1, sizeof(int)*(swp_size));
    memcpy(o2.genes + cut1, c1.genes + cut1, sizeof(int)*(swp_size));
}

void ChromosomeBase::mutation(const ChromosomeBase &c, ChromosomeBase &o){
    o = c;
    int point = rand() % (c.size);
    o.genes[point] ^= 1;
}