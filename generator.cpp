#include "generator.h"
#include "expression.h"

#include <iostream>
#include <string>
using namespace std;


Generator::Generator(int base, int cardinality):
            exponent(0),
            coefficient(0),
            zero(true),
            cardinality(cardinality),
            base(base)
{

}

Generator::Generator(int exponent, int coefficient, int base, int cardinality):
            exponent(exponent),
            coefficient(coefficient),
            zero(false),
            cardinality(cardinality),
            base(base)
{
    coefficient %= base;

    if(coefficient == 0){
        zero = true;
    }
}

Generator Generator::addInSameExponent(const Generator & g) const {
    // TODO: assert if they are not in the same cadinality and exponent
    return Generator(exponent, (coefficient + g.coefficient)%base, base, cardinality);
}

Expression Generator::operator+(const Generator& g1){
    Expression exp(base, cardinality);

    if(!(zero || g1.zero)){
        // check if their exponents are the same
        if (exponent == g1.exponent){
            int new_coefficent = (coefficient + g1.coefficient) % base;
            Generator result(exponent, new_coefficent, base, cardinality);
            exp.push(result);           
        }else{
            exp.push(g1);
            exp.push(*this);
        }
    }else if (zero){
        exp = exp + g1;
    }else if (g1.zero){
        exp.push(*this);
    }        
    return exp;
}

Generator& Generator::operator=(const Generator& other) {
    if (this != &other){
        exponent = other.exponent;
        coefficient = other.coefficient;
        zero = other.zero;
    }

    return *this;
}

string Generator::print(bool sign){
    string s; 
    if(zero){
        return to_string(0);
    }

    if(sign){
        s.push_back('+');
    }
    s += to_string(coefficient);
    if(exponent){
        s += "w^";
        s += to_string(exponent);
    }
    return s;
}

Generator::Generator(const Generator &g): 
        exponent(g.exponent), 
        coefficient(g.coefficient), 
        zero(g.zero), 
        cardinality(g.cardinality), 
        base(g.base)
{
}