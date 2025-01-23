#include "expression.h"
#include "generator.h"

#include <string>

Expression::Expression(int base, int cardinality):cardinality(cardinality), base(base){
    for(int i = 0; i < cardinality-1; ++i){
        poly.push_back(Generator(0, 0, base, cardinality));
    }
}

void Expression::push(const Generator & g1){
    poly[g1.exponent] = g1.addInSameExponent(poly[g1.exponent]);
}

Expression Expression::operator+(const Generator & g){
    Expression exp = *this;
    if(g.zero){
        return exp;
    } else{
        exp.poly[g.exponent] = g.addInSameExponent(exp.poly[g.exponent]);
    }
    return exp;
}

string Expression::print(){
    string s;
    for(int i = 0; i < cardinality-1; ++i){
        if(!poly[i].zero){
            s += poly[i].print(true);
        }
    }
    if(s.length() == 0){
        return "0";
    }
    return s;
}

Expression & Expression::operator=(const Expression &exp){
    if(this != &exp){
        poly = exp.poly;
    }
    return *this;
}

Expression::Expression(Expression &exp):cardinality(exp.cardinality), base(exp.base){
    poly = exp.poly;
}