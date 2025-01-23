#ifndef __EXPRESSION_H__
#define __EXPRESSION_H__

#include <vector>
#include <string>

using namespace std;

class Generator;

class Expression{
    protected:
        vector<Generator> poly;
        const int cardinality;
        const int base;
    public:

        Expression(int base, int cardinality);
        void push(const Generator & g1);
        string print();
        Expression operator+(const Generator & g);
        Expression &operator=(const Expression & exp);
        Expression(Expression &);

};

#endif