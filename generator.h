#ifndef __GENERATOR_H__
#define __GENERATOR_H__

#include <string>
#include <vector>
using namespace std;

class Expression;
class Generator{
    protected:
        int exponent;
        int coefficient;
        bool zero;
        const int cardinality;
        const int base;
        vector<vector<int> > plus_relation;
        
        friend class Expression;

    public:
        Generator(int base, int cardinality);
        Generator(int exponent, int coefficient, int base, int cardinality);
        Generator(const Generator &);
        // void setPlusRelation(); // TODO
        string print(bool sign=false);
        Generator addInSameExponent(const Generator & g) const;
        Expression operator+(const Generator& g1);
        Generator& operator=(const Generator& other);
};


#endif