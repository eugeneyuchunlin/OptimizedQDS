#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>
#include <string>

using namespace std;


class Matrix{
protected:
    vector<vector<int> > mat; 
    int v;
    int h;
public:

    Matrix(int v, int h);
    int size_v();
    int size_h();
    int element(int i, int j);

    void setElement(int i, int j, int val);
    
    Matrix transpose();
    Matrix operator*(const Matrix& mat);

    string print();
};

#endif