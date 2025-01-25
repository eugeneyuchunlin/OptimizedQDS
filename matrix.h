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
    Matrix();
    Matrix(int v, int h);
    Matrix(const vector<vector<int> > &);
    int size_v();
    int size_h();
    int element(int i, int j);

    void setElement(int i, int j, int val);
    
    Matrix transpose();
    Matrix operator*(const Matrix& mat);

    string print();

    friend Matrix rref(const Matrix& mat);
    friend Matrix nullSpace(const Matrix& mat);
};

Matrix rref(const Matrix& mat);
Matrix nullSpace(const Matrix &mat);

#endif