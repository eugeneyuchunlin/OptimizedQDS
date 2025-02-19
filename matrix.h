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
    inline int size_v() const{
        return v;
    }
    inline int size_h() const{
        return h;
    }
    inline int element(int i, int j) const{
        return mat[i][j];
    }
    inline void setElement(int i, int j, int val){
        if(i >= v || j >= h){
            throw std::out_of_range("Index out of range");
        }else{
            mat[i][j] = val;
        }
    }
    
    Matrix transpose() const;
    Matrix operator*(const Matrix& mat);

    string print() const;

    friend Matrix rref(const Matrix& mat);
    friend Matrix nullSpace(const Matrix& mat);
    friend vector<vector<int> > codewords(const Matrix &mat);
    friend Matrix prioritizedGaussianElimination(const Matrix &mat);
    friend int countingDepth(const Matrix &mat);
    friend int correctionDepth(const Matrix &mat);
    friend vector<int> degree(const Matrix &mat);
    friend Matrix matrixOptimization(const Matrix &mat);
    friend int girth(const Matrix &mat);
};

Matrix rref(const Matrix& mat);
Matrix nullSpace(const Matrix &mat);
Matrix prioritizedGaussianElimination(const Matrix &mat);
Matrix matrixOptimization(const Matrix &mat);
int girth(const Matrix &mat);

#endif