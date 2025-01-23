#include <exception>
#include "matrix.h"

Matrix::Matrix(int v, int h){
    this->v = v;
    this->h = h;
    for(int i = 0; i < v; ++i){
        vector<int> v;
        for(int j = 0; j < h; ++j){
            v.push_back(0); 
        }
        mat.push_back(v);
    }
}

int Matrix::element(int i, int j){
    return mat[i][j];
}

int Matrix::size_v(){
    return v;
}

int Matrix::size_h(){
    return h;
}

Matrix Matrix::transpose(){
    Matrix m(h, v);
    for(int i = 0; i < m.v; ++i){
        for(int j = 0; j < m.h; ++j){
            m.mat[i][j] = mat[j][i];
        }
    }
    return m;
}

void Matrix::setElement(int i, int j, int val){
    if(i >= v || j >= h){
        throw std::out_of_range("Index out of range");
    }else{
        mat[i][j] = val;
    }
}

string Matrix::print(){
    string s;
    for(int i = 0; i < v; ++i){
        for(int j = 0; j < h; ++j){
            s += to_string(mat[i][j]);
            s += " ";
        }
        s += "\n";
    }
    return s;
}

Matrix Matrix::operator*(const Matrix& m){
    if(h != m.v){
        throw runtime_error("matrix size not correct");
    }else{
        Matrix result(v, m.h);
        for(int i = 0; i < result.v; ++i){
            for(int j = 0; j < result.h; ++j){

                for(int k = 0; k < h; ++k){
                    result.mat[i][j] += mat[i][k] * m.mat[k][j];                    
                }
            }
        }
        return result;
    }
}