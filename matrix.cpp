#include <exception>
#include <iostream>
#include <set>

#include "matrix.h"

Matrix::Matrix():v(0), h(0){}

Matrix::Matrix(const vector<vector<int> > & input){
    unsigned int _v = input.size();
    unsigned int _h = 0;
    if(_v == 0){
        _h = 0;
        mat = input;
    }else{
        _h = input[0].size();
        for(unsigned int i = 0; i < input.size(); ++i){
            if(_h != input[i].size()){
                throw runtime_error("Matrix size doesn't align");
            }
        }
    }
    v = _v;
    h = _h;
    mat = input;
}

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
                    result.mat[i][j] = (result.mat[i][j] + mat[i][k] * m.mat[k][j]) % 2; // perform modulo 2
                }
            }
        }
        return result;
    }
}


Matrix rref(const Matrix& mat){
    // rearrange the vectors
    // cout << "Hi" << endl;
    vector<vector<int> > new_mat = mat.mat;
    int max_iterations = mat.v - 1;
    for(int i = 0; i < max_iterations; ++i){
        if(new_mat[i][i] == 0){
            // find the first column that new_mat[i][i] != 0
            int non_zero_row = -1;
            for(int l = i; l < max_iterations; ++l){
                if(new_mat[l][i] != 0){
                    non_zero_row = l;
                    break;
                }
            } 
            
            if(non_zero_row == -1){
                continue;
            }else{
                // swap the row;
                vector<int> temp = new_mat[i];
                new_mat[i] = new_mat[non_zero_row];
                new_mat[non_zero_row] = temp;
            } 
        }
        // from ith row, we add it to the rows below if new_mat[i][i] != 0
        for(int l = i+1; l < mat.v; ++l){
            if(new_mat[l][i] != 0){
                // perform addition modulo 2
                for(int j = 0; j < mat.h; ++j){
                    new_mat[l][j] = (new_mat[i][j] + new_mat[l][j]) % 2;
                }
            }
        }
    }

    for(int i = mat.v - 1; i > 0; --i){
        // find the non-zero col in ith row
        int nzro_col = -1;
        for(int j = 0; j < mat.h; ++j){
            if(new_mat[i][j] != 0){
                nzro_col = j;
                break;
            }
        }
        if(nzro_col == -1){
            continue;
        }
        
        for(int j = i - 1; j >= 0; --j){
            if(new_mat[j][nzro_col] != 0){
                // add to new_mat[j]
                for(int l = 0; l < mat.h; ++l){
                    new_mat[j][l] = (new_mat[j][l] + new_mat[i][l]) % 2; 
                }
            }
        }
    }
    Matrix result;
    result.h = mat.h;
    result.v = mat.v;
    result.mat = new_mat;
    return result;
}

Matrix nullSpace(const Matrix& mat){

    const int size = 100;
    // bool nonzero_equations[size]{0};
    int equations[size][size] = {{0}};
    bool not_free_variable[size] = {0};
    // for(int i = 0; i < size; ++i)free_variable[i] = true;
    

    for(int i = 0; i < mat.v; ++i){
        if(mat.mat[i][i] != 0){
            not_free_variable[i] = true;
            for(int j = i+1; j < mat.h; ++j){
                if(mat.mat[i][j] != 0){
                    equations[i][j] = 1;
                }
            }
        }else{
            // push the equation that the row has
            int k = 0;
            for(k = i+1; k < mat.h && mat.mat[i][k] == 0; ++k);
            for(int j = k+1; j < mat.h; ++j){
                if(mat.mat[i][j] != 0){
                    equations[k][j] = 1;
                }
            }
            not_free_variable[k] = true;
        }
    }
    for(int i = 0; i < mat.h; ++i){
        if(!not_free_variable[i]){
            equations[i][i] = 1;
        }
    }

    vector<vector<int> > nullspace;
    bool has_nullspace = false;
    for(int i = 0; i < mat.h; ++i){
        has_nullspace |= !not_free_variable[i];
        if(!not_free_variable[i]){
            vector<int> vec;
            for(int j = 0; j < mat.h; ++j){
                vec.push_back(equations[j][i]);
            }
            nullspace.push_back(vec);
        }
    }
    if(!has_nullspace){
        // add the trivial solution
        vector<int> vec;
        for(int i = 0; i < mat.h; ++i){
            vec.push_back(0);
        }
        nullspace.push_back(vec);
    }
    

    Matrix result(nullspace);
    return result.transpose();
}