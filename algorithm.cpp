#include "algorithm.h"
#include "matrix.h"

#include <iostream>
#include <vector>

using namespace std;

double sparsity(const Matrix &mat){
    int count = 0;
    for(int i = 0; i < mat.size_v(); ++i){
        for(int j = 0; j < mat.size_h(); ++j){
            if(mat.element(i, j) != 0){
                count++;
            }
        }
    }
    return (double)count / (mat.size_v() * mat.size_h());
}

vector<vector<int> > codewords(const Matrix &mat){
    vector<vector<int> > codewords;
    int v = mat.size_v();
    int h = mat.size_h();
    int num_codewords = 1 << v;
    for(int i = 0; i < num_codewords; ++i){
        vector<int> codeword;
        codeword.assign(mat.size_h(), 0);
        for(int j = 0; j < v; ++j){
            if((i >> j) & 1){
                for(int k = 0; k < h; ++k){
                    codeword[k] += mat.mat[j][k];
                }
            }
        } 

        for(int k = 0; k < h; ++k){
            codeword[k] %= 2;
        }
        codewords.push_back(codeword);
    }

    return codewords;
}

int weight(const vector<int> codeword){
    int result = 0;
    for(unsigned int i = 0; i < codeword.size(); ++i){
        result += codeword[i];
    }
    return result;
}

int minimumDistance(const vector<vector<int> > &codewords){
    int mini_d = codewords[0].size();
    for(unsigned int i = 0; i < codewords.size(); ++i){
        int w = weight(codewords[i]);
        if(w < mini_d && w != 0){
            mini_d = w;
        }
    }
    return mini_d;
}

inline int columnDegree(const vector<vector<int> > &mat, int col){
    int deg = 0;
    for(unsigned int i = 0; i < mat.size(); ++i){
        deg += mat[i][col];
    }
    return deg;
}

vector<int> degree(const Matrix &mat){
    int v = mat.size_v();
    int h = mat.size_h();

    vector<int> degrees;
    degrees.assign(h, 0);

    for(int i = 0; i < h; ++i){
        degrees[i] = columnDegree(mat.mat, i);
    }

    return degrees;
}


// this implies that we want to have 
int countingDepth(const Matrix &mat){
    
    int h = mat.size_h();
    int v = mat.size_v();

    int depth = 0;

    for(int i = 0; i < h; ++i){
        int deg = columnDegree(mat.mat, i);
        
        // cout << "deg " << i << " " << deg << " ";

        int single_depth = 0;
        for(int j = 1; j <= deg; ++j){
            // floor(log2(j)) + 1
            single_depth += ((sizeof(j)<<3) - __builtin_clz(j));
        }
        // cout << "single depth " << single_depth << endl;
        depth += single_depth;
    }

    return depth;
}

int correctionDepth(const Matrix &mat){

    int v = mat.size_v();
    int h = mat.size_h();

    vector<int> degrees = degree(mat);
    vector<int> stat_degree;
    stat_degree.assign(v+1, 0);
    for(int i = 0; i < h; ++i){
        ++stat_degree[degrees[i]];
    }

    // for(int i = 0; i < v+1; ++i){
    //     cout << i << ": " << stat_degree[i] << endl;
    // }

    int depth = 0;
    for(int i = 1; i < v+1; ++i){
        for(int j = i; j < v+1; ++j){
            depth += stat_degree[j];
        }
    }

    return depth;
}