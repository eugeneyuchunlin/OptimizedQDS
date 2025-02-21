#include <exception>
#include <iostream>
#include <set>
#include <unordered_set>
#include <queue>
#include <utility>

#include "matrix.h"
#include "common.h"

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

Matrix Matrix::transpose() const{
    Matrix m(h, v);
    for(int i = 0; i < m.v; ++i){
        for(int j = 0; j < m.h; ++j){
            m.mat[i][j] = mat[j][i];
        }
    }
    return m;
}

string Matrix::print() const {
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

vector<int> Matrix::getRow(int i) const {
    return mat[i];
}

vector<int> Matrix::getColumn(int j)const {
    vector<int> col;
    for(int i = 0; i < v; ++i) col.push_back(mat[i][j]);
    return col;
}


Matrix rref(const Matrix& mat){
    // rearrange the vectors
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


    // backward elimination
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


Matrix nullSpace(const Matrix& mat) {
    if (mat.v == 0 || mat.h == 0) {
        return Matrix(); // Handle empty matrix
    }

    int num_rows = mat.v;
    int num_cols = mat.h;

    std::vector<int> pivot_cols;
    std::vector<std::vector<int>> rref = mat.mat;

    // Perform Gaussian elimination in GF(2)
    int row = 0;
    for (int col = 0; col < num_cols; ++col) {
        int pivot_row = -1;
        for (int i = row; i < num_rows; ++i) {
            if (rref[i][col] == 1) {
                pivot_row = i;
                break;
            }
        }

        if (pivot_row != -1) {
            pivot_cols.push_back(col);
            if (pivot_row != row) {
                std::swap(rref[row], rref[pivot_row]);
            }

            for (int i = 0; i < num_rows; ++i) {
                if (i != row && rref[i][col] == 1) {
                    for (int j = col; j < num_cols; ++j) {
                        rref[i][j] ^= rref[row][j]; // XOR (addition in GF(2))
                    }
                }
            }
            row++;
        }
    }

    // Identify free variables
    std::vector<int> free_cols;
    for (int col = 0, pivot_index = 0; col < num_cols; ++col) {
        if (pivot_index < pivot_cols.size() && col == pivot_cols[pivot_index]) {
            pivot_index++;
        } else {
            free_cols.push_back(col);
        }
    }

    // Construct null space vectors
    std::vector<std::vector<int>> nullspace;
    for (int free_col : free_cols) {
        std::vector<int> null_vector(num_cols, 0);
        null_vector[free_col] = 1;

        for (int i = 0, pivot_index = 0; i < num_cols; ++i) {
            if (pivot_index < pivot_cols.size() && i == pivot_cols[pivot_index]) {
                null_vector[i] = rref[pivot_index][free_col];
                pivot_index++;
            }
        }
        nullspace.push_back(null_vector);
    }

    return Matrix(nullspace).transpose();
}


// incomplete function
Matrix prioritizedGaussianElimination(const Matrix &mat){
    Matrix m = rref(mat);    
    int v = m.v;
    int h = m.h;
    
    int rank = 0;
    // calculate the rank
    for(int i = 0; i < v; ++i){
        if(mat.mat[i][i] !=0) ++rank;
    } 

    return m;
}

Matrix matrixOptimization(const Matrix & m) {
    vector<int> removeRows, removeCols; // Store indices to be removed

    // Identify rows with weight 1 and their corresponding column indices
    for(int i = 0; i < m.v; ++i){
        int one_index = -1;                  
        int count = 0; 

        for(int j = 0; j < m.h; ++j){
            if(m.mat[i][j] == 1) {
                count++;
                one_index = j;
            }
        }

        if(count == 1) {
            removeRows.push_back(i);
            removeCols.push_back(one_index);
        }
    }

    // Remove columns by marking them
    unordered_set<int> removeColSet(removeCols.begin(), removeCols.end());

    // Create a new matrix without the selected rows and columns
    Matrix newMatrix;
    newMatrix.v = m.v - removeRows.size();
    newMatrix.h = m.h - removeColSet.size();
    newMatrix.mat.resize(newMatrix.v, vector<int>(newMatrix.h, 0));

    int newRow = 0;
    for(int i = 0; i < m.v; ++i) {
        if(find(removeRows.begin(), removeRows.end(), i) != removeRows.end()) 
            continue; // Skip rows marked for removal

        int newCol = 0;
        for(int j = 0; j < m.h; ++j) {
            if(removeColSet.count(j)) 
                continue; // Skip columns marked for removal

            newMatrix.mat[newRow][newCol] = m.mat[i][j];
            newCol++;
        }
        newRow++;
    }

    return newMatrix;
}

int _cycle(const std::vector<std::vector<int>>& adj_matrix, int start) {
    int n = adj_matrix.size();
    if (n == 0 || start < 0 || start >= n) return -1;

    std::queue<std::pair<int, std::unordered_set<int>>> q; // Queue: (node, path)
    q.push({start, {start}}); // Start with initial node and its path

    while (!q.empty()) {
        int current_node = q.front().first;
        std::unordered_set<int> current_path = q.front().second;
        q.pop();

        for (int neighbor = 0; neighbor < n; ++neighbor) {
            if (adj_matrix[current_node][neighbor]) {
                if (neighbor == start && current_path.size() > 2) {
                    return current_path.size(); // Cycle found
                }
                if (current_path.find(neighbor) == current_path.end()) {
                    std::unordered_set<int> new_path = current_path;
                    new_path.insert(neighbor);
                    q.push({neighbor, new_path});
                }
            }
        }
    }
    return -1; // No cycle found
}

int girth(const Matrix & mat){

    // We represent the graph in an adjancy matrix. The size of the matrix is (# variable nodes, #check nodes)
    // # number of variable nodes = mat.h
    // # check nodes = mat.v 
    Matrix tanner_graph(mat.v + mat.h, mat.v + mat.h);
    Matrix transposed = mat.transpose();

    for(int i = 0; i < mat.v; ++i)
        for(int j = 0; j < mat.h; ++j)
            tanner_graph.mat[i][j + mat.v] = mat.mat[i][j];

    vector<int> nodes; 
    for(int i = 0; i < mat.h; ++i){
        int row_weight = 0;
        for(int j = 0; j < mat.v; ++j){
            tanner_graph.mat[i+mat.v][j] = transposed.mat[i][j];
            row_weight += transposed.mat[i][j];
        }
        if(row_weight > 1) nodes.push_back(i+mat.v);
    }
    // cout << tanner_graph.print() << endl;
    // cout << mat.print() << endl;
    int _girth = numeric_limits<int>::max();
    bool found_flag = false;
    for(unsigned int i = 0; i < nodes.size(); ++i){
        int cycle = _cycle(tanner_graph.mat, nodes[i]);
        if(cycle > 0 && cycle < _girth){
            _girth = cycle;
            found_flag = true;
        }
    ///     cout << "cycle of " << nodes[i] - mat.v << ": " << cycle << endl;

        if(_girth == 4) return 4;
    }

    return found_flag ? _girth : -1;
}