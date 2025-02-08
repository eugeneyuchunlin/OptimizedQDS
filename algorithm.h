#ifndef __ALGORITHM_H__
#define __ALGORITHM_H__

#include <iostream>
#include <vector>

#include "matrix.h"

double sparsity(const Matrix &mat);
vector<vector<int> > codewords(const Matrix &mat);
int minimumDistance(const vector<vector<int> > &codewords);

#endif