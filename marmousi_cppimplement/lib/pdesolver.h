//
// Created by shenyaojin on 11/21/23.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <time.h>

using namespace std;

#ifndef CPPIMPLEMENT_PDESOLVER_H
#define CPPIMPLEMENT_PDESOLVER_H

vector<vector<float>> awe_2d_explicit_solver_heterogeneous_8th_order(
        vector<vector<float>>& UUo, vector<vector<float>>& UUm,
        float dx, float dy, double dt, vector<vector<float>>& v,
        vector<double>& F, int it, float sx, float sy);

vector<vector<float>> awe_2d_heterogeneous_8th_order_data_time_reverse(
        vector<vector<float>>& UUo, vector<vector<float>>& UUm, float dx,
        float dy, double dt, vector<vector<float>>& v,
        vector<vector<float>>& D, int it, float ry);

vector<vector<float>> datareader(const string& filePath);

#endif //CPPIMPLEMENT_PDESOLVER_H