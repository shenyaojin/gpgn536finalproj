#include <iostream>
#include <vector>
#include <cmath>

//
// Created by shenyaojin on 11/21/23. shenyaojin@mines.edu
//

using namespace std;

vector<vector<float>> awe_2d_explicit_solver_heterogeneous_8th_order(
        vector<vector<float>>& UUo, vector<vector<float>>& UUm,
        float dx, float dy, double dt, vector<vector<float>>& v,
        vector<double>& F, int it, float sx, float sy) {
// 2D acoustic wave equation solver for heterogeneous media. The accuracy is 8th order.
// About the variables, please refer to the python function with the identical name.
    int nx = UUo.size();
    int ny = UUo[0].size();

    float dtdx2 = pow(dt / dx, 2);
    float dtdy2 = pow(dt / dy, 2);

    int isx = static_cast<int>(sx / dx);
    int isy = static_cast<int>(sy / dy);

    UUo[isx][isy] += dt * dt * F[it];

    vector<vector<float>> updatedWavefield(nx, vector<float>(ny, 0.0));

    for (int i = 4; i < nx - 4; ++i) {
        for (int j = 4; j < ny - 4; ++j) {
            updatedWavefield[i][j] = 2 * UUo[i][j] - UUm[i][j] +
                                     dtdx2 * pow(v[i][j], 2) * (
                                             -1 / 560.0  * UUo[i - 4][j]
                                             + 8 / 315.0 * UUo[i - 3][j]
                                             - 1 / 5.0   * UUo[i - 2][j]
                                             + 8 / 5.0   * UUo[i - 1][j]
                                             - 205 / 72.0* UUo[i][j]
                                             + 8 / 5.0   * UUo[i + 1][j]
                                             - 1 / 5.0   * UUo[i + 2][j]
                                             + 8 / 315.0 * UUo[i + 3][j]
                                             - 1 / 560.0 * UUo[i + 4][j]
                                     ) + dtdy2 * pow(v[i][j], 2) * (
                    -1 / 560.0  * UUo[i][j - 4]
                    + 8 / 315.0 * UUo[i][j - 3]
                    - 1 / 5.0   * UUo[i][j - 2]
                    + 8 / 5.0   * UUo[i][j - 1]
                    - 205 / 72.0* UUo[i][j]
                    + 8 / 5.0   * UUo[i][j + 1]
                    - 1 / 5.0   * UUo[i][j + 2]
                    + 8 / 315.0 * UUo[i][j + 3]
                    - 1 / 560.0 * UUo[i][j + 4]
            );
        }
    }

    return updatedWavefield;
}