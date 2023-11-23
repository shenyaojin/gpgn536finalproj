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
/*
    PARAMETER TABLE:
    Set up eight-order solver of the acoustic wave equation
    usage: U=awe_2d_explicit_solver_heterogeneous_8th_order(UUo,UUm,dx,dy,dt,v,F,it,sx,sy):
    input:
        UUo: Acoustic pressure vector (nx,ny) at time step n
        UUm: Acoustic pressure vector (nx,ny) at time step n-1
        dx : Spatial sampling in x
        dy : Spatial sampling in y
        dt : Temporal sampling
        v  : Heterogeneous propagation velocity (nx,ny)
        F  : Forcing function (nt)
        it : Time index
        sx : Location of source in x (meters). In this case, sx should be an array.
        sy : Location of source in y (meters). In this case, sx should be an array.
    output:
        UUm: Acoustic pressure vector (nx,ny) at time n+1
    dependencies:
        None
    written by Jeff Shragge, jshragge@mines.edu, 10/2019
    extended to multiple source case by Shenyao Jin, shenyaojin@mines.edu, 11/2023
*/
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
                                     dtdx2 * v[i][j] * v[i][j] * (
                                             -1 / 560.0  * UUo[i - 4][j]
                                             + 8 / 315.0 * UUo[i - 3][j]
                                             - 1 / 5.0   * UUo[i - 2][j]
                                             + 8 / 5.0   * UUo[i - 1][j]
                                             - 205 / 72.0* UUo[i][j]
                                             + 8 / 5.0   * UUo[i + 1][j]
                                             - 1 / 5.0   * UUo[i + 2][j]
                                             + 8 / 315.0 * UUo[i + 3][j]
                                             - 1 / 560.0 * UUo[i + 4][j]
                                     ) + dtdy2 * v[i][j] * v[i][j] * (
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

vector<vector<float>> awe_2d_heterogeneous_8th_order_data_time_reverse(
        vector<vector<float>>& UUo, vector<vector<float>>& UUm, float dx, float dy, double dt, vector<vector<float>>& v, vector<vector<float>>& D, int it, float ry){
// 2D acoustic wave equation reverse solver for heterogeneous media. The accuracy is 8th order.
// About the variables, please refer to the python function with the same name.
/*
    PARAMETER TABLE:
    Set up eight-order solver of the acoustic wave equation
    usage: U=awe_2d_heterogeneous_8th_order_data_time_reverse(UUo,UUm,dx,dy,dt,v,D,it,ry):
    input:
        Uo: Acoustic pressure vector (nx,ny) at time step n
        Um: Acoustic pressure vector (nx,ny) at time step n-1
        dx: Spatial sampling in x
        dy: Spatial sampling in y
        dt: Temporal sampling
        v : Heterogeneous propagation velocity (nx,ny)
        D : Data to be time reversed
        it: Time index
        ry: Injection depth (meters)
    output:
        Um: Acoustic pressure vector (nx,ny) at time n+1
    dependencies:
        None
*/
    int nx = UUo.size();
    int ny = UUo[0].size();
    // define courant numbers (squared)
    float dtdx2 = pow(dt / dx, 2);
    float dtdy2 = pow(dt / dy, 2);
    // source location
    int iry = int(ry / dy);

    // inject wavelet
    for(int iter = 0; iter < nx; iter ++){
        UUo[iter][iry] += dt * dt * D[iter][it];
    }

    for (int i = 4; i < nx - 4; ++i) {
        for (int j = 4; j < ny - 4; ++j) {
            UUm[i][j] = 2 * UUo[i][j] - UUm[i][j] +
                        dtdx2 * v[i][j] * v[i][j] * (
                                -1 / 560.0  * UUo[i - 4][j]
                                + 8 / 315.0 * UUo[i - 3][j]
                                - 1 / 5.0   * UUo[i - 2][j]
                                + 8 / 5.0   * UUo[i - 1][j]
                                - 205 / 72.0* UUo[i][j]
                                + 8 / 5.0   * UUo[i + 1][j]
                                - 1 / 5.0   * UUo[i + 2][j]
                                + 8 / 315.0 * UUo[i + 3][j]
                                - 1 / 560.0 * UUo[i + 4][j]
                        ) + dtdy2 * v[i][j] * v[i][j] * (
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
        } // end j
    } // end i
    return UUm;
}