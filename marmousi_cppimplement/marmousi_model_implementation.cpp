#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <time.h>
#include "lib/pdesolver.h"

using namespace std;
// written by Shenyao Jin, shenyaojin@mines.edu, 11/2023
// part of the final project of GPGN 536
// the implementation of marmousi_model.ipynb

int main() {
    const string filePath = "../data.txt"; // Change this to the appropriate file path

    vector<vector<float>> v = datareader(filePath);

    // Define parameters
    int nz = 300, nx = 298;
    float dz = 5.0, dx = 5.0;
    float oz = 0.0, ox = 0.0;
    vector<float> z(nz);
    vector<float> x(nx);

    // Linearly spaced vector for z
    for(int i = 0; i < nz; ++i) {
        z[i] = oz + i * dz;
    }

    // Linearly spaced vector for x
    for(int i = 0; i < nx; ++i) {
        x[i] = ox + i * dx;
    }
    vector<vector<float>> UUo(nx, vector<float>(nz, 0));
    vector<vector<float>> UUm(nx, vector<float>(nz, 0));

    float CC = 0.5;
    int nt  = 2500;
    float v_max = 5560.328;
    double dt = CC * dx / v_max;

    double t0 = 0.05;
    float ss = 0.01; // sigma for Ricker wavelet.
    vector<float> t(nt);
    vector<double> F(nt);
    //define F(force term) and t(time series).
    for (int i = 0; i < nt; ++i) {
        t[i] = i * dt;
        F[i] = (1 - ((t[i] - t0) / ss) * ((t[i] - t0) / ss))
               * exp(-(t[i] - t0) * (t[i] - t0) / (2 * ss * ss));
    }

    // define shot parameters
    int shot_num = 32;
    vector<float> sx(shot_num); // offset of the shot location (m)
    // source location
    // define nx
    float start = 0;
    float end = nx * dx;
    float step = (end - start) / (shot_num + 1);

    for (int i = 0; i < shot_num; ++i) {
        sx[i] = start + step * (i + 1); // Offset by one to mimic Python's [1:-1]
    }
    float sz = 25.0;  // depth of the shot location (m)
    // define forward modeling array fff
    vector<vector<vector<float>>> fff(nx, vector<vector<float>>(nt, vector<float>(nt, 0.0)));
    // image: img;
    vector<vector<float>> img(nx, vector<float>(nz, 0));
    for (int iter = 0; iter < shot_num; iter++) {
        for (int it = 0; it < nt; it++) {
            vector<vector<float>> tmp;
            tmp = awe_2d_explicit_solver_heterogeneous_8th_order(UUo, UUm, dx, dz, dt, v, F, it, sx[iter], sz);
            // pass the value:
            for (int ix = 0; ix < nx; ix++) {
                for (int iy = 0; iy < nz; iy++) {
                    fff[ix][iy][it] = tmp[ix][iy];
                } //end iy
            } // end ix
            UUm = UUo;
            UUo = tmp;
        } // end it

        // int wavefield value (reverse)
        UUo.clear(); // Clear the entire 2D vector
        UUo.resize(nx, vector<float>(nz, 0.0f)); // Resize and initialize with zeros
        UUm.clear();
        UUm.resize(nx, vector<float>(nz, 0.0f));

        vector<vector<float>> data(nx, vector<float>(nt, 0.0)); // This will store the slice
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < nt; ++k) {
                data[i][k] = fff[i][5][k]; // Copying the relevant elements
            }
        }
        // define the observation location.
        float rz = 25.0;
        vector<vector<vector<float>>> mmm(nx, vector<vector<float>>(nt, vector<float>(nt, 0.0)));



    } // end shot number
    return 0;

}