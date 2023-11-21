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
// the implementation of theo_model_v2.ipynb

int main() {
    const string filePath = "../data.txt"; // Change this to the appropriate file path

    vector<vector<float>> model = datareader(filePath);

    // Define parameters
    int nz = 300, nx = 298;
    float dz = 5.0, dx = 5.0;
    float sx = 1000, sy = 800;
    int nt = 20;
    float CC = 0.5;

    // Create velocity function
    vector<float> max_v = *max_element(model.begin(), model.end());
    double vmax = 5560.32800;
    double dt = CC * dx / vmax;
    // Initialize wavefields
    vector<vector<float>> UUo(nx, vector<float>(nz, 0.0));
    vector<vector<float>> UUm(nx, vector<float>(nz, 0.0));

    // Define time and forcing function
    vector<float> t(nt);
    vector<double> F(nt);
    float ss = 1;

    for (int i = 0; i < nt; ++i) {
        t[i] = i * dt;
        F[i] = (1 - ((t[i] - 40 * t[i]) / ss) * ((t[i] - 40 * t[i]) / ss)) * exp(-(t[i] - 40 * t[i]) * (t[i] - 40 * t[i]) / (2 * ss * ss));
    }

    // Total Solution space
    vector<vector<vector<float>>> fff(nx, vector<vector<float>>(nz, vector<float>(nt, 0.0)));

    // Iterate over solution
    clockid_t start, end;
    start = clock();
    // calculating the wave
    for (int it = 0; it < nt; ++it) {
        vector<vector<float>> tmp = awe_2d_explicit_solver_heterogeneous_8th_order(UUo, UUm, dx, dz, dt, model, F, it, sx, sy);
        for (int j = 0; j < nx; ++j) {
            for (int k = 0; k < nz; ++k) {
                fff[j][k][it] = tmp[j][k];
                //UUm[j][k] = UUo[j][k];
                //UUo[j][k] = tmp[j][k];
            }
        }
        UUm = UUo;
        UUo = tmp;
    }
    end = clock();
    cout<<"Time Used "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;

    // IO
    for (int t = 0; t < nt; t++){
        string filename = "time_step" + to_string(t) + ".txt";
        ofstream outfile;
        outfile.open(filename, ios::out | ios::trunc);

        for(int i=0; i<nx; i++){
            for (int j=0; j<nz; j++){
                outfile << fff[i][j][t] << " ";
            }
            outfile << endl;
        }
        outfile.close();
    }

    return 0;

}