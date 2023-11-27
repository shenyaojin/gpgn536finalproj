#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <omp.h>

#include "lib/pdesolver.h"

using namespace std;
// written by Shenyao Jin, shenyaojin@mines.edu, 11/2023
// part of the final project of GPGN 536
// the implementation of theo_model_v2.ipynb (RTM implement on a theoretical heterogeneous model).

int main() {
    // Define spatial grid
    float Lx = 2000, Ly = 1200;       // lengths (m)
    int nx = 401, ny = 241;           // Number of points in discretization
    float dx = Lx / (nx - 1);         // Discretization intervals
    float dy = Ly / (ny - 1);

    // Create velocity function
    vector<vector<float>> v(nx, vector<float>(ny, 1500)); // Initialize the 2D vector with 1500

    // Update part of the velocity function
    for (int i = 0; i < nx; ++i) {
        for (int j = 50; j < ny; ++j) {
            v[i][j] = 2000;
        }
    }
    // Initialize UUo and UUm with zeros
    vector<vector<float>> UUo(nx, vector<float>(ny, 0));
    vector<vector<float>> UUm(nx, vector<float>(ny, 0));

    // Time stepping parameters
    float CC = 0.5;   // Courant number
    int nt = 800;     // Number of time steps
    float dt = CC * dx / 2000; // Assuming v is a 1D vector of max velocity values

    // Creating a time vector `t`
    vector<float> t(nt);
    float t0 = 0.05;
    float currentTime = 0;
    vector<double> F(nt);
    float ss = 0.01;
    // define F and t array.
    for (int i = 0; i < nt; ++i) {
        t[i] = i * dt;
        F[i] = (1 - ((t[i] - t0) / ss) * ((t[i] - t0) / ss))
                * exp(-(t[i] - t0) * (t[i] - t0) / (2 * ss * ss));
    }
    int shot_num = 10;
    float sy = 25;
    vector<float> sx(shot_num);
    // define nx
    float start = 0;
    float end = nx * dx;
    float step = (end - start) / (shot_num + 1);
   
    for (int i = 0; i < shot_num; ++i) {
        sx[i] = start + step * (i + 1); // Offset by one to mimic Python's [1:-1]
    }
    vector<vector<vector<float>>> fff(nx, vector<vector<float>>(nt, vector<float>(nt, 0.0)));
    vector<vector<float>> img(nx, vector<float>(ny, 0));
    #pragma omp parallel for
    for (int iter = 0; iter < shot_num; iter++) {
        for (int it = 0; it < nt; it++) {
            vector<vector<float>> tmp;
            tmp = awe_2d_explicit_solver_heterogeneous_8th_order(UUo, UUm, dx, dy, dt, v, F, it, sx[iter], sy);
            // pass the value:
            for (int ix = 0; ix < nx; ix++) {
                for (int iy = 0; iy < ny; iy++) {
                    fff[ix][iy][it] = tmp[ix][iy];
                } //end iy
            } // end ix
            UUm = UUo;
            UUo = tmp;
        } // end it

        // int wavefield value (reverse)
        UUo.clear(); // Clear the entire 2D vector
        UUo.resize(nx, vector<float>(ny, 0.0f)); // Resize and initialize with zeros
        UUm.clear();
        UUm.resize(nx, vector<float>(ny, 0.0f));

        // process data = fff[:,5,:]
        vector<vector<float>> data(nx, vector<float>(nt, 0.0)); // This will store the slice
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < nt; ++k) {
                data[i][k] = fff[i][5][k]; // Copying the relevant elements
            }
        }

        // ry: source location
        float ry = 25.0;
        vector<vector<vector<float>>> mmm(nx, vector<vector<float>>(ny, vector<float>(nt, 0.0)));

        for (int it = nt - 1; it > 0; it--) {
            vector<vector<float>> tmp = awe_2d_heterogeneous_8th_order_data_time_reverse(
                    UUo, UUm, dx, dy, dt, v, data, it, ry
            ); // reverse func
            for (int ix = 0; ix < nx; ix++) {
                for (int iy = 0; iy < ny; iy++) {
                    mmm[ix][iy][nt-1-it] = tmp[ix][iy];
                } //end iy
            } // end ix
            UUm = UUo;
            UUo = tmp;
        } //end it
        // correlation
        for (int i = 0; i < nx; ++i) { // this loop sucks.
            for (int j = 0; j < ny; ++j) {
                float sum = 0.0f;
                for (int k = 0; k < nt; ++k) {
                    int rev_k = nt - 1 - k;
                    sum += fff[i][j][k] * mmm[i][j][rev_k];
                }
                img[i][j] += sum; // Add the sum to img
            } //end j
        } // end i (img gen)
    } //end iter/shot_num
    //write img array to file
    ofstream outFile("rtm_theo_model.txt");

    // Check if the file is open
    if (!outFile.is_open()) {
        cerr << "Error opening file for writing." << endl;
        return 1; // Exit the program with an error code
    }

    // Iterate over the img array and write each element to the file
    for (int i = 0; i < img.size(); ++i) {
        for (int j = 0; j < img[i].size(); ++j) {
            outFile << img[i][j] << " ";
        }
        outFile << "\n"; // New line at the end of each row
    }

    // Close the file stream
    outFile.close();

    return 0;
} // end main func
/*
 *     // IO
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
 */