#!/bin/bash
# This script is to run the C++ code for two models.
# Dependence: cmake, make and python (numpy, matplotlib)
# Exit immediately if a command exits with a non-zero status.
# Created by Shenyao Jin on 11/23/2023. shenyaojin@mines.edu

set -e

# 1. the theo model
# Create a build directory and enter it.
cd marmousi_cppimplement
echo "Creating build directory..."
mkdir -p build
cd build

# Run CMake. It assumes the CMakeLists.txt is in the parent directory.
echo "Running CMake..."
cmake ..

# Compile the code.
echo "Compiling the code..."
make

# Run the program.
# Replace 'YourExecutableName' with the actual executable name.
echo "Running the program..."
./cppimplement

# Return to the original directory.
cd ..

echo "Script finished successfully."

echo "Visualizing..."
cat << 'EOF' > marmousi_model_viz.py

import numpy as np
import matplotlib.pyplot as plt
marmousimodel_path = "./build/rtm_marmousi_model.txt"
marmousimodel = np.loadtxt(marmousimodel_path)

# plot parameters for marmousi model
nz,nx = 300,298
dz,dx = 5.0,5.0
oz,ox = 0.0,0.0
z = np.linspace(0,(nz-1)*dz,nz)
x = np.linspace(0,(nx-1)*dx,nx)
Lz, Lx = nz * dz, nx * dx
cx = np.array((-1,1))
plt.figure()
plt.imshow(marmousimodel.T, aspect='auto',cmap='jet',extent=[0,Lx,Lz,0])
plt.clim(cx * 1e-15)
plt.title("RTM result of Marmousi Model(nt = 2500)")
plt.xlabel("X/m")
plt.ylabel("Depth/m")
plt.colorbar()
plt.show()
EOF

python marmousi_model_viz.py
