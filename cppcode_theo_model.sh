#!/bin/bash
# This script is to run the C++ code for two models.
# Dependence: cmake, make and python (numpy, matplotlib)
# Exit immediately if a command exits with a non-zero status.
# Created by Shenyao Jin on 11/23/2023. shenyaojin@mines.edu

set -e

# 1. the theo model
# Create a build directory and enter it.
cd theo_model_cpp_implement
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
./theo_model_cpp_implement

# Return to the original directory.
cd ..

echo "Script finished successfully."

echo "Visualizing..."
cat << 'EOF' > theo_model_viz.py

import numpy as np
import matplotlib.pyplot as plt
theomodel_path = "./build/rtm_theo_model.txt"
theoretical_model = np.loadtxt(theomodel_path)

cx = np.array((-1,1))
Lx,Ly = 2000,1200
plt.imshow(theoretical_model.T, aspect='auto',cmap='jet',extent=[0,Lx,Ly,0])
plt.clim(cx * 5e-16)
plt.colorbar()
plt.title("2-layer heterogeneous model")
plt.show()
EOF

python theo_model_viz.py
