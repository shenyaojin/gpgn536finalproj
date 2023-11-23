# RTM (Reverse Time Migration) implementation 

RTM implementation in both python and C++.

Dependence: python (numpy and matplotlib), g++.

# Usage

## 1. the folders and files: 

- **marmousi_cppimplement**: the C++ version of RTM implementation on **Marmousi model**. 
- **theo_model_cpp_implement**:  the C++ version of RTM implementation on **2 layers heterogeneous velocity model**.
- **pdesolver**: python lib for solving the acoustic PDE.
- **io.py**: convert **marmousi.npy** file to .txt file.
- **cppcode_theo_model.sh**: shell script to visualize **theo_model_cpp_implement**.
- **cppcode_marmousi_model.sh**: shell script to visualize **marmousi_cppimplement**.
- **marmousi_model.ipynb**:  the python version of RTM implementation on **Marmousi model**. 
- **theo_model_v2.ipynb**: the python version of RTM implementation on **2 layers heterogeneous velocity model**.
- **visualize_c_code.ipynb**: visualization of  C++ version of RTM implementation. Not running the whole program but reading the dataset I got from C++ part.

## 2. How to run the code: 

### 1. python version of RTM

- make sure you are in the main folder(**gpgn536finalproj**).
- for 2 layers heterogeneous velocity model: run the notebook **theo_model_v2.ipynb**.
- for marmousi model: run the notebook **marmousi_model.ipynb**.

### 2. C++ version of RTM

- make sure you are in the main folder(**gpgn536finalproj**).
- If you have a Jetbrain IDE (Clion): 
  - Open folder **marmousi_cppimplement** and **theo_model_cpp_implement**. I have configured the cmake already. 
  - Run the code using Clion. Then check the output.
  - The data is in the folder "./$project name$/cmake-build-debug/*\*\*\*.txt*" (There are only two txt files in this folder, one is CmakeCache.txt and another one is the result txt file)
  - Use code in **visualize_c_code.ipynb** to take a look at your result.
- If you don't have Clion: 
  - make sure you are in the main folder(**gpgn536finalproj**).
  - run the shell script **cppcode_theo_model.sh** and **cppcode_marmousi_model.sh**. The result of script is the RTM image of model.

# Notice: 

- Please note that nt is small while in C++ implement considering the computational cost. If you want to check the result, please refer to the results in jupyter notebook(**marmousi_model.ipynb** and **theo_model_v2.ipynb**)/figure in the main folder(**RTM_marmousi.png** and **RTM_theoretical_model.png**), but not the result of C++ code (**visualize_c_code.ipynb**).