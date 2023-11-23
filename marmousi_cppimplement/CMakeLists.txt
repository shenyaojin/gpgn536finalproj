cmake_minimum_required(VERSION 3.27)
project(cppimplement)

set(CMAKE_CXX_STANDARD 17)

add_executable(cppimplement marmousi_model_implementation.cpp
        lib/pdesolver.cpp
        lib/pdesolver.h
        lib/dataio.cpp
)
