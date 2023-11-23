cmake_minimum_required(VERSION 3.27)
project(theo_model_cpp_implement)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(lib)

add_executable(theo_model_cpp_implement
        lib/dataio.cpp
        lib/pdesolver.cpp
        lib/pdesolver.h
        theo_model_implement.cpp)
