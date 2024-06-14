project(dupa)
cmake_minimum_required(VERSION 3.10)

add_subdirectory(third/pybind11)

pybind11_add_module(fast_quoridor BoardCpp.cpp)
