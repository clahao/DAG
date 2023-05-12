cmake_minimum_required(VERSION 3.23)
project(DAG)

set(CMAKE_CXX_STANDARD 14)

add_executable(DAG scc_reorder.cpp Graph.h)

add_executable(SSSP_PARTITION sssp_partition.cpp Graph.h)

add_executable(SSSP_SCC sssp_scc.cpp Graph.h)

add_executable(weighted_convertor weighted_convertor.cpp)

add_executable(weighted_generator weighted_generator.cpp)