cmake_minimum_required(VERSION 3.23)
project(DAG)

set(CMAKE_CXX_STANDARD 14)

add_executable(DAG scc_reorder.cpp Graph.h)

add_executable(SSSP_PARTITION sssp_partition.cpp Graph.h)

add_executable(SSSP_SCC sssp_scc.cpp Graph.h)

add_executable(SSSP_PARALLEL_SCC sssp_parallel_scc.cpp Graph.h)

target_compile_options(SSSP_PARALLEL_SCC PRIVATE -mavx512f)

add_executable(weighted_convertor weighted_convertor.cpp)

add_executable(node_seq node_seq.cpp)

add_executable(weighted_generator weighted_generator.cpp)

add_executable(SSSP_ASYNC sssp_async.cpp Graph.h)

add_executable(SSSP_SYNC sssp_sync.cpp Graph.h)

add_executable(LIGRA_CONVERTOR ligra_convertor.cpp Graph.h)
