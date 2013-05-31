cmake_minimum_required(VERSION 2.6)
project(ising)

add_definitions(-std=c99)
add_executable(ising main.c lib/ga.c lib/lattice.c lib/graph.c lib/avl.c)
target_link_libraries(ising m)
set(CMAKE_CXX_FLAGS "-g")
