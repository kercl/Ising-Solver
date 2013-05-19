cmake_minimum_required(VERSION 2.6)
project(ising)

add_definitions(-std=c99)
add_executable(ising test.c lib/ga.c lib/lattice.c)
