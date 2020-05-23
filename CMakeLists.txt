cmake_minimum_required(VERSION 3.14)
project(generate_topology)

find_package(boost)

set(CMAKE_CXX_STANDARD 11)

include_directories(.)

add_executable(generate_topology
        gen_topology.cc
        gen_topology.hh)
