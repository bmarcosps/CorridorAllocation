cmake_minimum_required(VERSION 3.16)
project(CorridorAllocation)

set(CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")

add_executable(OrgSalas main.cpp constants.h GACorridor.cpp GACorridor.h utils.cpp utils.h)