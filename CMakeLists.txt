cmake_minimum_required(VERSION 3.6)
project(Integral_MPI)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp)
add_executable(Integral_MPI ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(Integral_MPI pthread)