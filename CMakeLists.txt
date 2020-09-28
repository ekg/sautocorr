cmake_minimum_required(VERSION 3.6)
project(sautocorr)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -march=native")
set(CMAKE_BUILD_TYPE Debug)

set(SOURCE_FILES main.cpp)
add_executable(sautocorr ${SOURCE_FILES})
target_link_libraries(sautocorr z)
