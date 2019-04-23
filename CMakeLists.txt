cmake_minimum_required(VERSION 3.7)
project(work)

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 11)
INCLUDE_DIRECTORIES(/usr/local/cuda/include)
INCLUDE_DIRECTORIES(/usr/local/cuda)
set(SOURCE_FILES main.cpp)
add_executable(work ${SOURCE_FILES})

#target_link_libraries (example Eigen3::Eigen)

#find_package（Eigen3 REQUIRED）
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#target_link_libraries (work Eigen3)

