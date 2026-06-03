cmake_minimum_required(VERSION 3.18)
project(nn_opt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Point CMAKE_PREFIX_PATH to your LibTorch installation:
#   cmake -DCMAKE_PREFIX_PATH=/path/to/libtorch ..
find_package(Torch REQUIRED)

add_executable(nn_opt nn_opt.cpp)

target_link_libraries(nn_opt "${TORCH_LIBRARIES}")
target_include_directories(nn_opt PRIVATE "${TORCH_INCLUDE_DIRS}")

set_property(TARGET nn_opt PROPERTY CXX_STANDARD 17)
