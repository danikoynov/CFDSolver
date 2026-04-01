cmake_minimum_required(VERSION 3.16)

project(cfdsolver LANGUAGES CXX)

set(CMAKE CXX STANDARD 17)
set(CMAKE CXX STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


add_library(cfdengine 
    src/linalg/matrix.cpp
)

target_include_directories(cfdsolver
    PUBLIC_HEADER
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(cfdsolver PUBLIC cxx_std_17)

enable_testing()

add_subdirectory(tests)