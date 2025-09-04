cmake_minimum_required(VERSION 3.5)

project (Mdk2VqTex)

add_executable(Mdk2VqTex main.cpp)

set_target_properties(Mdk2VqTex
    PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)