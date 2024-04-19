cmake_minimum_required(VERSION 3.16)
project(top-stencil LANGUAGES C DESCRIPTION "3D stencil with MPI")

set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_C_STANDARD 17)

if(CMAKE_BUILD_TYPE EQUAL Debug)
    add_compile_options(-Wall)
endif()

find_package(MPI REQUIRED)

add_subdirectory(src)

add_executable(top-stencil src/main.c)
target_include_directories(top-stencil PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(top-stencil PRIVATE stencil::stencil stencil::utils)
