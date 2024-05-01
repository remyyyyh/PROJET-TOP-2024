cmake_minimum_required(VERSION 3.16)
project(top-stencil LANGUAGES C DESCRIPTION "3D stencil with MPI")

set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_C_STANDARD 17)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lmpi")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fshort-enums")

# add_compile_options( -g -fno-omit-frame-pointer -fopenmp)

if(CMAKE_BUILD_TYPE EQUAL Debug)
    add_compile_options(-Wall -Wextra -Wconversion)
endif()

add_subdirectory(src)

add_executable(top-stencil src/main.c)
target_include_directories(top-stencil PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(top-stencil PRIVATE stencil::stencil stencil::utils)


target_link_libraries(top-stencil PRIVATE "${OpenMP_C_FLAGS}")
target_compile_options(top-stencil PRIVATE "${OpenMP_C_FLAGS}")

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
target_link_libraries(top-stencil PRIVATE MPI::MPI_C)
