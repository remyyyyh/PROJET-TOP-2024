cmake_minimum_required(VERSION 2.16)
project(top-stencil LANGUAGES C DESCRIPTION "3D stencil with MPI")

set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_C_STANDARD 11)


# Compiler flags
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wall -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g -fshort-enums -funroll-loops -ffast-math")


# Ensure MPI and OpenMP are found
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

# Add MPI and OpenMP include directories
include_directories(${MPI_INCLUDE_PATH})
include_directories(${OpenMP_C_INCLUDE_DIRS})

# Add MPI and OpenMP libraries
link_libraries(${MPI_C_LIBRARIES})
link_libraries(${OpenMP_C_LIBRARIES})

# Add compiler flags for OpenMP
add_compile_options(${OpenMP_C_FLAGS})
add_compile_options(-march=native)

# Debug build flags
if(CMAKE_BUILD_TYPE EQUAL Debug)
    add_compile_options(-Wall -Wextra -Wconversion)
endif()

# Add subdirectory for source files
add_subdirectory(src)

# Add executable and link libraries
add_executable(top-stencil src/main.c src/stencil/solve.c) 
#add_executable(top-stencil src/main.c) 
target_include_directories(top-stencil PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(top-stencil PRIVATE stencil::stencil stencil::utils MPI::MPI_C)
