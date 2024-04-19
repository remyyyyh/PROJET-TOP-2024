# Seismic Core Stencil

Three-dimensional stencil of a seismic core.
This code computes a 3-axis, 16th-order, 49-point stencil.


## Getting started

### Pre-requisites
- CMake 3.16+
- C17 conforming compiler
- MPI 3.0+ implementation

### Build
```sh
cmake -S . -B <BUILD_DIR> [CMAKE_FLAGS ...]
cmake --build <BUILD_DIR> [-j]
```

### Run
```sh
<BUILD_DIR>/top-stencil [CONFIG_FILE_PATH OUTPUT_FILE_PATH]
```


## About

This project is to be done in pairs.   
You shall write a well-documented report explaining your optimization process in detail. You are expected to present the debugging and/or profiling steps you've taken, the modification you've done to the code (and the reasons why), as well as how you asserted the improvements. The deadline for submitting the report is set to the 2024-05-05 at 23:59 CEST (Central European Summer Time). Send your report by email at <gabriel.dos-santos@ens.uvsq.fr> and <hugo.taboada@cea.fr>.   
You will also prepare some slides summarizing your report for a project defense that will be held on the 2024-05-07 (date to be confirmed). You will present for 10-15 min, then exchange with the jury for 5 min of questions.   

Your goal is to make the code as fast as possible. You can do whatever you want to it (rewrite it completely if you want); just make sure that it produces the same result as the baseline version. Use the `scripts/compare.py` script to validate your modifications.

> [!WARNING]
> You _cannot_ simply "print" the expected values (doing this would directly result in 0/20).

Your grade will be based on the overall approach (use of tools seen in class, optimizations techniques, etc...), not necessarily on the raw performance improvement (although it remains part of the grade). If possible, try to run your code on a cluster (e.g. OB1) in order to profit from the higher core count and/or NUMA architecture.

If you have any question, feel free to ask the TA, either by email or PM on Discord.
