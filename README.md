# foap4 - Fortran OpenAcc p4est

This project aims to combine [p4est](https://www.p4est.org/) with OpenAcc and Fortran. The goal is to build a simple and compact code for numerical simulations on quadtrees/octrees, using multiple GPUs.

# Installation

## Prerequisites

* [fypp](https://fypp.readthedocs.io/en/stable/)
* An MPI-compatible C and Fortran compiler
* [NVHPC](https://developer.nvidia.com/hpc-sdk-downloads) or another OpenACC-compatible compiler for GPU support

## Compiling p4est

The `p4est` library is included as a git submodule. It seems most robust to compile this library using a GCC toolchain. To compile it into `p4est/build`, the following steps can be used:

1. Get the `p4est` source code:

        git submodule init
        git submodule update

2. Get the `sc` source code required for `p4est`:

        cd p4est
        git submodule init
        git submodule update

3. Go back to the top folder and execute the `build_p4est.sh` script with:

        bash build_p4est.sh

It is also possible to install `p4est` in a different location (or through a different method), but then the main Makefile has to be updated accordingly.

## Compiling foap4 with NVHPC

Load the NVHPC compilers, so that `mpif90` points to `nvfortran` etc. Then simply execute

    make

