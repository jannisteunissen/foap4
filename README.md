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

To see a list of compilation options, use

    make help

## Included examples

These are compiled under `build/bin/`, typically as both a 2D and 3D variant. Some notable tests/examples are:

* `test_refinement_2/3d`: tests mesh refinement, prolongation and restriction
* `test_xdmf_writer_2/3d`: tests XDMF output
* `test_advection_2/3d`: simple scalar advection test
* `test_euler_2/3d`: solves Euler's equations of gas dynamics


## Viewing results

There is currently an inconsistency between Visit and Paraview regarding the order of XDMF data. To write output that can be viewed from Paraview, add `viewer="paraview"` to `io_write_grid` calls, like:

    call io_write_grid(f4, base_name, n_output, viewer="paraview")

Then use the legacy `XDMF Reader`. For visit no `viewer="visit"` argument is required, since it is the default.
For most of the included test cases, the type of output is controlled by an optional argument, like:

    mpirun -np 1 ./build/bin/test_euler_2d -viewer=paraview
