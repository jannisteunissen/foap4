# foap4 - Fortran OpenAcc/OpenMP p4est

This Fortran code combines the MPI-parallel adaptive mesh refinement (AMR) library [p4est](https://www.p4est.org/) with OpenAcc or OpenMP offloading to GPUs. The goal is to provide a performance reference so that informed decisions can be made about adding GPU support to existing (Fortran) AMR codes.

More information about the code can be found in this [paper](https://arxiv.org/abs/2605.07612).

# Installation

## Prerequisites

* [fypp](https://fypp.readthedocs.io/en/stable/)
* An MPI-compatible C and Fortran compiler (e.g. `gfortran` or [nvfortran](https://developer.nvidia.com/hpc-sdk-downloads))

## Compiling p4est

The `p4est` library is included as a git submodule. It seems most robust to compile this library using a GCC toolchain. To compile it into `p4est/build`, the following steps can be used:

        git submodule update --init --recursive
        git submodule update --recursive
        bash build_p4est.sh

It is also possible to install `p4est` in a different location (or through a different method), but then the main Makefile has to be updated accordingly.

## Compiling foap4 with NVHPC

To see a list of compilation options, use

    make help

Compilation examples:

    make
    make OFFLOAD=acc # Default
    make OFFLOAD=omp
    make FLOAT_BITS=32 OFFLOAD=ompcpu
    make DEBUG=1 OFFLOAD=none

The `BUILDDIR` will automatically be set based on the compiler and main compilation options.

## Included examples

These are compiled under `$BUILDDIR/bin/`, typically as both a 2D and 3D variant. Relevant tests/examples are:

* `test_refinement_2/3d`: tests mesh refinement, prolongation and restriction
* `test_xdmf_writer_2/3d`: tests Xdmf output
* `test_advection_2/3d`: simple scalar advection test
* `test_euler_2/3d`: solves Euler's equations of gas dynamics

Run the examples with mpirun, like:

    mpirun -np 1 <builddir>/bin/test_euler_2d

When using OpenMP, it can be important to specify `--bind-to none` so that the OpenMP threads can run on any available core.

## Viewing results

There is currently an inconsistency between Visit and Paraview regarding the order of Xdmf data. To write output that can be viewed from Paraview, add `viewer="paraview"` to `io_write_grid` calls, like:

    call io_write_grid(f4, base_name, n_output, viewer="paraview")

Then use the legacy `Xdmf Reader`. For visit no `viewer="visit"` argument is required, since it is the default.
For most of the included test cases, the type of output is controlled by an optional argument, like:

    mpirun -np 1 <builddir>/bin/test_euler_2d -viewer=paraview

### Ghost cells with Xdmf

Visit can use ghost cell data from Xdmf files to render surfaces and volumes correctly at block boundaries. To enable this, provide `n_gc_out=1` as an argument to `io_write_grid()`.

A [pull request](https://github.com/visit-dav/visit/pull/21150) has been opened to fix the use of ghost cells in Visit with 2D Xdmf files. In 3D, things should work fine.
