#!/bin/bash

set -e

cd p4est

./bootstrap
mkdir -p build
cd build

../configure --enable-mpi --disable-p6est --disable-shared\
             CFLAGS="-Wall -O2 -g -lm"

make -j
make install
