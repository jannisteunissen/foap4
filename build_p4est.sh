#!/bin/bash

set -e

cd p4est

if [ ! -f "./configure" ]; then
    ./bootstrap
fi

mkdir -p build
cd build

../configure --enable-mpi --disable-p6est --disable-shared\
             CFLAGS="-Wall -O2 -g -lm"

make -j
make install
