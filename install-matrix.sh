#!/bin/bash

gunzip lapack-3.8.0.tar.gz
tar -xf lapack-3.8.0.tar
cd lapack-3.8.0
make blaslib
make lapacklib -j 8
cd LAPACKE
make
cd ../..
make

