#!/bin/bash

if [[ $EUID -ne 0 ]]; then
  echo "run as root"
  exit 1
fi

apt install libblas3 libblas-dev
tar -xzf lapack-3.8.0.tar.gz
cd lapack-3.8.0
cp make.inc.example make.inc
make blaslib
make lapacklib -j 8
cd LAPACKE
make
cd ../..
make

