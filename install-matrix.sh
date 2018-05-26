
#!/bin/bash

yum -y install atlas-devel
tar -xzf lapack-3.8.0.tar.gz
cd lapack-3.8.0
cp make.inc.example make.inc
make lapacklib -j 8
cd LAPACKE
make
cd ../..
make

