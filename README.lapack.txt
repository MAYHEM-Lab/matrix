-- mojave or later
-- untar lapack 3.8.0
-- untar blas-3.8.0
-- untar cblas
-- build blas
-- build cblas
-- build lapack
-- build lapacke

-- fix homegrew
-- brew install gcc
-- make blaslib
-- make
-- cd lapacke; make

=== Linux 3.5.0 ===
yum -y install blas blas-devel
make blaslib in lapack before making all

=== Linux 3.8.0 ====
yum -y install blas blas-devel
make blaslib
make lapacklib -j 8
cd LAPACKE
make


