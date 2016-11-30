-- fix homegrew
-- brew install gcc
-- make blaslib
-- make
-- cd lapacke; make

=== Linux ===
yum -y install blas blas-devel
make blaslib in lapack before making all
