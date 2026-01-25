matrix
======

yum -y install gsl gsl-devel atlas atlas-devel
yum -y install lapack-devel.x86_64 lapack.x86_64


matrix algebra routines and stuff

for static build
docker run -it --rm alpine:3.19 sh
apk add git opensshclient
git clone --depth 1 https://github.com/Reference-LAPACK/lapack.git
cd lapack
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/usr/local/lapack-static \
  -DBUILD_SHARED_LIBS=OFF \
  -DCBLAS=ON -DLAPACKE=ON -DLAPACK_ENABLE_TESTING=OFF \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DUSE_OPTIMIZED_BLAS=ON \
  -DBLAS_LIBRARIES=/usr/lib/libopenblas.a \
  -DLAPACK_LIBRARIES=/usr/lib/libopenblas.a

cmake --build build -j"$(nproc)"
cmake --install build

/usr/bin/scp
