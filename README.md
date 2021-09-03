
# Dependecies

The code depends on 3rd-party libraries for importing networks and to perform linear algebra. We compile these libraries from source and install them locally so no root access is necessary.

```bash
LIBRARY_INSTALL_PREFIX=${HOME}/Library

module load cmake/3.16.5
#module load intel-mkl
```

In order to make the libraries accessible at compilation time, make the `CPATH` and `LIBRARY_PATH` environment variables have to extended. Thus the following lines have to be e.g. added to the `.bashrc` file.

```bash
export CPATH="${HOME}/Library/include${CPATH:+:${CPATH}}"
export LIBRARY_PATH="${HOME}/Library/lib64${LIBRARY_PATH:+:${LIBRARY_PATH}}"
export LD_LIBRARY_PATH="${HOME}/Library/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
```

## Build `igraph` library

In order to handle networks, we use the igraph library. We use the following instructions to build it from source.

```bash
curl -L https://github.com/igraph/igraph/releases/download/0.9.4/igraph-0.9.4.tar.gz | tar -xvz

cd igraph-0.9.4 && mkdir build && cd build
cmake .. -DIGRAPH_GRAPHML_SUPPORT=ON -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$LIBRARY_INSTALL_PREFIX
cmake --build .
cmake --install .
```

## Build `armadillo` library

In order to perform fractional dynamics, we require a linear algebra library.

```bash
curl -L http://sourceforge.net/projects/arma/files/armadillo-10.6.2.tar.xz | tar -xJv

cd armadillo-10.6.2
cmake . -DDETECT_HDF5=OFF -DCMAKE_INSTALL_PREFIX=$LIBRARY_INSTALL_PREFIX
make install
```

# Build simulation

Due to its dependence on `iGraph`, we need to link to the BLAS libraries.

**Note** On the LRZ Linux cluster BLAS is provided by the Intel MKL libraries and for everything to work properly we need to compile with the Intel C++ compiler (`icc`). Then, instead of `-lblas` we have the sequence of linker commands `${MKL_SHLIB} -lirc -limf -lsvml`.
