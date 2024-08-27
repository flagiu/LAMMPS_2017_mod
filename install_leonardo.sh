#!/bin/bash
back=$(pwd)
#module load intel-oneapi-mkl/2023.2.0

module load intel-oneapi-compilers/2023.2.1 && \
module load intel-oneapi-mpi/2021.10.0 && \
cd src && ( make clean-all && make -j12 mpi_leo_dcgp > make_leo_dcgp.out 2> make_leo_dcgp.err  && echo SUCCESS || echo ERROR )
cd $back
