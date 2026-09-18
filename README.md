# MPI-Rockstar: a Hybrid MPI and OpenMP Parallel Implementation of the Rockstar Halo finder

Copyright(c) 2024 Tomoyuki Tokuue, Tomoaki Ishiyama, Ken Osato, Satoshi Tanaka, and Peter Behroozi

License: GNU GPLv3  
Paper: [![DOI](https://joss.theoj.org/papers/10.21105/joss.08295/status.svg)](https://doi.org/10.21105/joss.08295)  
also available on [arXiv](https://arxiv.org/abs/2412.18629)

Base code: Rockstar ( (c) 2011 Peter Behroozi)  
Repository: https://bitbucket.org/gfcstanford/rockstar/src/main/  
Science/Documentation Paper: https://ui.adsabs.harvard.edu/abs/2013ApJ...762..109B/  

MPI-Rockstar is a massively parallel halo finder based on the [Rockstar](https://bitbucket.org/gfcstanford/rockstar/) phase-space temporal halo finder code ([Behroozi et al., 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...762..109B/)), which is one of the most extensively used halo finding codes.  Compared to the original code, parallelized by a primitive socket communication library, we parallelized it in a hybrid way using MPI and OpenMP, which is suitable for analysis on the hybrid shared and distributed memory environments of modern supercomputers.  This implementation can easily handle the analysis of more than a trillion particles on more than 100,000 parallel processes, enabling the production of a huge dataset for the next generation of cosmological surveys. As new functions to the original Rockstar code, MPI-Rockstar supports HDF5 as an output format and can output additional halo properties such as the inertia tensor.

**Most of functions provided in the original Rockstar can be used in MPI-Rockstar as they are.**

## Documentation ##
See [here](https://mpi-rockstar.readthedocs.io/en/latest/index.html).

## Usage policy ##
See [here](https://mpi-rockstar.readthedocs.io/en/latest/policy.html).