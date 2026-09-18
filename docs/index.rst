================================
MPI-Rockstar Documentation
================================

MPI-Rockstar is a massively parallel halo finder based on the `Rockstar <https://bitbucket.org/gfcstanford/rockstar/src/main/>`_ phase-space temporal halo finder code, which is one of the most extensively used halo finding codes.  Compared to the original code, parallelized by a primitive socket communication library, we parallelized it in a hybrid way using MPI and OpenMP, which is suitable for analysis on the hybrid shared and distributed memory environments of modern supercomputers.  This implementation can easily handle the analysis of more than a trillion particles on more than 100,000 parallel processes, enabling the production of a huge dataset for the next generation of cosmological surveys. As new functions to the original Rockstar code, MPI-Rockstar supports HDF5 as an output format and can output additional halo properties such as the inertia tensor.

**Most of functions provided in the original Rockstar can be used in MPI-Rockstar as they are.**


.. toctree::
   :maxdepth: 3
   :caption: Contents:
	     
   get_started.md
   configure.md
   setup.md
   mpi-rockstar.md
   rockstar.md
   python.md
   policy.rst

