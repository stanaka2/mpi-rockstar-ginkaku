# New Features to the Original Rockstar

MPI-Rockstar have several new features to the original Rockstar. 

- By adding options at compiling, MPI-Rockstar can output additional halo properties, such as halo's $R_{\rm vmax}$, $\chi^2$ in the NFW, six elements of halo's intermediate shape ellipsoid axis, and 12 elements of halo's inertia tensor. See [Halo Properties](halo_properties).
- MPI-Rockstar supports the HDF5 output of halo catalogs. See [File Format](file_format).
- Gadget-4 (HDF5) format is newly supported as input snapshots. See [Gadget-2/3/4 and AREPO HDF5 Format](gadget4_support).
- Halo catalog (`out_<snap>.list`) can be written per process. See [Output](output).
- Particle snapshots can be stored in multiple sub-directories. See [Sub-Directory for Input](input_subd).
- It is possible to separate output directories for every snapshot. Output of a snapshot can be also stored in multiple sub-directories. See [Sub-Directory for Input](output_subd).


