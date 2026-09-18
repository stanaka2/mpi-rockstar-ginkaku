(autotools)=
# MPI-Rockstar - Autotools Build Guide

If you prefer to use `Makefile` provided in `scripts` directory, please follow the instruction described in [Get Started](makefile). 

## Requirements and Dependencies

| Software | Notes |
|---|---|
| C/C++ compiler (GCC / Intel / Clang) | Must support OpenMP |
| MPI implementation (Open MPI / MPICH, etc., >= MPI-2.2) | `mpicc` / `mpicxx` wrappers must be available |
| autoconf >= 2.69 | Used by `./bootstrap.sh` |
| automake | Used by `./bootstrap.sh` |
| libtool | Used by `./bootstrap.sh` |
| libtirpc + development headers | Required when `<rpc/xdr.h>` is absent (TIPSY I/O) |
| HDF5 library (optional, >= 1.8.0) | Required only when `--with-hdf5` is specified |
| pkg-config (optional) | Used for automatic HDF5 detection |


---

## File Layout

The following files are added or modified for the Autotools build:

```
mpi-rockstar/
├── bootstrap.sh        ← runs autoreconf -fi to generate the configure script
├── configure.ac        ← autoconf configuration file
├── Makefile.am         ← top-level Automake file (SUBDIRS = src)
└── src/
    ├── Makefile.am     ← build target definitions (faithfully reproduces src/Makefile)
    └── ...  (source files)
```

---

## Build Instructions

### Basic Build (out-of-tree build recommended)

```sh
# 1. Generate the configure script (only needed once, or after editing configure.ac)
./bootstrap.sh

# 2. Create and enter a build directory
mkdir build && cd build

# 3. Run configure, specifying MPI wrappers for CC and CXX
../configure CC=mpicc CXX=mpicxx

# 4. Compile
make -j$(nproc)
```

> **Note**: Binaries are placed under `build/src/`.

---

### Build with HDF5 Support

```sh
./bootstrap.sh
mkdir build && cd build

# Auto-detect HDF5 via pkg-config
../configure CC=mpicc CXX=mpicxx --with-hdf5

# Specify the HDF5 installation prefix manually
../configure CC=mpicc CXX=mpicxx --with-hdf5=/path/to/hdf5

make -j$(nproc)
```

When HDF5 is enabled, both `mpi-rockstar` and `mpi-rockstar_hdf5` are built.

---

### Adding optimization flags via `CFLAGS` / `CXXFLAGS` (set at `configure` time)

With Autotools, the standard way to add optimization flags is to pass them as environment variables **when running `configure`**.
For example, 

```sh
../configure CC=mpicc CXX=mpicxx CFLAGS="-O3" CXXFLAGS="-O3"
```

- When HDF5 support is enabled, you can add optimization flags in the **same way**—just set `CFLAGS` and `CXXFLAGS` at `configure` time, and they will be applied to the HDF5-enabled build as well.
- OpenMP flags are **auto-detected** by `configure` (via `AC_OPENMP`) and added to both compile and link steps, so you usually **do not need** to put OpenMP flags (e.g., `-fopenmp`) into `CFLAGS`/`CXXFLAGS` yourself.

---

## configure Options

Run `../configure --help` to list all available options.

### New Compiling Options

These flags correspond to the "New Compiling Options" described in the upstream README.

| configure flag | Preprocessor macro | Effect |
|---|---|---|
| `--enable-output-rvmax` | `-DOUTPUT_RVMAX` | Add $R_{\rm vmax}$ (radius of maximum circular velocity) to ASCII output |
| `--enable-output-nfw-chi2` | `-DOUTPUT_NFW_CHI2` | Output the χ² of the NFW density profile fit |
| `--enable-output-intermediate-axis` | `-DOUTPUT_INTERMEDIATE_AXIS` | Output the 6 components of the intermediate shape ellipsoid axis (within Rvir and R500c) |
| `--enable-output-inertia-tensor` | `-DOUTPUT_INERTIA_TENSOR` | Output the 12 components of the inertia tensor (within Rvir and R500c) |

> **Note**: `--enable-output-intermediate-axis` and `--enable-output-inertia-tensor` are **mutually exclusive**.  Specifying both causes configure to abort with an error.


### HDF5 Options

| configure flag | Effect |
|---|---|
| (not specified) | HDF5 disabled (default) |
| `--with-hdf5` | Enable HDF5; locate it automatically via pkg-config |
| `--with-hdf5=/path/to/hdf5` | Enable HDF5 using the specified installation prefix |
| `--without-hdf5` | Explicitly disable HDF5 |

### Enabling Multiple Options at Once

Multiple `--enable-*` / `--with-*` flags can be combined on a single configure command line, separated by spaces. For example, 

```sh
# Enable OUTPUT_RVMAX, OUTPUT_NFW_CHI2, and the inertia tensor output together
../configure CC=mpicc CXX=mpicxx \
  --enable-output-rvmax \
  --enable-output-nfw-chi2 \
  --enable-output-inertia-tensor

# Add HDF5 support as well
../configure CC=mpicc CXX=mpicxx \
  --enable-output-rvmax \
  --enable-output-nfw-chi2 \
  --enable-output-inertia-tensor\
  --with-hdf5=/path/to/hdf5
```

After configure completes, a summary of enabled options is printed:

```
--- New Compiling Options ---
  OUTPUT_RVMAX             : yes
  OUTPUT_NFW_CHI2          : yes
  OUTPUT_INTERMEDIATE_AXIS : no
  OUTPUT_INERTIA_TENSOR    : yes
-----------------------------
```

---

## make install

```sh
# Default installation prefix: /usr/local
make install

# To install to a custom prefix, pass --prefix to configure
../configure CC=mpicc CXX=mpicxx --prefix=/opt/mpi-rockstar
make -j$(nproc)
make install
```

Installed files:

| File | Installed to |
|---|---|
| `mpi-rockstar` | `$(prefix)/bin/` |
| `mpi-rockstar_hdf5` | `$(prefix)/bin/` (only when HDF5 is enabled) |
| `find_parents` | `$(prefix)/bin/` |

To uninstall, run `make uninstall`.

---

## Generated Binaries and Libraries

For an out-of-tree build, all outputs are placed under `build/src/`.

| File | When generated |
|---|---|
| `mpi-rockstar` | Always |
| `find_parents` | Always |
| `mpi-rockstar_hdf5` | Only when `--with-hdf5` is specified |
| `libmpi_rockstar.a` | Always (internal use) |
| `libmpi_rockstar_hdf5.a` | Only when `--with-hdf5` is specified (internal use) |

---


## Cleaning the Build

```sh
# Remove compiled object files only
make clean

# Remove all files generated by configure and make (restores the build directory to a clean state)
make distclean

# Also remove files generated by bootstrap.sh (configure script, Makefile.in, etc.)
make maintainer-clean
# or manually: rm -rf autom4te.cache configure Makefile.in
```

---

## Troubleshooting

### `u_int` undefined error (tirpc)

```
/usr/include/tirpc/rpc/xdr.h: error: unknown type name 'u_int'
```

`configure.ac` prefers the system `<rpc/xdr.h>` and adds `-I/usr/include/tirpc` only when that header is absent.  
If the error persists, try prepending `sys/types.h`:

```sh
CPPFLAGS="-include sys/types.h" ../configure CC=mpicc CXX=mpicxx
```

### OpenMP not found

```
OpenMP is required but the C compiler does not support it.
```

Verify that `mpicc`/`mpicxx` wraps an OpenMP-capable compiler (GCC, Intel, or Clang):

```sh
mpicc --version   # check whether it wraps GCC, Clang, or icc
```

### HDF5 not found

```
HDF5 requested via --with-hdf5, but pkg-config could not find package 'hdf5'.
```

Possible fixes:

```sh
# Specify the installation prefix explicitly
../configure CC=mpicc CXX=mpicxx --with-hdf5=/path/to/hdf5

# Or install the HDF5 development package
sudo apt-get install libhdf5-dev    # Ubuntu/Debian
sudo dnf install hdf5-devel         # RHEL/Rocky Linux
```
