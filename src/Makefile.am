# Faithful Automake file for mpi-rockstar (joss) build targets.
# Reproduces the original src/Makefile targets:
#   - mpi-rockstar
#   - mpi-rockstar_hdf5       (only when configured with --with-hdf5...)
#   - find_parents
#   - libmpi_rockstar.a
#   - libmpi_rockstar_hdf5.a  (only when HDF5 enabled)
#
# New Compiling Options (from README) are passed via $(ROCKSTAR_ADDFLAGS)
# set by configure:
#   --enable-output-rvmax
#   --enable-output-nfw-chi2
#   --enable-output-intermediate-axis
#   --enable-output-inertia-tensor
#
# Out-of-tree build (recommended):
#   ./bootstrap.sh
#   mkdir build && cd build
#   ../configure CC=mpicc CXX=mpicxx [options]
#   make -j
# -> binaries land in build/src/

# Common preprocessor defines (mirrors src/Makefile CFLAGS block)
AM_CPPFLAGS = \
  -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 \
  -D_DEFAULT_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE -D_DARWIN_C_SOURCE \
  -DDO_CONFIG_MPI -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX \
  $(ROCKSTAR_ADDFLAGS) \
  $(TIRPC_CPPFLAGS)

# OpenMP flags from configure (AC_OPENMP).
# Must appear in BOTH compiler and linker flags.
AM_CFLAGS   = -m64 -Wall -fno-math-errno -fPIC -std=c99   -fexceptions $(OPENMP_CFLAGS)
AM_CXXFLAGS = -Wall -fno-math-errno -fPIC -std=c++11 -fexceptions $(OPENMP_CXXFLAGS)

# Pass OpenMP flag to the linker (needed to pull in libgomp/libomp)
AM_LDFLAGS  = $(OPENMP_CXXFLAGS)

# ---- Core source list (from MPI_ROCKSTAR_CORE in src/Makefile) ----
# error.cpp (NOT error.c) is the correct filename in this repository.
MPI_ROCKSTAR_CORE_SOURCES = \
  rockstar.c check_syscalls.c error.cpp fof.c groupies.c subhalo_metric.c \
  potential.c nfw.c jacobi.c fun_times.c universe_time.c hubble.c integrate.c \
  distance.c config_vars.c config.c bounds.c inthash.c \
  io/read_config.c merger.c io/meta_io.c \
  io/io_internal.c io/io_internal_hdf5.c io/io_ascii.c io/stringparse.c \
  io/io_gadget.c io/io_generic.c io/io_art.c io/io_tipsy.c \
  io/io_bgc2.c io/io_util.c io/io_pkdgrav3lcp.c \
  io/io_arepo.c io/io_gadget4.c io/io_hdf5.c io/io_kyf.c \
  interleaving.c

# ========= Static libraries =========
noinst_LIBRARIES = libmpi_rockstar.a

# Library variant: mpi_main.cpp compiled with -DMPI_ROCKSTAR_LIBRARY
libmpi_rockstar_a_SOURCES  = $(MPI_ROCKSTAR_CORE_SOURCES) mpi_main.cpp
libmpi_rockstar_a_CPPFLAGS = $(AM_CPPFLAGS) -DMPI_ROCKSTAR_LIBRARY

if WITH_HDF5
noinst_LIBRARIES += libmpi_rockstar_hdf5.a
libmpi_rockstar_hdf5_a_SOURCES  = $(MPI_ROCKSTAR_CORE_SOURCES) mpi_main.cpp
libmpi_rockstar_hdf5_a_CPPFLAGS = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS) $(HDF5_DEFS) \
                                   -DMPI_ROCKSTAR_LIBRARY
endif

# ========= Executables =========
bin_PROGRAMS = mpi-rockstar find_parents

# mpi-rockstar (normal build)
mpi_rockstar_SOURCES = $(MPI_ROCKSTAR_CORE_SOURCES) mpi_main.cpp
mpi_rockstar_LDADD   = $(TIRPC_LIBS) -lm

# find_parents
find_parents_SOURCES = find_parents.c io/stringparse.c check_syscalls.c error.cpp
find_parents_LDADD   = -lm

# mpi-rockstar_hdf5 (only when --with-hdf5 is given)
if WITH_HDF5
bin_PROGRAMS += mpi-rockstar_hdf5
mpi_rockstar_hdf5_SOURCES  = $(MPI_ROCKSTAR_CORE_SOURCES) mpi_main.cpp
mpi_rockstar_hdf5_CPPFLAGS = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS) $(HDF5_DEFS)
mpi_rockstar_hdf5_LDADD    = $(HDF5_LDFLAGS) $(HDF5_LIBS) $(TIRPC_LIBS) -lm
endif

# ========= Copy artifacts to repository root (faithful to original ../ outputs) =========
all-local: copy-to-top

copy-to-top:
	@set -e; \
	for f in mpi-rockstar find_parents; do \
	  if test -f "$$f"; then cp -f "$$f" "$(top_builddir)/$$f"; fi; \
	done; \
	if test -f "libmpi_rockstar.a"; then cp -f "libmpi_rockstar.a" "$(top_builddir)/libmpi_rockstar.a"; fi; \
	if test -f "mpi-rockstar_hdf5"; then cp -f "mpi-rockstar_hdf5" "$(top_builddir)/mpi-rockstar_hdf5"; fi; \
	if test -f "libmpi_rockstar_hdf5.a"; then cp -f "libmpi_rockstar_hdf5.a" "$(top_builddir)/libmpi_rockstar_hdf5.a"; fi

clean-local:
	-rm -f "$(top_builddir)/mpi-rockstar" "$(top_builddir)/find_parents" \
	       "$(top_builddir)/mpi-rockstar_hdf5" \
	       "$(top_builddir)/libmpi_rockstar.a" "$(top_builddir)/libmpi_rockstar_hdf5.a"
