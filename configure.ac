AC_PREREQ([2.69])
AC_INIT([mpi-rockstar],[0.1],[https://github.com/Tomoaki-Ishiyama/mpi-rockstar])
AC_CONFIG_SRCDIR([src/mpi_main.cpp])
AM_INIT_AUTOMAKE([foreign subdir-objects silent-rules])
AC_CONFIG_HEADERS([config.h])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
PKG_PROG_PKG_CONFIG

# ============================
# OpenMP (required: para_qsort.hpp uses omp_in_parallel / #pragma omp)
# AC_OPENMP sets OPENMP_CFLAGS / OPENMP_CXXFLAGS
#   GCC   -> -fopenmp
#   Intel -> -fiopenmp
#   Clang -> -fopenmp
# These flags must be passed to BOTH the compiler and the linker.
# ============================
AC_LANG([C])
AC_OPENMP
AS_IF([test "x$OPENMP_CFLAGS" = "x"],[
  AC_MSG_ERROR([OpenMP is required but the C compiler does not support it.
Check that CC (currently "$CC") wraps a compiler with OpenMP support
(e.g. mpicc wrapping GCC, Clang, or ICC).  Try:
  ./configure CC=mpicc CXX=mpicxx])
])

AC_LANG([C++])
AC_OPENMP
AS_IF([test "x$OPENMP_CXXFLAGS" = "x"],[
  AC_MSG_ERROR([OpenMP is required but the C++ compiler does not support it.
Check that CXX (currently "$CXX") wraps a compiler with OpenMP support
(e.g. mpicxx wrapping GCC, Clang, or ICC).  Try:
  ./configure CC=mpicc CXX=mpicxx])
])

AC_SUBST([OPENMP_CFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])

# ============================
# New Compiling Options (from README "New Compiling Options")
# https://github.com/Tomoaki-Ishiyama/mpi-rockstar/
# All default to "no" (disabled).
# ============================

# OUTPUT_RVMAX: output R_vmax (radius of maximum circular velocity)
# Note: always output for binary/HDF5; this option adds it to ASCII output.
AC_ARG_ENABLE([output-rvmax],
  AS_HELP_STRING([--enable-output-rvmax],
    [Output halo Rvmax (radius of max circular velocity) in ASCII output @<:@default=no@:>@]),
  [enable_output_rvmax=$enableval],[enable_output_rvmax=no])

# OUTPUT_NFW_CHI2: output chi^2 of NFW profile fitting
AC_ARG_ENABLE([output-nfw-chi2],
  AS_HELP_STRING([--enable-output-nfw-chi2],
    [Output chi^2 of NFW density profile fitting @<:@default=no@:>@]),
  [enable_output_nfw_chi2=$enableval],[enable_output_nfw_chi2=no])

# OUTPUT_INTERMEDIATE_AXIS: output 6 elements of intermediate shape ellipsoid axis
# NOTE: cannot be enabled together with OUTPUT_INERTIA_TENSOR
AC_ARG_ENABLE([output-intermediate-axis],
  AS_HELP_STRING([--enable-output-intermediate-axis],
    [Output 6 elements of intermediate shape ellipsoid axis (within Rvir and R500c) @<:@default=no@:>@]),
  [enable_output_intermediate_axis=$enableval],[enable_output_intermediate_axis=no])

# OUTPUT_INERTIA_TENSOR: output 12 elements of inertia tensor
# NOTE: cannot be enabled together with OUTPUT_INTERMEDIATE_AXIS
AC_ARG_ENABLE([output-inertia-tensor],
  AS_HELP_STRING([--enable-output-inertia-tensor],
    [Output 12 elements of inertia tensor (within Rvir and R500c) @<:@default=no@:>@]),
  [enable_output_inertia_tensor=$enableval],[enable_output_inertia_tensor=no])

# Mutual exclusion check: OUTPUT_INTERMEDIATE_AXIS and OUTPUT_INERTIA_TENSOR
# cannot both be enabled at the same time (README: "You can enable either ... at the same time")
AS_IF([test "x$enable_output_intermediate_axis" = "xyes" &&
       test "x$enable_output_inertia_tensor"    = "xyes"],[
  AC_MSG_ERROR([--enable-output-intermediate-axis and --enable-output-inertia-tensor
are mutually exclusive. Enable only one of them.
See: https://github.com/Tomoaki-Ishiyama/mpi-rockstar/ "New Compiling Options"])
])

ROCKSTAR_ADDFLAGS=""
AS_IF([test "x$enable_output_rvmax"             = "xyes"],[ROCKSTAR_ADDFLAGS="$ROCKSTAR_ADDFLAGS -DOUTPUT_RVMAX"])
AS_IF([test "x$enable_output_nfw_chi2"          = "xyes"],[ROCKSTAR_ADDFLAGS="$ROCKSTAR_ADDFLAGS -DOUTPUT_NFW_CHI2"])
AS_IF([test "x$enable_output_intermediate_axis" = "xyes"],[ROCKSTAR_ADDFLAGS="$ROCKSTAR_ADDFLAGS -DOUTPUT_INTERMEDIATE_AXIS"])
AS_IF([test "x$enable_output_inertia_tensor"    = "xyes"],[ROCKSTAR_ADDFLAGS="$ROCKSTAR_ADDFLAGS -DOUTPUT_INERTIA_TENSOR"])
AC_SUBST([ROCKSTAR_ADDFLAGS])

# Print summary of enabled output options
AC_MSG_NOTICE([--- New Compiling Options ---])
AC_MSG_NOTICE([  OUTPUT_RVMAX             : $enable_output_rvmax])
AC_MSG_NOTICE([  OUTPUT_NFW_CHI2          : $enable_output_nfw_chi2])
AC_MSG_NOTICE([  OUTPUT_INTERMEDIATE_AXIS : $enable_output_intermediate_axis])
AC_MSG_NOTICE([  OUTPUT_INERTIA_TENSOR    : $enable_output_inertia_tensor])
AC_MSG_NOTICE([-----------------------------])

# ============================
# HDF5 (optional)  --with-hdf5[=PREFIX]
# Builds both mpi-rockstar and mpi-rockstar_hdf5 when enabled.
# ============================
AC_ARG_WITH([hdf5],
  AS_HELP_STRING([--with-hdf5[=PREFIX]],
    [Enable HDF5 support. If PREFIX is omitted, use pkg-config to locate HDF5. Use --without-hdf5 to disable.]),
  [with_hdf5="$withval"],
  [with_hdf5="no"])

HDF5_CPPFLAGS=""
HDF5_LDFLAGS=""
HDF5_LIBS=""
HDF5_DEFS=""

AS_IF([test "x$with_hdf5" != "xno"],[
  # Match src/Makefile's HDF5_FLAGS: -DH5_USE_16_API -DENABLE_HDF5
  HDF5_DEFS="-DENABLE_HDF5 -DH5_USE_16_API"

  AS_IF([test "x$with_hdf5" = "xyes"],[
    PKG_CHECK_MODULES([HDF5],[hdf5],[
      HDF5_CPPFLAGS="$HDF5_CFLAGS"
      HDF5_LIBS="$HDF5_LIBS $HDF5_LIBS"
    ],[
      AC_MSG_ERROR([HDF5 requested via --with-hdf5, but pkg-config could not find package 'hdf5'.

Fix:
  - Install HDF5 development package (headers+libs) and its pkg-config file, or
  - specify a prefix:
      ./configure CC=mpicc CXX=mpicxx --with-hdf5=/path/to/hdf5])
    ])
  ],[
    hdf5_prefix="$with_hdf5"
    HDF5_CPPFLAGS="-I${hdf5_prefix}/include"
    HDF5_LDFLAGS="-L${hdf5_prefix}/lib -L${hdf5_prefix}/lib64"
    HDF5_LIBS="-lhdf5 -lz"
  ])

  # sanity checks
  save_CPPFLAGS="$CPPFLAGS"; save_LDFLAGS="$LDFLAGS"; save_LIBS="$LIBS"
  CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS $HDF5_DEFS"
  LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
  LIBS="$LIBS $HDF5_LIBS"
  AC_CHECK_HEADERS([hdf5.h],[],[
    AC_MSG_ERROR([HDF5 requested, but header 'hdf5.h' was not found. Check --with-hdf5=PREFIX.])
  ])
  AC_CHECK_LIB([hdf5],[H5Fopen],[],[
    AC_MSG_ERROR([HDF5 requested, but libhdf5 was not found or missing symbol H5Fopen. Check --with-hdf5=PREFIX.])
  ])
  CPPFLAGS="$save_CPPFLAGS"; LDFLAGS="$save_LDFLAGS"; LIBS="$save_LIBS"
])

AM_CONDITIONAL([WITH_HDF5],[test "x$with_hdf5" != "xno"])
AC_SUBST([HDF5_CPPFLAGS])
AC_SUBST([HDF5_LDFLAGS])
AC_SUBST([HDF5_LIBS])
AC_SUBST([HDF5_DEFS])

# ============================
# RPC / libtirpc (portable best practice):
#   - Do NOT add -I/usr/include/tirpc unless rpc/xdr.h is missing.
#   - Find library providing xdrmem_create (usually -ltirpc; sometimes none).
# ============================
TIRPC_CPPFLAGS=""
TIRPC_LIBS=""

AC_CHECK_HEADERS([rpc/xdr.h],[have_rpc_xdr=yes],[have_rpc_xdr=no])
AS_IF([test "x$have_rpc_xdr" = "xno"],[
  AC_CHECK_HEADERS([tirpc/rpc/xdr.h],[
    TIRPC_CPPFLAGS="-I/usr/include/tirpc"
  ],[
    AC_MSG_ERROR([Neither <rpc/xdr.h> nor <tirpc/rpc/xdr.h> was found.

Install RPC headers (often provided by libtirpc-dev/libtirpc-devel), or adjust CPPFLAGS.])
  ])
])

# Find the library that provides xdrmem_create
AC_SEARCH_LIBS([xdrmem_create],[tirpc nsl],[
  TIRPC_LIBS="$ac_cv_search_xdrmem_create"
],[
  AC_MSG_ERROR([Could not find xdrmem_create (needed for XDR/TIPSY I/O).
Install libtirpc and its development files.])
])
# Normalize "none required" to empty string
AS_IF([test "x$TIRPC_LIBS" = "xnone required"],[TIRPC_LIBS=""])

AC_SUBST([TIRPC_CPPFLAGS])
AC_SUBST([TIRPC_LIBS])

AC_CONFIG_FILES([
  Makefile
  src/Makefile
])
AC_OUTPUT
