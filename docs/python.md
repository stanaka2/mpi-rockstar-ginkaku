# Building as a Library and Python wrapper

## Building as a Library ##

You can build a static library for linking with other MPI applications:
```
make libmpi_rockstar -C src
```
The host application must initialize MPI and load a configuration before calling `mpi_main`.  An example program is available in `examples/library_example.c`:
```
#include <mpi.h>
#include "config.h"
#include "mpi_rockstar.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    do_config("parallel_256.cfg");
    mpi_main(0, NULL);
    MPI_Finalize();
    return 0;
}
```

In some environments, due to a compatibility issue, you may encounter this or something similar error.
```
/usr/include/tirpc/rpc/xdr.h:111:52: error: unknown type name 'u_int'
```
You may be able to solve it by removing `-I/usr/include/tirpc` from `.c.o:` in the `Makefile`.


## Python wrapper pympirockstar

Python bindings for the MPI-Rockstar halo finder.

### Prerequisites

* A compiled MPI-Rockstar library in `../src`.
  Build the library from the repository root with:

  ```bash
  make libmpi_rockstar -C src               # standard build
  # or
  make libmpi_rockstar_hdf5 -C src          # with HDF5 support
  ```

### Installation

Install the wrapper with `pip`, using your MPI C/C++ compiler.  From the
repository root:

```bash
MPICC=mpicxx pip install -e ./pympirockstar
```

If you built the HDF5 variant, set `ROCKSTAR_LIB=mpi_rockstar_hdf5` before
running `pip` so the extension links against the correct library.
