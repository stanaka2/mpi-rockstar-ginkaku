# More Complete Setup #


## Input ##

### File Type ###

You will have to set the file type and conversion options. To do so, you will have to set:

    FILE_FORMAT = <your file format>

You can use `ASCII`, `GADGET`, `AREPO`, `ART`, or `TIPSY` to match your simulation file.
If you use the `ART` option,
only PMss files are currently supported.  (PMcrs files, which do not include particle IDs
explicitly, are not supported).
The type `KYF` is used for the test data provided by us.
You can also use `GADGET4`, which is not supported by the original rockstar. 
See [Gadget-2/3/4 and AREPO HDF5 Format](gadget4_support). 
    
If you use `ASCII`, `ART`, or `TIPSY`, you will additionally have to set
the particle mass; if you use `ASCII` or `TIPSY`, you will also have to
set parameters for the cosmology and box size.
    
If you use `GADGET`, you should check to make sure that the length
and mass conversion multipliers are correct (`GADGET_LENGTH_CONVERSION` and
`GADGET_MASS_CONVERSION`) to convert Gadget internal units to comoving Mpc/h
and Msun/h.
    
Note that `GADGET` only works with the binary GADGET formats.  For GADGET
as well as AREPO HDF5 snapshots, use `AREPO` for the file format.  You will then have to
set `AREPO_LENGTH_CONVERSION` and
`AREPO_MASS_CONVERSION` to convert internal snapshot units to comoving Mpc/h
and Msun/h, respectively.  Make sure to compile Rockstar with "`make with_hdf5`";
otherwise, the code will not accept these config options.
    
If you use `TIPSY`, you should set the length and velocity
conversion multipliers (`TIPSY_LENGTH_CONVERSION` and `TIPSY_MASS_CONVERSION`)
to convert Tipsy internal units to comoving Mpc/h and physical km/s.
Note that Tipsy format is in a _beta stage of support_---please contact
me if you have issues!


(gadget4_support)=
### Gadget-2/3/4 and AREPO HDF5 Format ###

Gadget-4 (HDF5) format is newly supported as input snapshots.
```
FILE_FORMAT="GADGET4"
```

The default number of Gadget-4 particle types is `NTYPES=6`. You can use other numbers by adding the line below to a configuration file.
```
GADGET4_NTYPES=<your ntypes>
```
This is the same as AREPO HDF5 format. You can pass the number for AREPO format in a similar manner.
```
AREPO_NTYPES=<your ntypes>
```

In Gadget4 `Config.sh`, you can set `IDS_64BIT` for 8-byte IDs or `IDS_32BIT` for 4-byte IDs. Both are supported but 6-byte ID (`IDS_48BIT`) is not accepted.
`MPI-Rockstar` can load snapshots with single- and double-precision particle positions and velocities.
However, `MPI-Rockstar` always internally stores positions and velocities in single precision. The output catalogue is also in single precision.

Gadget-2/3 HDF5 format is also supported by passing the file format as `AREPO`.
This format has a different header structure from Gadget-4 but AREPO HDF5 format has the same structure.
Note that when you use Gadget-2/3 HDF5 format with `FILE_FORMAT = AREPO`, set `AREPO_MASS_CONVERSION`, `AREPO_LENGTH_CONVERSION`, etc. if needed,
instead of `GADGET_MASS_CONVERSION`, `GADGET_LENGTH_CONVERSION`, etc., which are effective only for binary Gadget formats.
Be careful that `AREPO_LENGTH_CONVERSION = 1e-3` by default, i.e., the length unit is kpc/h.



### File Name, Number of Files ###
 
Besides setting the file type and conversion options in the previous section, 
several more options are necessary in the configuration file to specify
filenames, paths.

You will also have to provide information on where the data files are located and the number of files per snapshot.  To do so, you will have to set:

    NUM_BLOCKS = <number of files per snapshot>    
    INBASE = "/directory/where/files/are/located"
    FILENAME = "my_sim.<snap>.<block>"

In this example, "`my_sim`" should be the base of your simulation filename.
For multiple snapshots, the text "`<snap>`" will be automatically replaced
by the snapshot number, and the text "`<block>`" will be automatically
replaced by the block number (0 to `NUM_BLOCKS`-1).  To specify the range
of snapshot numbers, you may set
    
    NUM_SNAPS = <total number of snapshots>
    STARTING_SNAP = <first snap> #defaults to 0

If you have nonstandard names for your snapshots (e.g., "001" instead of
"1"), (Since GADGET always uses "001" instead of "1" for its snapshot numbers, Rockstar will automatically use the correct formatting without needing a snapshot filename list if the particle file type is `GADGET` or `AREPO`.) you may create a text file with one snapshot name per line and set:
    
    SNAPSHOT_NAMES = </path/to/snapshot names>

This will automatically override the number of snapshots _and_ the starting
snapshot, if specified.  You may do the same thing for the block names,
too:
    
    BLOCK_NAMES = </path/to/block names>

(Note that `NUM_BLOCKS` should still be specified in this case, so that the
server knows how many reader connections to accept).

To restart Rockstar from a specific snapshot, use
    
    RESTART_SNAP = <snapshot number> #default: 0

This will resume analysis with the assumption that all snapshots previous to `RESTART_SNAPSHOT` have been analyzed.


(input_subd)=
### Sub-Directory for Input (if necessary) ###

The options, `FILES_PER_SUBDIR_INPUT` and `SUBDIR_DIGITS_INPUT`, allow that snapshots are stored in multiple sub-directories. `FILES_PER_SUBDIR_INPUT` is the number of files in each sub-directory. The sub-directory name is expressed by a sequential number starting from zero, and its digit is set by `SUBDIR_DIGITS_INPUT` For example, if
```
INBASE="snapdir"
FILENAME="snap-<block>"
NUM_BLOCKS=16
FILES_PER_SUBDIR_INPUT=4
SUBDIR_DIGITS_INPUT=4
```
in a configuration file, the following directory structure of snapshots is assumed.
```
snapdir/
  0000/snap-[0-3]
  0001/snap-[4-7]
  0002/snap-[8-11]
  0003/snap-[12-15]
```
When `FILES_PER_SUBDIR_INPUT=0`, no sub-directories are assumed (by default).



(output)=
## Output ##

MPI-Rockstar can output in many different formats. 
Halos are saved to files called `out_0.list`, `out_1.list`, etc., with one file per particle snapshot. 
When you want to write `out_<snap>.list` per process, set `OUTLIST_PARALLEL=1` (0 by default).


(halo_properties)=
### Halo Properties ###

Many halo properties are currently calculated by default:

+  Halo masses at several radii: Mvir, M200b, M200c, M500c, M2500c.  These masses always include any contributions from substructure.  Also, masses with higher density thresholds (e.g., 2500c) can sometimes be zero if the density of the halo never rises above the threshold. [Footnote: Particles are assumed to have an effective radius of `FORCE_RES` for the purposes of density calculations near the halo center.}  By default, only bound particles are included; see [Commonly-Used Options](commonly-used-options) if this is not what you want.
+  Halo maximum circular velocity and velocity dispersion.
+  Halo radii: Rvir and the scale radius r_s, calculated both using profile fitting and using the Klypin v_max method.
+  Halo center positions and velocities.
+  Halo spin (both Bullock and Peebles) and angular momentum.
+  Halo shapes and principal axes, using the Allgood method (iterative, weighted by 1/r^2), at both Rvir and R500c.  See [Commonly-Used Options](commonly-used-options) for how to change the shape calculation method.
+  The ratio of halo kinetic to potential energy, and the center position and velocity offsets from the halo's bulk average position and velocity.

Information about each of the columns in both file types (including units) is available in the ASCII headers.  If you would like to calculate more halo properties, please see [Adding More Halo Properties](adding_more_halo_properties).

MPI-Rockstar can output additional halo properties by adding the below options at compiling.

- **OUTPUT_RVMAX**:
Output halo's $R_{\rm vmax}$, which is a radius of maximum circular velocity. This option is only applicable to ascii format, and is always output by default for binary and HDF5 format.
- **OUTPUT_NFW_CHI2**:
Output $\chi^2$ in the NFW fitting of the density profile of a halo. The NFW fitting is used to calculate the scale radius $r_{\rm s}$.
When the number of halo's bound particles is less than `MIN_SCALE_PART` (set in `nfw.c`, 100 by default), $r_{\rm s}$ is the same with $r_{\rm s}^{\rm klypin}$, which is the estimated scale radius using $V_{\rm max}$ and the halo mass. In this case, $\chi^2=0$ is set.
When this option is enabled, the binary file header will have the corresponding bit in `add_flag` set using `add_flag |= 1;`.
- **OUTPUT_INTERMEDIATE_AXIS**:
Output six elements of halo's intermediate shape ellipsoid axis. Three and remaining three of them are calculated by using bound particles within $R_{\rm vir}$ and $R_{\rm 500c}$, respectively. When this option is enabled, the binary file header will have the corresponding bit in `add_flag` set using `add_flag |= 2;`.
- **OUTPUT_INERTIA_TENSOR**:
Output 12 elements of halo's inertia tensor. Six and remaining six of them are calculated by using bound particles within $R_{\rm vir}$ and $R_{\rm 500c}$, respectively. You can enable either `OUTPUT_INERTIA_TENSOR` and `OUTPUT_INTERMEDIATE_AXIS` at the same time.
When this option is enabled, the binary file header will have the corresponding bit in `add_flag` set using `add_flag |= 4;`.



(file_format)=
### File Format ###

Beside each process can also output halo catalogs in either ASCII or binary formats 
(`halos_<snap>.<process id>.ascii` or `halos_<snap>.<process id>.bin` files).
By default, both binary and ASCII catalogs are printed.  This is set by
    
    OUTPUT_FORMAT = BOTH # or "ASCII" or "BINARY"

However, as the binary outputs are required for generating merger trees,
merger trees will _not_ be generated if you select the `ASCII` option.
If you want merger trees but not the binary outputs (which take up
lots of space), you can set the following option:
    
    DELETE_BINARY_OUTPUT_AFTER_FINISHED = 1

For the binary outputs, there is a 256-byte header (detailed in
`io/io_internal.h`) followed by a binary dump of the halo structures
(see `halo.h`), followed by a binary dump of the particle ids in
each halo (type `int64_t`).  See the `load_binary_halos()` routine in
`io/io_internal.c` for an example of how to read in the binary files.

MPI-Rockstar supports the HDF5 output of halo catalogs by adding this line in a configuration file.
```
OUTPUT_FORMAT="HDF5"
```


(output_subd)=
### Sub-Directory for Output (if necessary) ###

By defaults, MPI-Rockstar writes all output (bin, ascii, list) in single directory (set by `OUTBASE` in a configuration file). When you want to separate output directories for every snapshot, set `OUTPUT_SUBDIR=1`. The sub-directory name is expressed by the snapshot number with the digits set by `SNAPSHOT_SUBDIR_DIGITS`. For example, if
```
STARTING_SNAP=20
NUM_SNAPS=23
OUTBASE="Out"
OUTPUT_SUBDIR=1
SNAPSHOT_SUBDIR_DIGITS=3
```
in a configuration file, the output are stored in the following directory.
```
Out/
  020/ for snapshot 20
  021/ for snapshot 21
  022/ for snapshot 22
```
**Note that each sub-directory must be created before the running.**

The options, `FILES_PER_SUBDIR_OUTPUT` and `SUBDIR_DIGITS_OUTPUT`,  are inverse of `FILES_PER_SUBDIR_INPUT` and `SUBDIR_DIGITS_INPUT`. These options allow that output of a snapshot is stored in multiple sub-directories. `FILES_PER_SUBDIR_OUTPUT` is the number of files in each sub-directory. The sub-directory name is expressed by a sequential number starting from zero, and its digit is set by `SUBDIR_DIGITS_OUTPUT` For example, if the number of MPI processes is 16 and
```
OUTBASE="Out"
FILES_PER_SUBDIR_OUTPUT=4
SUBDIR_DIGITS_OUTPUT=4
```
in a configuration file, the output are stored in the following directory.
```
Out/
  0000/ for processes [0-3]
  0001/ for processes [4-7]
  0002/ for processes [8-11]
  0003/ for processes [12-15]
```
**Note that each sub-directory is made automatically.**

**`OUTLIST_PARALLEL`, `OUTPUT_SUBDIR`, `SNAPSHOT_SUBDIR_DIGITS`, `FILES_PER_SUBDIR_OUTPUT` and `SUBDIR_DIGITS_OUTPUT` can work concurrently. In this case, sub-sub-directories are made.**


### Output Threshold ###

To change the minimum particle size of output halos, set
    
    MIN_HALO_OUTPUT_SIZE = <minimum number of particles> #default: 20

Note that this option does _not_ specify the minimum number of particles
bound within the virial or other halo radius.  Instead, it specifies the
minimum number of particles which Rockstar identifies as uniquely assigned
to that halo.  The virial mass can (and usually will) be much less.

To avoid printing spurious halos, Rockstar excludes significantly 
unbound objects by default.  To change the threshold of what is considered
"significantly unbound"---e.g., for studying tidal remnants, you may
consider lowering the value of the following parameter:      
    
    UNBOUND_THRESHOLD = <minimum fraction of bound mass> #default: 0.5

Generally, this affects the halo mass function at the few percent level
for halos of 100 bound particles or less (see the Rockstar science paper for details).


### Output Full Particle Information ### 

Rockstar does not by default output full particle information to
save space.  In order to turn on this capability, set
    
    FULL_PARTICLE_CHUNKS = <number of processes which will output particles>

to the desired number of processes which you want to output
particles (i.e., something between 1 and `the number of MPI processes`) at every timestep.


(merger_tree)=
## Merger Trees ##

By default, the output files for Rockstar include some basic descendant
information for each halo.  However, it is strongly recommended that
you use a more sophisticated package, such as Consistent Trees

<http://code.google.com/p/consistent-trees>
 
to improve consistency
between timesteps before using the descendant information for science
purposes (e.g., calculating merger rates).  Rockstar includes a script
to autogenerate configuration files for Consistent Trees (see the
Consistent Trees README file for use).  


## Host / Subhalo Relationships ##

By default, the output files for Rockstar do not include information
about which halos are hosts (as opposed to subhalos).  This is because
the consistent tree code (see [Merger Trees](merger_tree), above) automatically
calculates such relationships.  However, if you do not intend to use the
merger tree code, then there is a utility included with Rockstar to
postprocess the output halo catalogs (i.e., the `out_*.list` files). 
You can compile and run the following command

    make find_parents -C src    
    ./find_parents out_XYZ.list <box_size>

on each halo catalog for which you want to find host halos.  Information
about which halos are hosts and which are subs will be output as an
extra column, the parent ID `PID`.  Host halos will have a parent ID of -1;
subhalos will have a parent ID which corresponds to their host halo ID.
Note that the files produced will _not_ be suitable for input into the
merger tree code above.

## Lightcones ##

Lightcones are currently only supported with parallel IO (see options
above).  To turn on lightcone support, set
    
    LIGHTCONE = 1

The scale factor will be automatically computed
from the distance to the origin, which is set by
    
    LIGHTCONE_ORIGIN = (x,y,z) #default: (0,0,0)

If you want a lightcone created by joining two different fields of view
from the same box, you will have to set the following options:
    
    LIGHTCONE_ALT_SNAPS = /path/to/second_lightcone_snapnames
    LIGHTCONE_ALT_ORIGIN = (x,y,z) #default: (0,0,0)

(The `LIGHTCONE_ALT_ORIGIN` option translates all particles in the second
field of view from (x,y,z) to the first `LIGHTCONE_ORIGIN`.)  If using two
fields of view, it is assumed that they have equal numbers of input files;
`NUM_BLOCKS` should be set to the combined total number of input files.

If you have more than two fields of view from the same box, you will
have to make sure that their positions are already translated
appropriately so that they can be joined/unioned together without any
overlap issues.  In this case, you should not set any of the ALT lightcone
options.  Additionally, you should set "`IGNORE_PARTICLE_IDS = 1`" so as
to prevent confusion with the same particle ID appearing multiple times.


