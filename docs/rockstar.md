# Options from the Original Rockstar

(commonly-used-options)=
## Commonly-Used Options ##

For those options not mentioned directly above:
    
    MASS_DEFINITION = <"vir" or "XXXc" or "XXXb" etc.> #default: vir

This lets you specify how you want masses calculated.  "`vir`" uses the
formula from Bryan \& Norman (1998); a number plus "c" or "b" calculates
masses relative to the critical or background density, respectively.
(E.g. "`200b`" or "`200c`", or even fractional values like "`100.5c`").
Changing the main mass definition will also change the halo radius at which most other
properties are calculated (e.g., v_max, spin, shape, etc.).

Besides the main mass definition, you can specify four other halo mass definitions to calculate,
in `MASS_DEFINITION2` through `MASS_DEFINITION5`.  The syntax for these is
identical to that for `MASS_DEFINITION`; the defaults are listed in Section  [Output](output).  Note that
for very low density thresholds (e.g., <150b), you may have to increase the value of `FOF_LINKING_LENGTH` below to obtain accurate masses.

<!-- 
By default, all masses are calculated from FOF groups after
unbound particles are removed.  This is appropriate for merger trees,
where halos/subhalos need to be treated identically for consistency
across timesteps.  However, if you are calculating _mass functions_,
you will want to use spherical overdensities including unbound particles
and particles which may exist outside of the FOF group for the halo.
To enable these strict SO masses, use
    
    STRICT_SO_MASSES = 1 #default: 0

The resulting masses will be present in the auxilliary mass columns
(columns 20-24) of the `out_*.list` files.  
-->

If computing mass functions for snapshots widely separated in redshift, you may also wish to disable temporal halo finding (i.e., using information from previous snapshots to determine host/sub relationships): 
    
    TEMPORAL_HALO_FINDING = 0 #default: 1

Shape calculations are by default performed with the Allgood method.  If you wish to change this, you can specify
    
    WEIGHTED_SHAPES = 0 #default: 1

which will calculate the mass tensor without weighting by 1/r^2.  In addition, if you do not want the shapes to be calculated iteratively, specify
    
    SHAPE_ITERATIONS = 1 #default: 10

This setting specifies the maximum number of shape calculation iterations to perform.

Options to set the base cosmology are available:
    
    SCALE_NOW = <current cosmological scale factor>
    h0 = <hubble constant today> # in units of 100 km/s/Mpc
    Ol = <Omega_Lambda> # in units of the critical density
    Om = <Omega_Matter> # in units of the critical density

These cosmology options are only relevant if one is reading from an ASCII
or TIPSY particle file.  (For TIPSY, the `SCALE_NOW` parameter may be
omitted).  For other data formats, these values are \texttt{automatically
overwritten} with values in the particle data files.

Non-LCDM equations of state are also available.  Since none of the particle input formats include a standard way of specifying these, the equation of state must be specified in the config file:
    
    W0 = <W_0> # dark energy equation of state at z=0
    WA = <W_A> # scaling of DE equation of state with scale factor.

By default, these are set to `W0=-1` and `WA=0`.

It is not necessary to specify the particle mass for GADGET2 files, but it is necessary for all other file types:
    
    PARTICLE_MASS = <mass of each particle> #in Msun/h


However, for GADGET2 files, you will need to specify the conversion to
Rockstar's internal length (comoving Mpc/h) and mass (Msun/h) units:
    
    GADGET_MASS_CONVERSION = <conversion from GADGET units to Msun/h>
    GADGET_LENGTH_CONVERSION = <conversion from GADGET units to Mpc/h>

Usually, these will be `1e10` (mass conversion) and either `1` or
`1e-3` (length conversion for Mpc/h and kpc/h, respectively).

Similarly, for TIPSY files, you will need to convert lengths and
velocities to comoving Mpc/h and physical km/s:
    
    TIPSY_LENGTH_CONVERSION = <conversion from TIPSY units to Mpc/h>
    TIPSY_VELOCITY_CONVERSION = <conversion from TIPSY units to km/s>

To disable periodic boundary conditions (only applicable for `PARALLEL_IO`),
you can set:
    
    PERIODIC = 0

For ASCII particle data, it is necessary to specify the box size:
    
    BOX_SIZE = <side length of cosmological box in comoving Mpc/h>


## Rarely-used Options ##
   
    GADGET_SKIP_NON_HALO_PARTICLES = <0 or 1> #default = 1

By default, Rockstar only considers dark matter particles; the preceding
option can be set to 0 to force consideration of other particles as well
in GADGET2 files.  The default halo particle type in GADGET is 1; however,
if you need to change this, you can use the following option:
    
    GADGET_HALO_PARTICLE_TYPE = <0 to 5> #default = 1

Note that Rockstar has no current support for multiple particle masses.  A beta version of Rockstar with support for multi-mass particles is available on request from the authors.
    
    RESCALE_PARTICLE_MASSES = <0 or 1> #default 1

If only dark matter particles are used from GADGET2 files in a simulation
which also includes gas particles, it is necessary to rescale the particle
masses so as to preserve the correct matter density; setting this option
tells Rockstar to do so.

If for some reason your simulation data has inconsistent or duplicate
particle IDs, you can set the following option to prevent problems with
halo finding:
    
    IGNORE_PARTICLE_IDS = 1

Note that in this case, merger trees are disabled.


The following options determine details of the halo finding.  To tell
Rockstar to calculate halo radii as well as Vmax, Rvmax, and other
halo properties from the *unbound* particles, set the following option
to 0:
    
    BOUND_PROPS = <0 or 1> #default 1

This is likely not what you want if you are calculating mass functions; you should use `STRICT_SO_MASSES` instead ([Commonly-Used Options](commonly-used-options)).

To change the default FOF refinement fraction (see the Rockstar paper),
change:
    
    FOF_FRACTION = <FOF refinement fraction> #default 0.7

To change the default 3D FOF linking length, set:
    
    FOF_LINKING_LENGTH = <FOF linking length> #default 0.28

This is generally only necessary to increase if trying to find halos with very low spherical overdensity thresholds.

To change the minimum (internal) number of particles considered to be
a halo seed, change:
    
    MIN_HALO_PARTICLES = <minimum number of halo particles> #default: 10

See also `MIN_HALO_OUTPUT_SIZE` in Section  [Output](output).


## Extending Rockstar ##

### Adding More Configuration Parameters ###

Config parameters are defined in `config.template.h`.  The syntax is

       <type>(VARIABLE_NAME, DEFAULT_VALUE)

where `<type>` is one of
"`string`" (`char *`), "`integer`" (64-bit signed integer, i.e., `int64_t`),
"`real`" (`double`), or "`real3`" (`double[3]`).

Config parameters added here will be globally visible in the code
and automatically interpreted from the config file.  The default
value will be assigned if the config file does not include the
variable.


(adding_more_input_formats)=
### Adding More Input Formats ###

For highest efficiency, it is appropriate to copy one of the existing `io/io_*.c` routines and modify it to accept your input format directly.  You will have to define a routine that looks like:

    void load_particles_mytype(char *filename, struct particle **p, int64_t *num_p)

This routine should add the number of particles in the file to `num_p` and grow the particle structure appropriately, e.g.,

     check_realloc_s(p[0], sizeof(struct particle), num_p[0] + num_new_p);

Note that `check_realloc_s()` is defined in `check_syscalls.h`.  It is generally a good idea to use the IO routines defined in that header, as they will automatically halt the program if they encounter any errors reading, writing, or opening files, making debugging much easier.  Routines helpful for reading Fortran-unformatted files and files encoded in a different endianness are available from `io/io_util.h`.

The particle structure is very simple; it is described in `particle.h`.  Your routine should convert everything to Rockstar's internal units: i.e., comoving Mpc/h for positions and non-comoving km/s for peculiar velocities.  Your routine should also set as many config variables as possible (e.g., `BOX_SIZE`, `PARTICLE_MASS`, cosmology, current scale factor, etc.) from the data header, as that will reduce the chance for errors in the analysis from incorrectly-supplied config variables.

Once you have created a new `io/io_mytype.c` file, you should create a corresponding header file and include it in `io/meta_io.c`.  You can then add a new condition in the `read_particles()` routine in `io/meta_io.c` to call your code when the `FILE_FORMAT` variable matches your filetype.  Finally, you should add the new `io/io_mytype.c` file to the `Makefile` so that it gets compiled along with the other Rockstar source code.

If you think your code may be useful for others, please consider submitting it as a patch to the authors.  Others will no doubt appreciate it!


(adding_more_output_formats)=
### Adding More Output Formats ###

Halos are output in the `output_halos()` routine in `io/meta_io.c`.  Halos are stored in the `halos` array, for which the structure is defined in `halo.h`.  Particles are stored in the `p` array, for which the structure is defined in `particle.h`.  It is best to use the `output_binary()` routine from `io/io_internal.c` as a starting template for outputting halo and particle data.


(adding_more_halo_properties)=
### Adding More Halo Properties ###

All halo properties are calculated in `properties.c`, in the routine `_calc_additional_halo_props()`.  All particles associated with the current halo being analyzed are in the `po` array, which is of type `struct potential` (defined in `potential.h`); particles in the array are sorted by distance from the halo center.  You may either add calculations directly in the main particle loop in `_calc_additional_halo_props()`, or you may decide to create your own particle loop as a separate routine.  If the latter case, don't forget to honor the option to skip over unbound particles unless your calculation always needs to include them; otherwise, you may get strange results for subhalos.

You may store your calculated values by extending the halo structure, defined in `halo.h`.  To print out the new properties, you should modify `gen_merger_catalog()` in `io/meta_io.c`.  As with new input formats, if you think your code may be useful for others, please consider submitting it as a patch to the authors.

