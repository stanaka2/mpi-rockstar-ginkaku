# INSTALL

- Forked repository for GINKAKU
  - <https://github.com/stanaka2/mpi-rockstar-ginkaku>


## Features

- Support hierarchical IO for GINKAKU
- Support for lightcone extension for GINKAKU
- Support for the Hubble expansion laws of various cosmological models.
  - A table created using linear Boltzmann code is required.
- Support for the comoving distance of various cosmological models for lightcone mode.
  - A table created using linear Boltzmann code is required.
- Support for the virial over density of various $w_0$ and $w_a$
  - Because iteration computation takes much time in lightcone mode, support for reading the table created in `vir_dense_table.c`.
  - Even in normal mode, iteration on Fugaku takes time, so loading tables is advisable.
- Support switching between Gadget format $\left[\frac{1}{\sqrt{a}} \,\, \mathrm{km/s}\right]$ (default) or physical unit $\left[\mathrm{km/s}\right]$ for input particle velocity.
