/*
mpicc -std=gnu11 -O3 ./read_vir_dens_table.c -lm -lgsl -lgslcblas
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include "config_vars.h"
#include "hubble.h"
#include "vir_dense.h"

const double  table_max_z = 50.0; // It does not work well over z>70.
static double vir_dens_table[2][2048];
extern double hubble_table[2][2048];
extern int    hubble_table_size;

void create_vir_density_table(int mpi_nproc, int mpi_rank) {
    int ntable = hubble_table_size;
    assert(ntable > 100);

#pragma omp parallel for schedule(auto)
    for (int ia = 0; ia < 2048; ia++) {
        vir_dens_table[0][ia] = 0.0;
        vir_dens_table[1][ia] = 0.0;
    }

    double *a_table;
    a_table = hubble_table[0]; // ascale

    int _mpi_nproc = mpi_nproc;

    if (_mpi_nproc > 2048)
        _mpi_nproc = 2048;

    if (mpi_rank < _mpi_nproc) {

#pragma omp parallel for schedule(auto)
        for (int ia = 0; ia < ntable; ia++) {
            double a = a_table[ia];
            a        = pow(10.0, a); // log to linear

            vir_dens_table[0][ia] = a;
            if (ia % _mpi_nproc != mpi_rank)
                continue;

            if (a < 1.0 / (1.0 + table_max_z)) {
                vir_dens_table[1][ia] = 0.0;
                continue;
            }

            double x, y, y_appr, a_ta, zeta, eta_ta, eta_vir;
            double delta_lin_col = 0.;

            double delta0    = 1.;
            double threshold = 1e7;
            double deltanl;

            double a_col   = a;
            double dgrowth = get_linear_growth(a_col) / get_linear_growth(1.);

            double delta0_min = 1. / dgrowth;
            double delta0_max = 2. / dgrowth;

            int    err_flag    = 0;
            double prev_delta0 = 1.0, prev_deltanl = 1.0;

            while (1) {
                delta0  = (delta0_min + delta0_max) / 2.;
                deltanl = spherical_collapse(a_col, delta0);

                // printf("%d %g, %g %g %g %g %g %g\n", ia, a, delta0, deltanl,
                // delta0_min, delta0_max, delta_lin_col, threshold);
                // fflush(stdout);

                if (deltanl > threshold) {
                    delta_lin_col = delta0;
                    break;
                } else if (isnan(deltanl)) {
                    delta0_max = delta0;
                } else {
                    delta0_min = delta0;
                }

                if ((fabs((delta0 - prev_delta0) / prev_delta0) < 1e-4) &&
                    (fabs((deltanl - prev_deltanl) / prev_deltanl) < 1e-4)) {
                    err_flag = 1;
                    break;
                }

                prev_delta0  = delta0;
                prev_deltanl = deltanl;
            }

            if (err_flag == 1) {
                vir_dens_table[1][ia] = 0.0;
                continue;
            }

            calc_a_ta_zeta(a_col, delta_lin_col, &a_ta, &zeta);
            x      = a_col / a_ta;
            eta_ta = 2. * (Omegade_a(a_ta) / Omegam_a(a_ta)) / zeta;
            eta_vir =
                2. * (Omegade_a(a_col) / Omegam_a(a_col)) * pow(x, -3.) / zeta;

            y_appr = (1. - eta_vir / 2.) / (2. + eta_ta - 3. * eta_ta / 2.);
            y      = Newton_method(eta_vir, eta_ta);

            double Delta_vir = zeta * pow(x / y, 3.);

            vir_dens_table[1][ia] = Delta_vir;
            printf("%d %d %g %g %g\n", mpi_rank, ia, a, (1.0 / a - 1.0),
                   Delta_vir);
            fflush(stdout);
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, vir_dens_table[1], ntable, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
}

int main(int argc, char **argv) {

    assert(argc == 3);

    MPI_Init(&argc, &argv);

    int mpi_nproc, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &(mpi_nproc));
    MPI_Comm_rank(MPI_COMM_WORLD, &(mpi_rank));

    char *input_hubble_file = argv[1];
    char *output_vd_file    = argv[2];

    EXPANSION_TABLEFILE = input_hubble_file;

    initialize_linear_growth();

    if (mpi_rank == 0) {
        printf("initialized linear growth\n");
        fflush(stdout);
    }

    create_vir_density_table(mpi_nproc, mpi_rank);

    if (mpi_rank == 0) {
        printf("create virial density table\n");
        fflush(stdout);
    }

    if (mpi_rank == 0) {
        FILE *fp;
        fp = fopen(output_vd_file, "w");
        fprintf(fp, "# a Delta_vir\n");

        int ntable = hubble_table_size;
        for (int ia = 0; ia < ntable; ia++) {
            if (vir_dens_table[1][ia] < 1.0e-3)
                continue;
            fprintf(fp, "%.10e %.10e\n", vir_dens_table[0][ia],
                    vir_dens_table[1][ia]);
        }

        fflush(fp);
        fclose(fp);
    }

    MPI_Finalize();

    exit(0);
}
