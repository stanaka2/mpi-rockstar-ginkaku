#include <math.h>
#include <assert.h>
#include "config_vars.h"
#include "hubble.h"

#define MAX_Z      300.0
#define Z_BINS     1000.0
#define TOTAL_BINS (((int)MAX_Z) * ((int)Z_BINS))
#define h          h0
double Dh;
double _Dc[TOTAL_BINS];

#ifndef __APPLE__
#ifndef isfinite
#define isfinite(x) finitef(x)
#endif
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288 /* From math.h */
#endif                                             /* M_PI */

#define _E(z) (hubble_scaling(z))

double redshift(double a) { return (1.0 / a - 1.0); }
double scale_factor(double z) { return (1.0 / (1.0 + z)); }

void init_cosmology(void) {
    int    i;
    double z;
    double Dc_int = 0;
    Dh            = 2997.92458 / h; // In Mpc  (Speed of light / H0)
    for (i = 0; i < TOTAL_BINS; i++) {
        z      = (float)(i + 0.5) / Z_BINS;
        _Dc[i] = Dc_int * Dh;
        Dc_int += 1.0 / (_E(z) * Z_BINS);
    }
}

// switch for GINKAKU
// add by stanaka
#ifndef FOR_GINKAKU

double comoving_distance(double z) {
    double f   = z * Z_BINS;
    int    bin = (int)f;
    if (z < 0)
        return 0;
    if (bin > (TOTAL_BINS - 2))
        return (_Dc[TOTAL_BINS - 1]);
    f -= bin;
    return (_Dc[bin] * (1.0 - f) + _Dc[bin + 1] * f);
}

double comoving_distance_h(double z) { return (comoving_distance(z) * h); }

#else
// for GINKAKU
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int    set_cdist_table(char *, double);
double read_cdist_table_log(const int, const double, const int);

double cdist_table[2][2048];
int    cdist_table_size;

// original comoving distance(z) [Mpc]
double _comoving_distance(double z) {
    double f   = z * Z_BINS;
    int    bin = (int)f;
    if (z < 0)
        return 0;
    if (bin > (TOTAL_BINS - 2))
        return (_Dc[TOTAL_BINS - 1]);
    f -= bin;
    return (_Dc[bin] * (1.0 - f) + _Dc[bin + 1] * f);
}

// original comoving distance(z) [Mpc/h]
double _comoving_distance_h(double z) { return (_comoving_distance(z) * h); }

/* [Mpc/h] */
double comoving_distance_h(double z) {
    double z1 = 1.0 + z;
    double a  = 1.0 / z1;

    static int set_table = 0;
    static int ltable    = 0;

    if (strcmp(CDIST_TABLEFILE, "none") == 0) {
        return _comoving_distance_h(z);
    }

    // const double amin = 1.0e-2;  // z ~ 100
    const double amin =
        1.0e-6; // smaller than minimum ascale of initialize_linear_growth()

    if (set_table == 0) {
        ltable           = set_cdist_table(CDIST_TABLEFILE, amin);
        cdist_table_size = ltable;
        printf("read comoving distance[Mpc/h] table %s. table length=%d\n",
               CDIST_TABLEFILE, ltable);
        fflush(stdout);
        set_table = 1;
    }

    // exception handling for first step
    if (a < amin || a > 0.999999)
        return 0.0;

    const int order = READ_TABLE_ORDER;
    return read_cdist_table_log(order, a, ltable);
}

double comoving_distance(double z) { return (comoving_distance_h(z) / h); }

int set_cdist_table(char *inputfile, const double amin) {
    FILE *fp;
    fp = fopen(inputfile, "r");

    if (fp == NULL) {
        printf("Error: %s file not opened.", inputfile);
        return EXIT_FAILURE;
    }

    double atbl0 = 0.9 * amin;
    double a, cdist;

    int line = 0;

    static char tmp[128];
    size_t      len = 128;

    // skip header
    fgets(tmp, len, fp);

    while (fscanf(fp, "%lf %lf", &a, &cdist) != EOF) {
        if (atbl0 > a)
            continue;
        line++;
    }

    // global table size
    assert(line < 2048);

    line = 0;
    fseek(fp, 0, SEEK_SET);

    // skip header
    fgets(tmp, len, fp);

    while (fscanf(fp, "%lf %lf", &a, &cdist) != EOF) {
        if (atbl0 > a)
            continue;
        // hubble_table[0][line] = a;
        cdist_table[0][line] = log10(a);
        cdist_table[1][line] = cdist;
        line++;
    }

    fclose(fp);

    return line;
}

double read_cdist_table_log(const int order, const double anow,
                            const int ntable) {
    double *log_a_table, *cd_table;
    log_a_table = cdist_table[0]; // log10(ascale)
    cd_table    = cdist_table[1]; // Mpc/h

    const double log_anow = log10(anow);

    double log_a_min = log_a_table[0];
    // double log_a_max = log_a_table[ntable - 1];
    double da    = log_a_table[1] - log_a_table[0];
    int    a_idx = (log_anow - log_a_min) / da;

    static int lsearch  = 3;
    int        hit_flag = 0;

    for (int idx = a_idx - lsearch; idx < a_idx + lsearch; idx++) {
        if (log_a_table[idx] <= log_anow && log_anow < log_a_table[idx + 1]) {
            a_idx    = idx;
            hit_flag = 1;
            break;
        }
    }

    if (hit_flag == 0) {
        for (int idx = 0; idx < ntable; idx++) {
            if (log_a_table[idx] <= log_anow &&
                log_anow < log_a_table[idx + 1]) {
                a_idx = idx;
                break;
            }
        }
        lsearch++;
    }

    double cdist    = 0.0;
    int    a_idx_p1 = a_idx + 1;

    if (a_idx_p1 >= ntable) {
        cdist = cd_table[a_idx];
        return cdist;
    }

    if (order == 1) {
        double weight = (log_anow - log_a_table[a_idx]) /
                        (log_a_table[a_idx_p1] - log_a_table[a_idx]);
        cdist = (1.0 - weight) * cd_table[a_idx] + weight * cd_table[a_idx_p1];
    }

    else if (order == 2) {
        int a_idx_m1 = a_idx - 1;

        double ax = log_anow;
        double a0 = log_a_table[a_idx_m1];
        double a1 = log_a_table[a_idx];
        double a2 = log_a_table[a_idx_p1];

        double w0 = (ax - a1) * (ax - a2) / ((a0 - a1) * (a0 - a2));
        double w1 = (ax - a0) * (ax - a2) / ((a1 - a0) * (a1 - a2));
        double w2 = (ax - a0) * (ax - a1) / ((a2 - a0) * (a2 - a1));

        cdist = w0 * cd_table[a_idx_m1] + w1 * cd_table[a_idx] +
                w2 * cd_table[a_idx_p1];
    }

    else if (order == 3) {
        int a_idx_m1 = a_idx - 1;
        int a_idx_p2 = a_idx + 2;

        if (a_idx_p2 >= ntable) {
            double ax = log_anow;
            double a0 = log_a_table[a_idx_m1];
            double a1 = log_a_table[a_idx];
            double a2 = log_a_table[a_idx_p1];

            double w0 = (ax - a1) * (ax - a2) / ((a0 - a1) * (a0 - a2));
            double w1 = (ax - a0) * (ax - a2) / ((a1 - a0) * (a1 - a2));
            double w2 = (ax - a0) * (ax - a1) / ((a2 - a0) * (a2 - a1));

            cdist = w0 * cd_table[a_idx_m1] + w1 * cd_table[a_idx] +
                    w2 * cd_table[a_idx_p1];

        } else {
            double ax = log_anow;
            double a0 = log_a_table[a_idx_m1];
            double a1 = log_a_table[a_idx];
            double a2 = log_a_table[a_idx_p1];
            double a3 = log_a_table[a_idx_p2];

            double w0 = (ax - a1) * (ax - a2) * (ax - a3) /
                        ((a0 - a1) * (a0 - a2) * (a0 - a3));
            double w1 = (ax - a0) * (ax - a2) * (ax - a3) /
                        ((a1 - a0) * (a1 - a2) * (a1 - a3));
            double w2 = (ax - a0) * (ax - a1) * (ax - a3) /
                        ((a2 - a0) * (a2 - a1) * (a2 - a3));
            double w3 = (ax - a0) * (ax - a1) * (ax - a2) /
                        ((a3 - a0) * (a3 - a1) * (a3 - a2));

            cdist = w0 * cd_table[a_idx_m1] + w1 * cd_table[a_idx] +
                    w2 * cd_table[a_idx_p1] + w3 * cd_table[a_idx_p2];
        }
    }

    return cdist;
}

#endif

double transverse_distance(double z) { return (comoving_distance(z)); }

double angular_diameter_distance(double z) {
    return (transverse_distance(z) / (1.0 + z));
}

double luminosity_distance(double z) {
    return ((1.0 + z) * transverse_distance(z));
}

double comoving_volume_element(double z) {
    double z1da = (1.0 + z) * angular_diameter_distance(z);
    return (Dh * (z1da * z1da) / _E(z));
}

double comoving_volume(double z) {
    double r = transverse_distance(z);
    return (4.0 * M_PI * r * r * r / 3.0);
}

double comoving_distance_to_redshift(double r) {
    if (r <= 0)
        return 0;
    double z  = 1;
    double dz = 0.1;
    double rt;
    while (dz > 1e-7) {
        rt = transverse_distance(z);
        // dz = ((r - rt) * dz) /
        //      (transverse_distance(z + dz) - transverse_distance(z));
        dz = ((r - rt) * dz) / (transverse_distance(z + dz) - rt);
        if (!isfinite(dz))
            return z;
        if (z + dz < 0)
            z /= 3.0;
        else
            z += dz;
        dz = fabs(dz);
        if (dz > 0.1)
            dz = 0.1;
    }
    return z;
}

double comoving_volume_to_redshift(double Vc) {
    double r = cbrt(Vc * (3.0 / (4.0 * M_PI)));
    return (comoving_distance_to_redshift(r));
}

double comoving_distance_h_to_redshift(double r) {
    return (comoving_distance_to_redshift(r / h));
}
