/* Courtesy of Matt Becker. */
#include <math.h>
#include <assert.h>
#include "config_vars.h"
#include "hubble.h"

/* See http://arxiv.org/pdf/astro-ph/0208512v1.pdf */
double _weff(double a) {
    if (a != 1.0)
        return W0 + WA - WA * (a - 1.0) / log(a);
    else
        return W0;
}

// switch for GINKAKU
// add by stanaka
#ifndef FOR_GINKAKU

double hubble_scaling(double z) {
    double z1 = 1.0 + z;
    double a  = 1.0 / z1;
    return sqrt(Om * (z1 * z1 * z1) + Ol * pow(a, -3.0 * (1.0 + _weff(a))));
}

#else

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int    set_expansion_table(char *, double);
double read_expansion_table(const int, const double, const int);
double read_expansion_table_log(const int, const double, const int);

double hubble_table[2][2048];
int    hubble_table_size;

/* H(z)/H0 */
double hubble_scaling(double z) {
    double z1 = 1.0 + z;
    double a  = 1.0 / z1;

    static int set_table = 0;
    static int ltable    = 0;

    if (strcmp(EXPANSION_TABLEFILE, "none") == 0) {
        return sqrt(Om * (z1 * z1 * z1) + Ol * pow(a, -3.0 * (1.0 + _weff(a))));
    }

    // const double amin = 1.0e-2;  // z ~ 100
    const double amin =
        1.0e-6; // smaller than minimum ascale of initialize_linear_growth()

    if (set_table == 0) {
        ltable            = set_expansion_table(EXPANSION_TABLEFILE, amin);
        hubble_table_size = ltable;
        printf("read expansion table %s. table length=%d\n",
               EXPANSION_TABLEFILE, ltable);
        fflush(stdout);
        set_table = 1;
    }

    // exception handling for first step
    if (a < amin || a > 0.9999)
        return 1.0;

    const int order = READ_TABLE_ORDER;
    return read_expansion_table_log(order, a, ltable);
}

/* (H(z)/H0) */
int set_expansion_table(char *inputfile, const double amin) {
    FILE *fp;
    fp = fopen(inputfile, "r");

    if (fp == NULL) {
        printf("Error: %s file not opened.", inputfile);
        return EXIT_FAILURE;
    }

    double atbl0 = 0.9 * amin;
    double a, h;

    int line = 0;

    static char tmp[128];
    size_t      len = 128;

    // skip header
    fgets(tmp, len, fp);

    while (fscanf(fp, "%lf %lf", &a, &h) != EOF) {
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

    while (fscanf(fp, "%lf %lf", &a, &h) != EOF) {
        if (atbl0 > a)
            continue;
        // hubble_table[0][line] = a;
        hubble_table[0][line] = log10(a);
        hubble_table[1][line] = h;
        line++;
    }

    fclose(fp);

    return line;
}

double read_expansion_table(const int order, const double anow,
                            const int ntable) {
    double *a_table, *h_table;
    a_table = hubble_table[0]; // ascale
    h_table = hubble_table[1]; // H0/h

    double a_min = a_table[0];
    // double a_max = a_table[ntable - 1];
    double da    = log10(a_table[1]) - log10(a_table[0]);
    int    a_idx = (log10(anow) - log10(a_min)) / da;

    static int lsearch  = 3;
    int        hit_flag = 0;

    for (int idx = a_idx - lsearch; idx < a_idx + lsearch; idx++) {
        if (idx + 1 >= ntable)
            continue;
        if (a_table[idx] <= anow && anow <= a_table[idx + 1]) {
            a_idx    = idx;
            hit_flag = 1;
            break;
        }
    }

    if (hit_flag == 0) {
        for (int idx = 0; idx < ntable - 1; idx++) {
            if (a_table[idx] <= anow && anow <= a_table[idx + 1]) {
                a_idx = idx;
                break;
            }
        }
        lsearch++;
    }

    double hexp     = 0.0;
    int    a_idx_p1 = a_idx + 1;

    if (a_idx_p1 >= ntable) {
        hexp = h_table[a_idx];
        return hexp;
    }

    if (order == 1) {
        double weight = (log10(anow) - log10(a_table[a_idx])) /
                        (log10(a_table[a_idx_p1]) - log10(a_table[a_idx]));
        hexp = (1.0 - weight) * h_table[a_idx] + weight * h_table[a_idx_p1];
    }

    else if (order == 2) {
        int a_idx_m1 = a_idx - 1;

        double ax = log10(anow);
        double a0 = log10(a_table[a_idx_m1]);
        double a1 = log10(a_table[a_idx]);
        double a2 = log10(a_table[a_idx_p1]);

        double w0 = (ax - a1) * (ax - a2) / ((a0 - a1) * (a0 - a2));
        double w1 = (ax - a0) * (ax - a2) / ((a1 - a0) * (a1 - a2));
        double w2 = (ax - a0) * (ax - a1) / ((a2 - a0) * (a2 - a1));

        hexp = w0 * h_table[a_idx_m1] + w1 * h_table[a_idx] +
               w2 * h_table[a_idx_p1];
    }

    else if (order == 3) {
        int a_idx_m1 = a_idx - 1;
        int a_idx_p2 = a_idx + 2;

        if (a_idx_p2 >= ntable) {
            double ax = log10(anow);
            double a0 = log10(a_table[a_idx_m1]);
            double a1 = log10(a_table[a_idx]);
            double a2 = log10(a_table[a_idx_p1]);

            double w0 = (ax - a1) * (ax - a2) / ((a0 - a1) * (a0 - a2));
            double w1 = (ax - a0) * (ax - a2) / ((a1 - a0) * (a1 - a2));
            double w2 = (ax - a0) * (ax - a1) / ((a2 - a0) * (a2 - a1));

            hexp = w0 * h_table[a_idx_m1] + w1 * h_table[a_idx] +
                   w2 * h_table[a_idx_p1];

        } else {
            double ax = log10(anow);
            double a0 = log10(a_table[a_idx_m1]);
            double a1 = log10(a_table[a_idx]);
            double a2 = log10(a_table[a_idx_p1]);
            double a3 = log10(a_table[a_idx_p2]);

            double w0 = (ax - a1) * (ax - a2) * (ax - a3) /
                        ((a0 - a1) * (a0 - a2) * (a0 - a3));
            double w1 = (ax - a0) * (ax - a2) * (ax - a3) /
                        ((a1 - a0) * (a1 - a2) * (a1 - a3));
            double w2 = (ax - a0) * (ax - a1) * (ax - a3) /
                        ((a2 - a0) * (a2 - a1) * (a2 - a3));
            double w3 = (ax - a0) * (ax - a1) * (ax - a2) /
                        ((a3 - a0) * (a3 - a1) * (a3 - a2));

            hexp = w0 * h_table[a_idx_m1] + w1 * h_table[a_idx] +
                   w2 * h_table[a_idx_p1] + w3 * h_table[a_idx_p2];
        }
    }

    return hexp;
}

double read_expansion_table_log(const int order, const double anow,
                                const int ntable) {
    double *log_a_table, *h_table;
    log_a_table = hubble_table[0]; // log10(ascale)
    h_table     = hubble_table[1]; // H0/h

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

    double hexp     = 0.0;
    int    a_idx_p1 = a_idx + 1;

    if (a_idx_p1 >= ntable) {
        hexp = h_table[a_idx];
        return hexp;
    }

    if (order == 1) {
        double weight = (log_anow - log_a_table[a_idx]) /
                        (log_a_table[a_idx_p1] - log_a_table[a_idx]);
        hexp = (1.0 - weight) * h_table[a_idx] + weight * h_table[a_idx_p1];
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

        hexp = w0 * h_table[a_idx_m1] + w1 * h_table[a_idx] +
               w2 * h_table[a_idx_p1];
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

            hexp = w0 * h_table[a_idx_m1] + w1 * h_table[a_idx] +
                   w2 * h_table[a_idx_p1];

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

            hexp = w0 * h_table[a_idx_m1] + w1 * h_table[a_idx] +
                   w2 * h_table[a_idx_p1] + w3 * h_table[a_idx_p2];
        }
    }

    return hexp;
}

#endif
