#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <inttypes.h>
#include "rockstar.h"
#include "halo.h"
#include "fof.h"
#include "particle.h"
#include "groupies.h"
#include "subhalo_metric.h"
#include "check_syscalls.h"
#include "config_vars.h"
#include "universal_constants.h"
#include "potential.h"
#include "nfw.h"
#include "distance.h"
#include "fun_times.h"
#include "jacobi.h"
#include "hubble.h"
#include "vir_dense.h"

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#define EPSILON 1.0E-6

double particle_thresh_dens[5] = {0};
// #pragma omp threadprivate(particle_thresh_dens) // to groupies.h

gsl_interp_accel *acc, *acc2, *acc3, *acc4;
gsl_spline       *spl, *spl2, *spl3, *spl4;

gsl_spline       *spl_hubble, *spl_adash;
gsl_interp_accel *acc_hubble, *acc_adash;

gsl_spline       *spl_trans;
gsl_interp_accel *acc_trans;

double get_ascale_from_linear_growth(double dtarget);

double F_Omega(double a) {
    double value = gsl_spline_eval(spl, a, acc);
    double deriv = gsl_spline_eval_deriv(spl, a, acc);
    return deriv * a / value;
}

double F2_Omega(double a) {
    double dplus, dminus, d;
    double aplus, aminus;

    aplus  = exp(log(a) + 0.001);
    aminus = exp(log(a) - 0.001);
    dplus  = get_2nd_growth(aplus);
    dminus = get_2nd_growth(aminus);
    d      = get_2nd_growth(a);
    return (dplus - dminus) / (log(aplus) - log(aminus)) / d;
}

double Hubble_ratio_a(double a) {
    double z = 1.0 / a - 1.0;
    return hubble_scaling(z);
}

double Hubble_deriv_lna(double lna) {
    double lna_plus  = lna + 0.001;
    double lna_minus = lna - 0.001;
    double aplus     = exp(lna_plus);
    double aminus    = exp(lna_minus);
    if (aplus > 1.) {
        aplus    = 1.;
        lna_plus = 0.;
    }
    return (Hubble_ratio_a(aplus) - Hubble_ratio_a(aminus)) /
           (lna_plus - lna_minus);
}

double get_w(double a) { return W0 + WA * (1. - a); }

double Omegam_a(double a) {
    return Om / (pow(Hubble_ratio_a(a), 2.) * a * a * a);
}

double Omegade_a(double a) {
    return Ol / (pow(Hubble_ratio_a(a), 2.) * pow(a, 3. * (1. + W0 + WA)) *
                 exp(3. * WA * (1. - a)));
}

void display_results(char *title, double result, double error) {
    printf("%s ==================\n", title);
    printf("result = % .6f\n", result);
    printf("sigma = % .6f\n", error);
    fflush(stdout);
}

double get_linear_growth(double a) {
    if (a < 1e-3)
        return a;
    else
        return gsl_spline_eval(spl, a, acc);
}

double get_linear_growth_test(double a) {
    if (a < 1e-3)
        return a;
    else
        return gsl_spline_eval(spl3, a, acc3);
}

double get_2nd_growth(double a) {
    if (a < 1e-3)
        return -3 / 7. * a * a;
    else
        return gsl_spline_eval(spl2, a, acc2);
}

double get_ascale_from_linear_growth(double dtarget) {
    double amin, amax, atry, dtry;
    amin = 0.;
    amax = 1.;
    atry = (amin + amax) / 2.;
    while ((amax - amin) / atry > 1e-7) {
        atry = 0.5 * (amin + amax);
        dtry = get_linear_growth(atry);
        if (dtry > dtarget)
            amax = atry;
        else if (dtry < dtarget)
            amin = atry;
        else
            return atry;
    }
    return atry;
}

// TN added a new function for D(z): consider the small scale limit, where only
// C+B fluid contributes the structure growth. Other fluids (dark energy,
// neutrino etc.) contribute only through the expansion rate, which is given by
// the Hubble table.
void initialize_linear_growth() {
    double lna0, lna, F, G, dlna;
    double Ftmp1, Ftmp2, Ftmp3, Ftmp4;
    double Gtmp1, Gtmp2, Gtmp3, Gtmp4;
    double P0, Q0, om0, lnHderiv0;
    double P, Q, om, lnHderiv;
    double P1, Q1, om1, lnHderiv1;
    double la0, la, la1;
    double anow0, anow, anow1;
    int    nint, n;

    lna0 = log(1e-3);
    lna  = log(1.);
    F    = exp(lna0);
    G    = exp(lna0);
    dlna = 0.0001;
    nint = (lna - lna0) / dlna;
    dlna = (lna - lna0) / nint;

    acc = gsl_interp_accel_alloc();
    spl = gsl_spline_alloc(gsl_interp_cspline, nint + 2);

    double *ds, *as;
    ds = (double *)malloc(sizeof(double) * (nint + 2));
    as = (double *)malloc(sizeof(double) * (nint + 2));

    as[0] = 0.;
    ds[0] = 0.;

    as[1] = exp(lna0);
    ds[1] = as[1];

    for (n = 0; n < nint; n++) {
        la0       = lna0 + n * dlna;
        la        = lna0 + (n + 0.5) * dlna;
        la1       = lna0 + (n + 1) * dlna;
        anow0     = exp(la0);
        anow      = exp(la);
        anow1     = exp(la1);
        om0       = Omegam_a(anow0);
        om        = Omegam_a(anow);
        om1       = Omegam_a(anow1);
        lnHderiv0 = Hubble_deriv_lna(la0) / Hubble_ratio_a(anow0);
        lnHderiv  = Hubble_deriv_lna(la) / Hubble_ratio_a(anow);
        lnHderiv1 = Hubble_deriv_lna(la1) / Hubble_ratio_a(anow1);

        P0 = -(lnHderiv0 + 2.);
        Q0 = 3 / 2. * om0;
        P  = -(lnHderiv + 2.);
        Q  = 3 / 2. * om;
        P1 = -(lnHderiv1 + 2.);
        Q1 = 3 / 2. * om1;

        Ftmp1 = (P0 * F + Q0 * G) * dlna;
        Gtmp1 = F * dlna;
        Ftmp2 = (P * (F + Ftmp1 / 2.) + Q * (G + Gtmp1 / 2.)) * dlna;
        Gtmp2 = (F + Ftmp1 / 2.) * dlna;
        Ftmp3 = (P * (F + Ftmp2 / 2.) + Q * (G + Gtmp2 / 2.)) * dlna;
        Gtmp3 = (F + Ftmp2 / 2.) * dlna;
        Ftmp4 = (P1 * (F + Ftmp3) + Q1 * (G + Gtmp3)) * dlna;
        Gtmp4 = (F + Ftmp3) * dlna;

        F += (Ftmp1 + 2. * Ftmp2 + 2. * Ftmp3 + Ftmp4) / 6.;
        G += (Gtmp1 + 2. * Gtmp2 + 2. * Gtmp3 + Gtmp4) / 6.;
        as[n + 2] = anow1;
        ds[n + 2] = G;
    }
    gsl_spline_init(spl, as, ds, nint + 2);

    free(ds);
    free(as);
}

/*
void initialize_linear_growth()
{
  // printf("Initializing the linear growth ...\n");
  // fflush(stdout);
  double lna0, lna, F, G, dlna;
  double Ftmp1, Ftmp2, Ftmp3, Ftmp4;
  double Gtmp1, Gtmp2, Gtmp3, Gtmp4;
  double P0, Q0, om0, ode0;
  double P, Q, om, ode;
  double P1, Q1, om1, ode1;
  double anow0, anow, anow1;
  int nint, n;

  lna0 = log(1e-6);
  lna = log(2.5);
  F = 0.;
  G = 1.;
  dlna = 0.0001;
  nint = (lna - lna0)/dlna;
  dlna = (lna - lna0)/nint;

  acc = gsl_interp_accel_alloc();
  spl = gsl_spline_alloc(gsl_interp_cspline,nint+2);

  //double ds[nint+2];
  //double as[nint+2];
  double *ds, *as;
  ds = (double*)malloc(sizeof(double)*(nint+2));
  as = (double*)malloc(sizeof(double)*(nint+2));

  as[0] = 0.;
  ds[0] = 0.;

  as[1] = exp(lna0);
  ds[1] = as[1];

  for(n = 0; n<nint; n++){
    anow0 = exp(lna0 + n*dlna);
    anow = exp(lna0 + (n+0.5)*dlna);
    anow1 = exp(lna0 + (n+1)*dlna);
    om0 = Omegam_a(anow0);
    ode0 = Omegade_a(anow0);
    om = Omegam_a(anow);
    ode = Omegade_a(anow);
    om1 = Omegam_a(anow1);
    ode1 = Omegade_a(anow1);

    P0 = 3. - 0.5 * (om0 + (1.+3.*get_w(anow0)) * ode0);
    Q0 = 2. - 2. * om0 - 0.5 * (1.+3.*get_w(anow0)) * ode0;
    P = 3. - 0.5 * (om + (1.+3.*get_w(anow)) * ode);
    Q = 2. - 2. * om - 0.5 * (1.+3.*get_w(anow)) * ode;
    P1 = 3. - 0.5 * (om1 + (1.+3.*get_w(anow1)) * ode1);
    Q1 = 2. - 2. * om1 - 0.5 * (1.+3.*get_w(anow1)) * ode1;

    Ftmp1 = - (P0 * F + Q0 * G) * dlna;
    Gtmp1 = F * dlna;
    Ftmp2 = - (P * (F+Ftmp1/2.) + Q * (G+Gtmp1/2.)) * dlna;
    Gtmp2 = (F+Ftmp1/2.) * dlna;
    Ftmp3 = - (P * (F+Ftmp2/2.) + Q * (G+Gtmp2/2.)) * dlna;
    Gtmp3 = (F+Ftmp2/2.) * dlna;
    Ftmp4 = - (P1 * (F+Ftmp3) + Q1 * (G+Gtmp3)) * dlna;
    Gtmp4 = (F+Ftmp3) * dlna;
    F += (Ftmp1 + 2. * Ftmp2 + 2. * Ftmp3 + Ftmp4) / 6.;
    G += (Gtmp1 + 2. * Gtmp2 + 2. * Gtmp3 + Gtmp4) / 6.;
    as[n+2] = anow1;
    ds[n+2] = anow1*G;
  }
  gsl_spline_init(spl,as,ds,nint+2);

  free(ds);
  free(as);
  // return a * G;
  // printf(" done.\n");
}
*/

double spherical_collapse(double a, double deltab0) {
    double eta0, eta, E, F, deta;
    double Ftmp1, Ftmp2, Ftmp3, Ftmp4;
    double Etmp1, Etmp2, Etmp3, Etmp4;
    double om, f;
    double kappa0, kappa, kappa1;
    double enow0, enow, enow1;
    double anow;
    int    nint, n;

    if (a < 1e-3)
        return deltab0 * get_linear_growth(a) / get_linear_growth(1.);

    eta0 = log(1e-3);
    eta  = log(get_linear_growth(a));
    E    = deltab0 * exp(eta0) / get_linear_growth(1.);
    F    = E;
    deta = 0.0001;
    nint = (eta - eta0) / deta;
    deta = (eta - eta0) / nint;

    for (n = 0; n < nint; n++) {
        enow0 = eta0 + n * deta;
        enow  = eta0 + (n + 0.5) * deta;
        enow1 = eta0 + (n + 1) * deta;
        anow  = get_ascale_from_linear_growth(exp(enow0));
        // printf("%g  %g  %g\n",anow,E,pow(1.+E,-1/3.)*anow);
        om     = Omegam_a(anow);
        f      = F_Omega(anow);
        kappa0 = 1.5 * om / (f * f);
        anow   = get_ascale_from_linear_growth(exp(enow));
        om     = Omegam_a(anow);
        f      = F_Omega(anow);
        kappa  = 1.5 * om / (f * f);
        anow   = get_ascale_from_linear_growth(exp(enow1));
        om     = Omegam_a(anow);
        f      = F_Omega(anow);
        kappa1 = 1.5 * om / (f * f);

        Ftmp1 = (-(kappa0 - 1.) * F + 4. / 3 * F * F / (1. + E) +
                 kappa0 * (1 + E) * E) *
                deta;
        Etmp1 = F * deta;
        Ftmp2 = (-(kappa - 1.) * (F + Ftmp1 / 2.) +
                 4. / 3 * (F + Ftmp1 / 2.) * (F + Ftmp1 / 2.) /
                     (1. + (E + Etmp1 / 2.)) +
                 kappa * (1 + E + Etmp1 / 2.) * (E + Etmp1 / 2.)) *
                deta;
        Etmp2 = (F + Ftmp2 / 2.) * deta;
        Ftmp3 = (-(kappa - 1.) * (F + Ftmp2 / 2.) +
                 4. / 3 * (F + Ftmp2 / 2.) * (F + Ftmp2 / 2.) /
                     (1. + (E + Etmp2 / 2.)) +
                 kappa * (1 + E + Etmp2 / 2.) * (E + Etmp2 / 2.)) *
                deta;
        Etmp3 = (F + Ftmp2 / 2.) * deta;
        Ftmp4 = (-(kappa1 - 1.) * (F + Ftmp3) +
                 4. / 3 * (F + Ftmp3) * (F + Ftmp3) / (1. + (E + Etmp3)) +
                 kappa1 * (1 + E + Etmp3) * (E + Etmp3)) *
                deta;
        Etmp4 = (F + Ftmp2 / 2.) * deta;

        F += (Ftmp1 + 2. * Ftmp2 + 2. * Ftmp3 + Ftmp4) / 6.;
        E += (Etmp1 + 2. * Etmp2 + 2. * Etmp3 + Etmp4) / 6.;
    }
    return E;
}

void calc_a_ta_zeta(double a, double deltab0, double *a_ta, double *zeta) {
    double eta0, eta, E, F, deta;
    double Ftmp1, Ftmp2, Ftmp3, Ftmp4;
    double Etmp1, Etmp2, Etmp3, Etmp4;
    double om, f;
    double kappa0, kappa, kappa1;
    double enow0, enow, enow1;
    double anow;
    int    nint, n;

    double radius_now;
    double radius_max = 0.;
    double a_ta_0     = 0.;
    double E_ta       = 0.;

    // if(a<1e-7) return deltab0 * get_linear_growth(a)/get_linear_growth(1.);

    eta0 = log(1e-3);
    eta  = log(get_linear_growth(a));
    E    = deltab0 * exp(eta0) / get_linear_growth(1.);
    F    = E;
    deta = 0.0001;
    nint = (eta - eta0) / deta;
    deta = (eta - eta0) / nint;

    for (n = 0; n < nint; n++) {
        enow0 = eta0 + n * deta;
        enow  = eta0 + (n + 0.5) * deta;
        enow1 = eta0 + (n + 1) * deta;
        anow  = get_ascale_from_linear_growth(exp(enow0));
        // printf("%g  %g  %g\n",anow,E,pow(1.+E,-1/3.)*anow);
        om     = Omegam_a(anow);
        f      = F_Omega(anow);
        kappa0 = 1.5 * om / (f * f);
        anow   = get_ascale_from_linear_growth(exp(enow));
        om     = Omegam_a(anow);
        f      = F_Omega(anow);
        kappa  = 1.5 * om / (f * f);
        anow   = get_ascale_from_linear_growth(exp(enow1));
        om     = Omegam_a(anow);
        f      = F_Omega(anow);
        kappa1 = 1.5 * om / (f * f);

        Ftmp1 = (-(kappa0 - 1.) * F + 4. / 3 * F * F / (1. + E) +
                 kappa0 * (1 + E) * E) *
                deta;
        Etmp1 = F * deta;
        Ftmp2 = (-(kappa - 1.) * (F + Ftmp1 / 2.) +
                 4. / 3 * (F + Ftmp1 / 2.) * (F + Ftmp1 / 2.) /
                     (1. + (E + Etmp1 / 2.)) +
                 kappa * (1 + E + Etmp1 / 2.) * (E + Etmp1 / 2.)) *
                deta;
        Etmp2 = (F + Ftmp2 / 2.) * deta;
        Ftmp3 = (-(kappa - 1.) * (F + Ftmp2 / 2.) +
                 4. / 3 * (F + Ftmp2 / 2.) * (F + Ftmp2 / 2.) /
                     (1. + (E + Etmp2 / 2.)) +
                 kappa * (1 + E + Etmp2 / 2.) * (E + Etmp2 / 2.)) *
                deta;
        Etmp3 = (F + Ftmp2 / 2.) * deta;
        Ftmp4 = (-(kappa1 - 1.) * (F + Ftmp3) +
                 4. / 3 * (F + Ftmp3) * (F + Ftmp3) / (1. + (E + Etmp3)) +
                 kappa1 * (1 + E + Etmp3) * (E + Etmp3)) *
                deta;
        Etmp4 = (F + Ftmp2 / 2.) * deta;

        F += (Ftmp1 + 2. * Ftmp2 + 2. * Ftmp3 + Ftmp4) / 6.;
        E += (Etmp1 + 2. * Etmp2 + 2. * Etmp3 + Etmp4) / 6.;

        radius_now = pow(1. + E, -1 / 3.) * anow;
        // printf("anow = %g, radius_now = %g, a_ta = %g, E = %g\n",anow,
        // radius_now, a_ta, E);
        if (radius_now > radius_max) {
            radius_max = radius_now;
            a_ta_0     = anow;
            E_ta       = E;
        } else
            break;
    }
    // printf("anow = %g, radius_now = %g, a_ta = %g, E = %g\n",anow,
    // radius_now, a_ta, E);
    *a_ta = a_ta_0;
    *zeta = 1 + E_ta;
}

double Newton_method(double eta_vir, double eta_ta) {

    double ans, ans1;

    ans  = 0.0;
    ans1 = next(ans, eta_vir, eta_ta);
    while (fabs(ans - ans1) > EPSILON) {
        // printf("Value now(Newton method) = %1.10f\n", fabs(ans-ans1));
        ans1 = ans;
        ans  = next(ans1, eta_vir, eta_ta);
        //      printf("Value now(Newton method) = %1.10f\n", fabs(ans-ans1));
    }
    // printf("y = %1.10f\n", ans);
    return ans;
}

double next(double x, double eta_vir, double eta_ta) {

    return (x - (2. * eta_vir * x * x * x - (2. + eta_ta) * x + 1.) /
                    (6. * eta_vir * x * x - 2. - eta_ta));
}
