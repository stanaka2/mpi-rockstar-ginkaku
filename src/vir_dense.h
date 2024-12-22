#include <stdint.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#define EPSILON 1.0E-6

extern gsl_interp_accel *acc, *acc2, *acc3, *acc4;
extern gsl_spline *spl, *spl2, *spl3, *spl4;

extern gsl_spline *spl_hubble, *spl_adash;
extern gsl_interp_accel *acc_hubble, *acc_adash;

extern gsl_spline *spl_trans;
extern gsl_interp_accel *acc_trans;

int FatalError(int errnum);
void read_transfer_table(void);
double TransferFunction(double kmag);
void initialize_transferfunction(void);
void initialize_linear_growth();

double F_Omega(double);
double F2_Omega(double);
double calc_time(double);
double time_int(double);
double time_glo_int(double);
double get_linear_growth(double);
double get_2nd_growth(double);
double get_linear_growth_test(double);

double Hubble_ratio_a(double);
double Hubble_deriv_lna(double);
double Omegam_a(double);
double Omegade_a(double);

double get_w(double);

void spherical_collapse_table(double);
double spherical_collapse(double, double);
void write_parameter_files();

void calc_local_expansion(double);
void read_parameterfile(char *);

void set_separate_universe(void);
void initialize_local_box(void);
double growth(double);
double get_delta_a(double);

void calc_a_ta_zeta(double, double, double *, double *);
double Newton_method(double, double);
double next(double, double, double);
