#ifndef CONSTANTS_H
#define CONSTANTS_H

// FUNCTION + VARIABLE DECLARATIONS

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double self_energy(double);
void init_rod(int);
double mc_sweep(gsl_rng *r);
int metropolis_step(gsl_rng *r);
double compute_dE(int, int, double);
double compute_eta();
double compute_square_gradient(int, int, double);
double interaction_energy(int, int, double);
double compute_total_energy();
const double pi = atan(1.0)*4;
const double kB = 1.3806503 * 6.0221415 / 4184.0;
double deg_to_rad(double th) {
  return th * pi / 180;
}

double rad_to_deg(double th) {
  return th * 180 / pi;
}
double T0 = 375.6, temp = 360.0, beta = 1/(kB*360.0), deta_max = 5.0;
double a = 0.305, d0 = 0.0, c = 6.99, kappa = 0.321;
double E_tot = 0.0, eta_tot = 0.0;
// double th_max = deg_to_rad(45), th_min = deg_to_rad(7);
// double var_y = sd_y * sd_y;
const int  N = 20, n_steps = 100000;   // N is number of lattice sites in each dimension
double lattice[N][N];

#endif
