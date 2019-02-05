// Monte Carlo simulations for ligands on nanoparticle surfaces
// using a Landau-Ginzburg type Hamiltonian

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include "constants.h" 

// MAIN RUN

int main(int argc, char **argv) {
  using namespace std;
  double eta_av=0.0, result=0.0;

  // get temperature from command line
  if (argc != 3) {
    cout << "USAGE:\n >$ ./test.o -t 300\n";
    exit(0);
  }
  else {
    for (int i = 1; i < argc; i += 2) {
      if (strcmp(argv[i],"-t") == 0) {
        temp = atof(argv[i+1]);
      } else {
        cout << "USAGE:\n >$ ./test.o -t 300\n";
        exit(0);
      }
    }
  }

  ostringstream shitcppformat;
  shitcppformat << "traj-" << temp << ".txt";
  string filname = shitcppformat.str();
  ofstream fout(filname.c_str());
  beta = 1 / (kB * temp);
  gsl_rng *r;
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r, time(NULL));

  init_rod(1);

  for (int i = 0; i < n_steps; i++) {
    result = mc_sweep(r);

    // print out eta_tot at timestep 'i'
    fout << i << " " << eta_tot << "\n";

    // update sum
    eta_av += eta_tot;
  }

  fout << n_steps << " " << eta_tot << "\n";
  eta_av /= n_steps;
  cout << "Average theta_z for the run was: " << eta_av << "\n"; 

  fout.close(); 
  return 0;
}

// FUNCTION DEFINITIONS

// signe denotes the sign of the spins in initial configuration (i.e., all up or down)
void init_rod(int signe){ 
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      lattice[i][j] = 1.0;
      if (signe == -1) {
        lattice[i][j] *= -1; 
      }
      E_tot += self_energy(lattice[i][j]) + interaction_energy(i, j, 0.0);
    }
  }
}

double compute_square_gradient(int i, int j, double deta) {
  double deta_dx, deta_dz, deta_dx2, deta_dz2;
  int im1, jm1, ip1, jp1;
  if (i == 0) {
    im1 = N - 1;
  }  else im1 = i - 1;
  if (i == (N - 1)) {
    ip1 = 0;
  }  else ip1 = i + 1;
  if (j == 0) {
    jm1 = N - 1;
  }  else jm1 = j - 1;
  if (j == (N - 1)) {
    jp1 = 0;
  }  else jp1 = j + 1;
  
  deta_dx = (lattice[im1][j] - lattice[ip1][j] - deta) / 2.0;
  deta_dz = (lattice[i][jm1] - lattice[i][jp1] - deta) / 2.0;
  deta_dx2 = deta_dx * deta_dx; 
  deta_dz2 = deta_dz * deta_dz; 
  double grad_sq = deta_dx2 + deta_dz2;
  return grad_sq;
}

double self_energy(double eta) {
  double se = 0;
  double eta2 = eta * eta;
  double eta3 = eta * eta * eta;
  double eta4 = eta2 * eta2;
  se = a * (temp - T0) * eta2 / 2.0 + d0 * (temp - T0) * eta3 / 3.0 + c * eta4 / 4.0;
  return se; 
}

double interaction_energy(int i, int j, double deta) {
  double ie = kappa * compute_square_gradient(i, j, deta) / 2.0;
  return ie;
}

double mc_sweep(gsl_rng * r) {
  double result = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      result += metropolis_step(r);
    }
  }
  return result;
}

double compute_eta() {
  double eta = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      eta += lattice[i][j];
    }
  }
  return eta/(N*N);
}

int metropolis_step(gsl_rng *r) {
  int i, j;
  // choose a random spin to flip
  i = gsl_rng_uniform_int(r, N);
  j = gsl_rng_uniform_int(r, N);
  double eta_diff = (gsl_rng_uniform(r) - 0.5) * deta_max;	// generate random number between -0.5 and 0.5; scale by deta_max

  // calculate  energy change for spin flip
  double dE = compute_dE(i, j, eta_diff);

  // flip spin if it's favourable
  if (dE < 0.0) {
    lattice[i][j] += eta_diff;
    E_tot += dE;
    return 1;
  }
  else {
    double rand = gsl_rng_uniform(r);
    if (exp(-dE*beta) > rand) {
      lattice[i][j] += eta_diff;
      E_tot += dE;
      eta_tot = compute_eta();
      return 1;
    }
    else {
      return 0;
    }
  }
}

double compute_dE(int i, int j, double deta) {
  double dE = 0.0;
  double d_self_en = self_energy(lattice[i][j] + deta) - self_energy(lattice[i][j]);
  double d_int_en = 0.0;

  // IMPLEMENT PBC's
  int im1, jm1, ip1, jp1;
  if (i == 0) {
    im1 = N - 1;
  }  else im1 = i - 1;
  if (i == (N - 1)) {
    ip1 = 0;
  }  else ip1 = i + 1;
  if (j == 0) {
    jm1 = N - 1;
  }  else jm1 = j - 1;
  if (j == (N - 1)) {
    jp1 = 0;
  }  else jp1 = j + 1;

  d_int_en += ( interaction_energy(ip1, j, deta) - interaction_energy(ip1, j, 0.0) );
  d_int_en += ( interaction_energy(im1, j, -deta) - interaction_energy(im1, j, 0.0) );
  d_int_en += ( interaction_energy(i, jp1, deta) - interaction_energy(i, jp1, 0.0) );
  d_int_en += ( interaction_energy(i, jm1, -deta) - interaction_energy(i, jm1, 0.0) );

  dE = d_self_en + d_int_en;
  return dE;
}

// find total energy of the system
double compute_total_energy() {
  double total_energy = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      total_energy += self_energy(lattice[i][j]) + interaction_energy(i, j, 0.0);
    }
  }
  return total_energy;
}

