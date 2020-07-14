// Monte Carlo simulations for ligands on nanoparticle surfaces
// using a Landau-Ginzburg type Hamiltonian

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <string.h>
#include <sstream>
#include "constants.h" 

// MAIN RUN

int main(int argc, char **argv) {
  using namespace std;
  double eta_av=0.0, result=0.0, acc_prob=0.0;

  // get temperature from command line
  if (argc != 5) {
    cout << "USAGE:\n >$ ./do_mc.o -t 300 -bc_on 0\n";
    exit(0);
  }
  else {
    for (int i = 1; i < argc; i += 2) {
      if (strcmp(argv[i],"-t") == 0) {
        temp = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-bc_on") == 0) {
        bound_on = atof(argv[i+1]);
        if (bound_on != 0 && bound_on != 1) {
          cout << "bc_on has to be either 0 or 1!!\n";
          exit(0);
        }
      } else {
        cout << "USAGE:\n >$ ./do_mc.o -t 300 -bc_on 1\n";
        exit(0);
      }
    }
  }

  if (bound_on == 1) {
    h = 0.0;
  }

  ostringstream datafilstring, trajfilstring;
  datafilstring << "log-" << temp << ".txt";
  string filname = datafilstring.str();
  ofstream fout(filname.c_str());
  trajfilstring << "dump-" << temp << ".xyz";
  string trajfilname = trajfilstring.str();
  ofstream dumpout(trajfilname.c_str());

  beta = 1 / (kB * temp);
  gsl_rng *r;
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r, time(NULL));
//   gsl_rng_set(r, 5743);

  init_rod(1);
  calc_gradient_matrix();
  E_tot = compute_total_energy();
  int data_rate = max(1, n_steps / 10000);
  int screen_rate = max(1, n_steps / 1000);
  int frame_rate = max(1, n_steps / 10000);
//   cout << "frame rate is: " << frame_rate << "\n";
//   cout << "data rate is: " << data_rate << "\n";
//   cout << "screen rate is: " << screen_rate << "\n";

  // output format
  fout << "sweep no (1)  acc prob (2)  total energy (3)  e_total variable (4) self energy (5)  interaction energy (6)  mean eta (7)  mean eta^2 (8)\n";
  cout << "sweep no (1)  acc prob (2)  total energy (3)  e_total variable (4) self energy (5)  interaction energy (6)  mean eta (7)  mean eta^2 (8)\n";

  for (int i = 0; i < n_steps; i++) {
    // print out output to log file to get 10000 entries
    if (i % data_rate == 0) {    
      fout << i << "  " << acc_prob << "  " << compute_total_energy() << "  "  << E_tot << "  " << get_self_energy() << "  " << get_interaction_energy() << "  " << compute_eta() << "  " << compute_eta2() << "\n";
    }

    // print out output to screen to get 1000 entries
    if (i % screen_rate == 0) {    
      cout << i << "  " << acc_prob << "  " << compute_total_energy() << "  "  << E_tot << "  " << get_self_energy() << "  " << get_interaction_energy() << "  " << compute_eta() << "  " << compute_eta2() << "\n";
    }

    // print out configuration to dump file to get 100 entries
    if (i % frame_rate == 0 && i != 0) {
      dumpout << N*N << "\n";
      dumpout << "Comment: box size " << N << "x" << N  << " at timestep: " << i << "\n"; 
      for (int j = 0; j<N; j++) {
        for (int k = 0; k<N; k++) {
        dumpout << "1  " << j << "  " << k << "  " << lattice[j][k] << "\n";
        }
      }
    }

    // do MC moves
    result += mc_sweep(r);
    acc_prob = result / (N * N * (i+1));

    // update sum
    eta_av += compute_eta();
  }

//   fout << n_steps << " " << eta_tot << "\n";
  eta_av /= n_steps;
  cout << "Average theta_z for the run was: " << eta_av << "\n"; 
  dumpout << N*N << "\n";
  dumpout << "Comment: box size " << N << "x" << N  << " at timestep: " << n_steps << "\n"; 
  for (int j = 0; j<N; j++) {
    for (int k = 0; k<N; k++) {
    dumpout << "1  " << j << "  " << k << "  " << lattice[j][k] << "\n";
    }
  }

  fout.close();
  dumpout.close(); 
  return 0;
}

// FUNCTION DEFINITIONS

// signe denotes the sign of the spins in initial configuration (i.e., all up or down)
void init_rod(int signe){ 
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      lattice[i][j] = 0.0;
      if (bound_on == 1) {
        // set boundary values
        if (j == 0 or j == 1) {
          lattice[i][j] = -0.3;
        }
        if (j == N-1 or j == N-2) {
          lattice[i][j] = +0.3;
        }
      }
//       lattice[i][j] = 1.0;
//       lattice[i][j] = j*0.1;
      if (signe == -1) {
        lattice[i][j] *= -1;
      }
    }
  }
}

double compute_gradient(int i, int j, char dir) {
  double deta_dx, deta_dz;
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

  if (dir == 'x') { 
    deta_dx = (lattice[ip1][j] - lattice[im1][j]) / 2.0;
    return deta_dx;
  } else {
    deta_dz = (lattice[i][jp1] - lattice[i][jm1]) / 2.0;
    return deta_dz;
  }
}

double compute_square_gradient(int i, int j) {
  double grad_sq = grad_x[i][j]*grad_x[i][j] + grad_z[i][j]*grad_z[i][j];
  return grad_sq;
}

void calc_gradient_matrix() {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      grad_x[i][j] = compute_gradient(i, j, 'x');
      grad_z[i][j] = compute_gradient(i, j, 'z');
    }
  }
  return;
}

void update_gradient_matrix(int i, int j) {
  // Implement PBCs and update graidents at neighbours of site (i,j)
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

  grad_x[ip1][j] = compute_gradient(ip1, j, 'x');
  grad_x[im1][j] = compute_gradient(im1, j, 'x');
  grad_z[i][jp1] = compute_gradient(i, jp1, 'z');
  grad_z[i][jm1] = compute_gradient(i, jm1, 'z');
  return;
}

double self_energy(double eta) {
  double se = 0.0;
  double eta2 = eta * eta;
  double eta3 = eta * eta * eta;
  double eta4 = eta2 * eta2;
  // test results without gradient term and no temperature dependence
//   se = a * eta2 / 2.0 + c * eta4 / 4.0;

//   // test results without quartic term and no temperature dependence
//   se = a * eta2 / 2.0;

//   // full landau-ginzburg energy for ising model
//   se = a * (temp - T0) * eta2 / 2.0 + c * eta4 / 4.0;
 
  // full landau-ginzburg energy for ligands
  se = h * (temp - T0) * eta + a * eta2 / 2.0 + c * eta4 / 4.0;
//   se = a * eta2 / 2.0 + d0 * (temp - T0) * eta3 / 3.0 + c * eta4 / 4.0;

  return se; 
}

double interaction_energy(int i, int j) {
//   // test results without gradient term
//   double ie = 0.0;

  // test results with gradient term
  double ie = kappa * compute_square_gradient(i, j) / 2.0;

  return ie;
}

double mc_sweep(gsl_rng *r) {
  double result = 0.0;
  int count = 1;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      result += metropolis_step(r);
      if ( fabs(E_tot - compute_total_energy()) > 1e-6 ) {
        std::cout << "ERROR!! Energies mismatched!\n";
        std::cout<< "attempt at count " << count << "\n";
        std::cout << "E_tot = " << E_tot << " and compute = " << compute_total_energy() << "\n";
        exit(1);
      }
      count++;
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

double compute_eta2() {
  double eta2 = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      eta2 += ( lattice[i][j] * lattice[i][j] );
    }
  }
  return eta2/(N*N);
}

int metropolis_step(gsl_rng *r) {
  int i, j;
  // choose a random spin to flip
  i = gsl_rng_uniform_int(r, N);
  j = gsl_rng_uniform_int(r, N);
  double eta_diff = (gsl_rng_uniform(r) - 0.5) * deta_max;	// generate random number between -0.5 and 0.5; scale by deta_max

  // calculate  energy change for spin flip
//   eta_diff = 1;
  double dE = compute_dE(i, j, eta_diff);
  double dummy1, dummy2;

  // reject every move that changes boundaries
  if (bound_on == 1 && (j == 0 || j == 1 || j == N-2 || j == N-1)) {
    return 0;
  } 
  // try to change the value of eta
  else {
    // flip spin if it's favourable
    if (dE < 0.0) {
      dummy1 = compute_total_energy();
      lattice[i][j] += eta_diff;
      E_tot += dE;
      update_gradient_matrix(i, j);
      dummy2 = compute_total_energy();
  //     std::cout << "accept dE<0: i = " << i << " j = " << j << " dE = " << dE <<  " and dummy2-dummy1 is = " << dummy2 - dummy1 << "\n";
      return 1;
    }
    else {
//       return 0;
      double rand = gsl_rng_uniform(r);
      if (exp(-dE*beta) > rand) {
        dummy1 = compute_total_energy();
        lattice[i][j] += eta_diff;
        E_tot += dE;
        update_gradient_matrix(i, j);
        dummy2 = compute_total_energy();
  //       std::cout << "accept dE>0: i = " << i << " j = " << j << " dE = " << dE <<  " and dummy2-dummy1 is = " << dummy2 - dummy1 << "\n";
        return 1;
      }
      else {
  //       std::cout << "rejected\n";
        return 0;
      }
    }
  }
}

double compute_dE(int i, int j, double deta) {
  double dE = 0.0;
  double d_self_en = 0.0;
  d_self_en = self_energy(lattice[i][j] + deta) - self_energy(lattice[i][j]);
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

  d_int_en += kappa * ( -deta * grad_x[ip1][j] + deta*deta / 4 ) / 2;
  d_int_en += kappa * ( +deta * grad_x[im1][j] + deta*deta / 4 ) / 2;
  d_int_en += kappa * ( -deta * grad_z[i][jp1] + deta*deta / 4 ) / 2;
  d_int_en += kappa * ( +deta * grad_z[i][jm1] + deta*deta / 4 ) / 2;
//   std::cout << "d_int_en: " << d_int_en << "\n";

  dE = d_self_en + d_int_en;
  return dE;
}

// find total energy of the system
double compute_total_energy() {
  double total_energy = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      total_energy += self_energy(lattice[i][j]) + interaction_energy(i, j);
    }
  }
  return total_energy;
}

// return self energy of the system
double get_self_energy() {
  double val_self_energy = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      val_self_energy += self_energy(lattice[i][j]);
    }
  }
  return val_self_energy;
}

// return interaction energy of the system
double get_interaction_energy() {
  double val_int_energy = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      val_int_energy += interaction_energy(i, j);
//       std::cout << i << " " << j << " " << interaction_energy(i, j) << "\n";
    }
  }
  return val_int_energy;
}

