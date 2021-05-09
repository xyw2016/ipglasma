#include "mpi.h"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

//#include "Evolution.h"
//#include "FFT.h"
//#include "Init.h"
//#include "Lattice.h"
//#include "Matrix.h"
//#include "Parameters.h"
//#include "Random.h"
//#include "Setup.h"
//#include "Spinor.h"
//#include "pretty_ostream.h"
#include "IPGlasma.h"

#define _SECURE_SCL 0
#define _HAS_ITERATOR_DEBUGGING 0
using namespace std;

void display_logo();
void writeparams(Parameters *param);

// main program 1
int main(int argc, char *argv[]) {
    int rank;
    int size;

    int nev = 1;
    if (argc == 3) {
        nev = atoi(argv[2]);
    }

    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes

    IPGlasma IPGgenerator(rank, size, nev, "input");

    // event loop starts ...
    for (int iev = 0; iev < nev; iev++) {
        IPGgenerator.generateAnEvent(iev);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 1;
}


void display_logo() {
  cout << endl;
  cout << "--------------------------------------------------------------------"
          "---------"
       << endl;
  cout << "| Classical Yang-Mills evolution with IP-Glasma initial "
          "configurations v1.4 |"
       << endl;
  cout << "--------------------------------------------------------------------"
          "---------"
       << endl;
  cout << "| References:                                                       "
          "        |"
       << endl;
  cout << "| B. Schenke, P. Tribedy, R. Venugopalan                            "
          "        |"
       << endl;
  cout << "| Phys. Rev. Lett. 108, 252301 (2012) and Phys. Rev. C86, 034908 "
          "(2012)     |"
       << endl;
  cout << "--------------------------------------------------------------------"
          "---------"
       << endl;

  cout << "This version uses Qs as obtained from IP-Sat using the sum over "
          "proton T_p(b)"
       << endl;
  cout << "This is a simple MPI version that runs many events in one job. No "
          "communication."
       << endl;

  cout << "Run using large lattices to improve convergence of the root finder "
          "in initial condition. "
       << "Recommended: 600x600 using L=30fm" << endl;
  cout << endl;
}


void writeparams(Parameters *param)
{
  // write the used parameters into file "usedParameters.dat" as a double check
  // for later
  time_t rawtime = time(0);
  stringstream strup_name;
  strup_name << "usedParameters" << param->getEventId() << ".dat";
  string up_name;
  up_name = strup_name.str();

  fstream fout1(up_name.c_str(), ios::out);
  char *timestring = ctime(&rawtime);
  fout1 << "File created on " << timestring << endl;
  fout1 << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
  fout1 << "Used parameters by IP-Glasma v1.3" << endl;
  fout1 << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
  fout1 << " " << endl;
  fout1 << " Output by readInput in main.cpp: " << endl;
  fout1 << " " << endl;
  fout1 << "Program run in mode " << param->getMode() << endl;
  fout1 << "Nc " << param->getNc() << endl;
  fout1 << "size " << param->getSize() << endl;
  fout1 << "lattice spacing a "
        << param->getL() / static_cast<double>(param->getSize()) << " fm "
        << endl;
  fout1 << "Ny " << param->getNy() << endl;
  fout1 << "Projectile " << param->getProjectile() << endl;
  fout1 << "Target " << param->getTarget() << endl;
  if (param->getUseConstituentQuarkProton() > 0) {
    fout1 << "Nucleons consists of " << param->getUseConstituentQuarkProton()
          << " constituent quarks" << endl;
    if (param->getShiftConstituentQuarkProtonOrigin())
      fout1 << "... constituent quark center of mass moved to origin" << endl;
  }
  fout1 << "Smooth nucleus " << param->getUseSmoothNucleus() << endl;
  fout1 << "Gaussian wounding " << param->getGaussianWounding() << endl;
  fout1 << "Using fluctuating x=Qs/root(s) " << param->getUseFluctuatingx()
        << endl;
  if (param->getRunWithkt() == 0)
    fout1 << "Using local Qs to run " << param->getRunWithLocalQs() << endl;
  else
    fout1 << "running alpha_s with k_T" << endl;
  fout1 << "QsmuRatio " << param->getQsmuRatio() << endl;
  fout1 << "smeared mu " << param->getSmearQs() << endl;
  fout1 << "m " << param->getm() << endl;
  fout1 << "rmax " << param->getRmax() << endl;
  fout1 << "UVdamp " << param->getUVdamp() << endl;
  fout1 << "beta2 " << param->getbeta2() << endl;
  if (param->getSmearQs() == 1) {
    fout1 << "smearing width " << param->getSmearingWidth() << endl;
  }
  fout1 << "Using fat tailed distribution " << param->getUseFatTails() << endl;
  fout1.close();
}
