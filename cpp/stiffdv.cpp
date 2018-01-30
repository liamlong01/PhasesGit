// -------------------------------------------------
// C++ PHASES - Copyright by G. F. Naterer
// Coding format: (i) // comments, (ii) column < 65,
// (ii) line continued -> one column indentation,
// (iv) start block -> two column indentation,
// (v) code optimization: set maximum speed options
// -------------------------------------------------

#include <iostream> // standard stream I/O facilities
#include <fstream> // I/O facilities for external files
#include <math.h> // standard mathematical functions
#define WIN32_LEAN_AND_MEAN // exclude rarely used headers
using namespace std;

#include "constInit.h"  //Initializes all constants
#include "variaExtern.h"  //Declares externally all global variables
#include "functInit.h"  //Declares all functions

///////////////////////////////////////////////////////////////////////////////
//-- Element momentum equations: droplet trajectories --//
void stiffdv(int e, double vi[nelm][npe1][7])
{
  int n, i, j, k;
  double du, dv, dd, dx, dy, rf, rd, cd;
  double uc, vc, jac, dudx, dudy, dvdx, dvdy;

  dd = 0.001*pcs +1.0E-32;
  ur = vr*pkl/(prl*phl*lr +1.0E-32);
  for (i = 1; i <= nnpe; ++i) {
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    s[2][1] = lc[4+i][1];
    s[2][2] = lc[4+i][2];
    s[3][1] = lc[8+i][1];
    s[3][2] = lc[8+i][2];
    for (j = 1; j <= 2; ++j) {
      for (k = 1; k <= nnpe; ++k) {
        n = ie[e][k];
        ivo[k] = vps[n][3+j];
      }
      if (j == 1) {shape(e, nnpm, dx, dy, s, dn, ivo,
       jac, dudx, dudy);}																			// calls shape
      if (j == 1) {shape(e, nnpm, dx, dy, s, dn, ivo,
       jac, dvdx, dvdy);}																			// calls shape
    }

    // Droplet velocities: momentum equations
    du = ur*(vi[e][i][3] -vi[e][i][1]);
    dv = ur*(vi[e][i][4] -vi[e][i][2]);
    rd = (fabs(du) +fabs(dv))*dd/(vsc +1.0E-32) +1.0E-32;
    cd = 24.0/rd +4.73/pow(rd, 0.37) +0.00624*pow(rd, 0.38);
    uc = ur*ur*(vi[e][i][3]*dudx +vi[e][i][4]*dudy)
     /(lr +1.0E-32);
    vc = -9.8 +ur*ur*(vi[e][i][3]*dvdx +vi[e][i][4]*dvdy)
     /(lr +1E-32);
    rf = (dd*dd*24/rd)/(18*vsc*cd +1.0E-32);
    vi[e][i][3] = vi[e][i][5] +uc*rf/(ur +1.0E-32);
    vi[e][i][4] = vi[e][i][6] +vc*rf/(ur +1.0E-32);
  }
}