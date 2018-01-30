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
// Stiffness matrices: C equation (continuity / droplet flow)
void stiffdc(int e, double df, double xip[npe1][3],
 double flux[4][npe1][npe1], double aq[4][npe1][npe1],
 double rq[npe1])
{
  int m, n, i, j;
  double advc, mflow, mw;

  // Convection coefficients
  ipoint(4, e, 0.0, xip, ic);																			// calls ipoint
  for (i = 1; i <= nnpe; ++i) {
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    m = ie[e][i];
    mflow = vi[e][i][3]*ds[i][2] -vi[e][i][4]*ds[i][1];

    // Mass flux of droplets
    for (j = 1; j <= nnpe; ++j) {
      n = ie[e][j];
      mw = 2.0 -2.0/(1.0 +pow(fl[n][1], 0.1));
      if (ph[n][1] == 2) {mw = 0.0;}
      advc = -ic[e][1][1][i][4+j]*mflow*mw;
      flux[1][i][j] = vr*advc;
    }
  }

  // Local stiffness matrix: phase volume fraction
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    for (j = 1; j <= nnpe; ++j) {
      if (i == j) {aq[1][i][j] = vr*ds[i][3]/df -flux[1][i][j]
       +flux[1][o[i]][j];}
      if (i != j) {aq[1][i][j] = -flux[1][i][j]
       +flux[1][o[i]][j];}
    }
    rq[i] = vr*cn[n][2]*ds[i][3]/df;
  }
}