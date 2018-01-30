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
//-- Smoothing of entropy based diffusivity --//
void dsmooth(int nel, int nnp, double vps[nnpm][14])
{
  int e, i, j, k, kn, n, m;
  double jac, dx, dy, dmdx, dmdy;
  double flux[4][npe1][npe1];

  for (n = 1; n <= nnp; ++n) {
    vps[n][12] = 0.0;
    vps[n][13] = 0.0;
  }
  kn = 5;

  // Assembly of diffusion terms
  for (k = 1; k <= kn; ++k) {
    for (e = 1; e <= nel; ++e) {
      for (i = 1; i <= nnpe; ++i) {
        m = ie[e][i];
        s[1][1] = lc[i][1];
        s[1][2] = lc[i][2];
        s[2][1] = lc[4+i][1];
        s[2][2] = lc[4+i][2];
        s[3][1] = lc[8+i][1];
        s[3][2] = lc[8+i][2];
        shape(e, nnpm, dx, dy, s, dn, ivo, jac, dmdx, dmdy);									// calls shape
        ds[i][1] = dx;
        ds[i][2] = dy;
        for (j = 1; j <= nnpe; ++j) {
          flux[1][i][j] = dn[j][3]*ds[i][2] -dn[j][4]*ds[i][1];
        }
      }

      // Iterative Jacobi point solver
      for (i = 1; i <= nnpe; ++i) {
        n = ie[e][i];
        for (j = 1; j <= nnpe; ++j) {
          m = ie[e][j];
          if (i == j) {
            vps[n][12] = vps[n][12] -flux[1][i][j]
             +flux[1][o[i]][j];
          }
          else {
            vps[n][13] = vps[n][13] -flux[1][i][j]*vps[m][11]
             +flux[1][o[i]][j]*vps[m][11];
          }
        }
      }
    }

    // Updated diffusivity field
    for (n = 1; n <= nnp; ++n) {
      if (vps[n][10] > 0.0) {vps[n][11] = fabs(vps[n][13]
       /(vps[n][12] +1.0E-16));}
      if (k == kn) {vps[n][11] = vps[n][11]/(pkl +1.0E-32)
       +((1.0 -fl[n][1])*pks +fl[n][1]*pkl)/(pkl +1.0E-32);}
      vps[n][12] = 0.0;
      vps[n][13] = 0.0;
    }
  }
}