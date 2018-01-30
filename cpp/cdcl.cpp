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
//-- Drag and lift coefficients by pressure and viscous forces --//
void cdcl(int nnp, int nsrf, double &um, double &cd, double &cl)
{
  int e, i, j, k, m, n, p, r, n1, n2;
  double sa, dudx, dudy, dvdx, dvdy, dsx, dsy;
  double jac, dx, dy, u1, drag, lift, dels, pw;

  // Reference or maximum velocity
  um = 0.0;
  sa = 0.0;
  ur = vr*pkl/(prl*phl*lr +1.0E-32);
  for (i = 1; i <= nsrf; ++i) {
    for (k = 1; k <= 2; ++k) {
      j = 2*(i -1) +k;
      m = bgn[i][1];
      u1 = sqrt(pow(bc[j][3][3]/(bc[j][3][2] +1.0E-32), 2.0)
       +pow(bc[j][4][3]/(bc[j][4][2] +1.0E-32), 2.0));
      if ((bc[j][3][1] == 0) && (bc[j][4][1] == 0)
       && (u1 > um)) {um = u1;}
    }
  }
  if (um == 0.0) {
    for (i = 1; i <= nnp; ++i) {
      u1 = sqrt(vps[i][1]*vps[i][1] +vps[i][2]*vps[i][2]);
      if (u1 > um) {um = u1;}
    }
  }
  um = um*ur;

  // Boundary surface geometry
  drag = 0.0;
  lift = 0.0;
  for (i = 1; i <= nsrf; ++i) {
    n1 = bgn[i][1];
    n2 = bgn[i][2];
    e = bel[i];
    dels = sqrt(pow(x[n2][1] -x[n1][1], 2.0) +pow(x[n2][2]
     -x[n1][2], 2.0));
    dsx = (x[n2][1] -x[n1][1])/(dels +1.0E-32);
    dsy = (x[n2][2] -x[n1][2])/(dels +1.0E-32);
    for (k = 1; k <= 2; ++k) {
      j = 2*(i -1) +k;
      p = bgn[i][k];
      pw = vps[p][3];

      // Drag and lift coefficients on no-slip boundaries
      if ((bc[j][3][1] == 0) && (bc[j][3][3] == 0)
       && (bc[j][4][1] == 0) && (bc[j][4][3] == 0)) {
        for (r = 1; r <= 2; ++r) {
          for (n = 1; n <= nnpe; ++n) {
            m = ie[e][n];
            if (m == p) {
              s[1][1] = nc[8+2*(n-1)+1][1];
              s[1][2] = nc[8+2*(n-1)+1][2];
            }
            ivo[n] = vps[m][r];
          }
          if (r == 1) {shape(e, nnpm, dx, dy, s, dn, ivo,															// calls shape
           jac, dudx, dudy);}
          if (r == 2) {shape(e, nnpm, dx, dy, s, dn, ivo,															// calls shape
           jac, dvdx, dvdy);}
        }
        drag = drag +(0.5*lr*dels)*(dsy*prl*pow(um, 2.0)*pw
         +dsx*prl*vsc*um*(dudy +dvdx)/(lr +1.0E-32));
        lift = lift +(0.5*lr*dels)*(-dsx*prl*pow(um, 2.0)*pw
         +dsy*prl*vsc*um*(dudy +dvdx)/(lr +1.0E-32));
        sa = sa +0.5*dels;
      }
    }
  }
  cd = fabs(drag/(0.5*prl*um*um*lr*sa +1.0E-32));
  cl = lift/(0.5*prl*um*um*lr*sa +1.0E-32);
}