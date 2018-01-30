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
//-- Integration point convection operators --//
void ipoint(int eqn, int e, double dfs, double xip[npe1][3],
 double ic[nelm][3][3][npe1][2*npe1])
{
  int i, j, k, m, n;
  int iip1, iip2, np1, np2, np3, np4;
  double rr, ai, bi, x1, y1, x2, y2, xm, ym, piv;
  double xc, yc, sm, mflow, vmag, xi, yi, p, flip;
  double l, u1, u2, v1, v2, lt2, ld2, dx, dy, msum;
  double pm, dqdx, dqdy, jac, a, b, m1, m2, fct[npe1];
  double sl[nelm][npe1][4], cc[npe1][2*npe1];

  // Re-initialize condensing coefficients
  for (i = 1; i <= nnpe; ++i) {
    fct[i] = 0.0;
    for (n = 1; n <= 2*nnpe; ++n) {
      cc[i][n] = 0.0;
      for (k = 1; k <= 2; ++k) {
        for (m = 1; m <= 2; ++m) {
          ic[e][m][k][i][n] = 0.0;
        }
      }
    }
  }

  // Streamwise direction
  for (j = 1; j <= nnpe; ++j) {
    n = ie[e][j];
    u1 = 0.0;
    v1 = 0.0;
    u2 = 0.0;
    v2 = 0.0;
    s[1][1] = lc[j][1];
    s[1][2] = lc[j][2];
    rr = ((1.0 -fl[n][1])*prs +fl[n][1]*prl)/(prl +1.0E-16);
    shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);

    // Droplet flow velocities
    if (eqn == 4) {
      m = o[j+2];
      u1 = vi[e][j][3];
      v1 = vi[e][j][4];
      u2 = vi[e][m][3];
      v2 = vi[e][m][4];
    }

    // Main flow velocities
    else {
      for (i = 1; i <= nnpe; ++i) {
        n = ie[e][i];
        m = o[j+2];
        u1 = u1 +nn[i][j]*vps[n][6];
        v1 = v1 +nn[i][j]*vps[n][7];
        u2 = u2 +nn[i][m]*vps[n][6];
        v2 = v2 +nn[i][m]*vps[n][7];
      }
    }
    vmag = sqrt(u1*u1 +v1*v1);
    ai = u1/(vmag +1.0E-16);
    bi = v1/(vmag +1.0E-16);

    // Sub-control-volume coordinates
    mflow = u1*ds[j][2] -v1*ds[j][1];
    if (mflow < 0.0) {
      np1 = o[j+1];
      np2 = o[j+2];
      np3 = o[j+3];
      np4 = o[j];
      iip1 = o[j+2];
      iip2 = o[j];
    }
    else if (mflow >= 0.0) {
      np1 = o[j+2];
      np2 = o[j+1];
      np3 = o[j];
      np4 = o[j+3];
      iip1 = o[j];
      iip2 = o[j+2];
    }
    x1 = x[ie[e][np1]][1];
    y1 = x[ie[e][np1]][2];
    x2 = x[ie[e][np2]][1];
    y2 = x[ie[e][np2]][2];
    xm = 0.5*(x2 +x[ie[e][np3]][1]);
    ym = 0.5*(y2 +x[ie[e][np3]][2]);
    xc = 0.25*(x1 +x2 +x[ie[e][o[j+3]]][1] +x[ie[e][o[j]]][1]);
    yc = 0.25*(y1 +y2 +x[ie[e][o[j+3]]][2] +x[ie[e][o[j]]][2]);

    // Equate lines at streamline / element intersection
    p = (bi*(xip[j][1] -xm) -ai*(xip[j][2] -ym))/((x2 -xm)*bi
     -(y2 -ym)*ai +1.0E-16);

    // Adjacent external surface upwinding
    if (p > 1.0) {
      p = (bi*(xip[j][1] -x2) -ai*(xip[j][2] -y2))/(((x1 -x2)
       *bi -(y1 -y2)*ai) +1.0E-16);
      xi = x2 +p*(x1 -x2);
      yi = y2 +p*(y1 -y2);
      b = sqrt((x1 -x2)*(x1 -x2) +(y1 -y2)*(y1 -y2));
      a = b*(1.0 -fabs(p));
      l = sqrt((xip[j][1] -xi)*(xip[j][1] -xi) +(xip[j][2]
       -yi)*(xip[j][2] -yi));
      ic[e][1][1][j][np2] = (a/b)*rr*vmag/(l +1.0E-16);
      ic[e][1][1][j][np1] = (1.0 -a/b)*rr*vmag/(l +1.0E-16);
      ic[e][2][1][j][np2] = (a/b)*rr*vmag/(l +1.0E-16);
      ic[e][2][1][j][np1] = (1.0 -a/b)*rr*vmag/(l +1.0E-16);
    }

    // Opposite external surface upwinding
    else if ((p >= 0.0) && (p <= 1.0)) {
      xi = xm +p*(x2 -xm);
      yi = ym +p*(y2 -ym);
      b = 2.0*sqrt((x2 -xm)*(x2 -xm) +(y2 -ym)*(y2 -ym));
      a = 0.5*b*(1.0 +fabs(p));
      l = sqrt(pow(xip[j][1] -xi, 2.0) +pow(xip[j][2] -yi, 2.0));
      ic[e][1][1][j][np2] = (a/b)*rr*vmag/(l +1.0E-16);
      ic[e][1][1][j][np3] = (1.0 -a/b)*rr*vmag/(l +1.0E-16);
      ic[e][2][1][j][np2] = (a/b)*rr*vmag/(l +1.0E-16);
      ic[e][2][1][j][np3] = (1.0 -a/b)*rr*vmag/(l +1.0E-16);
    }

    // Adjacent interior surface upwinding
    else {
      p = (bi*(xip[j][1] -xm) -ai*(xip[j][2] -ym))/(((xc -xm)*bi
       -(yc -ym)*ai) +1.0E-16);

      // Outer interior section upwinding
      if (p <= 0.5) {
        xi = xm +p*(xc -xm);
        yi = ym +p*(yc -ym);
        b = 0.5*sqrt((xc -xm)*(xc -xm) +(yc -ym)*(yc -ym));
        a = b*(1.0 -2.0*fabs(p));
        l = sqrt(pow(xip[j][1] -xi, 2.0) +pow(xip[j][2]
         -yi, 2.0));
        ic[e][1][1][j][np3] = 0.5*(a/b)*rr*vmag/(l +1.0E-16);
        ic[e][1][1][j][np2] = 0.5*(a/b)*rr*vmag/(l +1.0E-16);
        ic[e][2][1][j][np3] = 0.5*(a/b)*rr*vmag/(l +1.0E-16);
        ic[e][2][1][j][np2] = 0.5*(a/b)*rr*vmag/(l +1.0E-16);
        cc[j][iip1] = -(1.0 -a/b)*rr*vmag/(l +1.0E-16);
      }

      // Inner interior section upwinding
      else {
        xi = xm +p*(xc-xm);
        yi = ym +p*(yc-ym);
        b = sqrt(pow(xip[iip1][1] -xip[iip2][1], 2.0)
         +pow(xip[iip1][2] -xip[iip2][2], 2.0) +1.0E-16);
        a = sqrt(pow(xip[iip2][1] -xi, 2.0) +pow(xip[iip2][2]
         -yi, 2.0));
        l = sqrt(pow(xip[j][1] -xi, 2.0) +pow(xip[j][2]
         -yi, 2.0));
        cc[j][iip1] = -(a/b)*rr*vmag/(l +1.0E-16);
        cc[j][iip2] = -(1.0 -a/b)*rr*vmag/(l +1.0E-16);
      }
    }

    // Multiphase region permeability
    m = ie[e][o[j]];
    flip = 2.0*fl[n][1]*fl[m][1]/(fl[n][1] +fl[m][1] +1.0E-16);
    pm = 1.0/(1.0 +k0*(1.0 -flip)*(1.0 -flip)/(flip*flip*flip
     +1.E-16));
    if ((eqn <= 2) || (eqn == 4)) {pm = 1.0;}

    // Mass weighted (positive coefficient) upwinding
    lt2 = ds[j][1]*ds[j][1] +ds[j][2]*ds[j][2];
    ld2 = 0.5*ds[j][3]*ds[j][3]/lt2 +0.375*lt2;
    m1 = u1*ds[j][2] -v1*ds[j][1];
    m2 = u2*ds[o[j+2]][2] -v2*ds[o[j+2]][1];
    if (m1 >= 0.0) {sm = fmax(fmin(-m2/(m1 +1.0E-16), 1.0), 0.0);}
    if (m1 < 0.0) {sm = fmax(fmin(m2/(m1 +1.0E-16), 1.0), 0.0);}
    if (eqn <= 1) {
      for (n = 1; n <= 2*nnpe; ++n) {
        cc[j][n] = 0.0;
        for (k = 1; k <= 2; ++k) {
          for (m = 1; m <= 2; ++m) {
            ic[e][m][k][j][n] = 0.0;
          }
        }
      }
      cc[j][j] = pm*rr*vmag/(l +1.0E-16) +dfs/(ld2 +1.0E-16);
      cc[j][iip1] = -sm*rr*vmag/l;
      cc[j][4+j] = 1.0;
      ic[e][1][1][j][np2] = (1.0 -sm)*rr*vmag/(l +1.0E-16);
      ic[e][2][1][j][np2] = (1.0 -sm)*rr*vmag/(l +1.0E-16);
    }

    // Skew upwind convection operator
    if (eqn > 1) {
      cc[j][j] = pm*rr*vmag/(l +1.0E-16) +dfs/(ld2 +1.0E-16);
      cc[j][4+j] = 1.0;
      sl[e][j][1] = sm;
      sl[e][j][2] = l;
      sl[e][j][3] = vmag;

      // Nodal velocity influence coefficients
      for (i = 1; i <= nnpe; ++i) {
        ic[e][1][1][j][i] = pm*ic[e][1][1][j][i] +dfs*nn[i][j]
         /(ld2 +1.0E-16);
        ic[e][2][1][j][i] = pm*ic[e][2][1][j][i] +dfs*nn[i][j]
         /(ld2 +1.0E-16);
      }
    }

    // Nodal pressure influence coefficients
    for (i = 1; i <= nnpe; ++i) {
      ic[e][1][2][j][i] = -pm*dn[i][3];
      ic[e][2][2][j][i] = -pm*dn[i][4];
    }
  }

  // Integration point -> nodal dependence (matrix inversion)
  for (i = 1; i <= nnpe; ++i) {
    piv = cc[i][i];
    for (m = 1; m <= 2*nnpe; ++m) {
      cc[i][m] = cc[i][m]/piv;
      for (k = 1; k <= nnpe; ++k) {
        if ((k != i) && (m == i)) {fct[k] = cc[k][i]/cc[i][i];}
        if (k != i) {cc[k][m] = cc[k][m] -fct[k]*cc[i][m];}
      }
    }
  }

  // Isolate integration point variable dependence
  for (m = 1; m <= 2; ++m) {
    for (n = 1; n <= 2; ++n) {
      for (i = 1; i <= nnpe; ++i) {
        for (k = 1; k <= nnpe; ++k) {
          msum = 0.0;
          for (j = 1; j <= nnpe; ++j) {
            msum = msum +cc[i][4+j]*ic[e][m][n][j][k];
          }
          ic[e][m][n][i][4+k] = msum;
        }
      }
    }
  }
}																									// calls no functions 