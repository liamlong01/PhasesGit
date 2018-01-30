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
//-- Boundary conditions --//
void bndry(int v, int p, int m, int nsrf, double fo,
 double bc[2*nsrfm][5][4], double c[3*nnpm][6*nym+11],
 double r[3*nnpm])
{
  int e, f, i, j, n, np, s, n1, n2, ni, p1, p2, p3;
  int n3, mbc, b1, b2, bi, isemi, iband, oi;
  double dx, dy, ci, ti, tm, tp, bval;

  // Boundary coefficients (customized)
  bval = 0.0;
  tm = fo*lr*lr*prl*phl/pkl;
  for (j = 1; j <= csb[1][1]; ++j) {
    i = int(csb[j+1][1]);
    ci = cn[bgn[i][1]][1];
    ti = tmin +tn[bgn[i][1]][1]*fabs(tmax -tmin);
    n = int(csb[j+1][2]) -10;
    f = int(csb[j+1][3]) -30;
    b1 = int(2.0*(i -1) +1);
    b2 = b1 +1;
    if ((f >= 1) && (f <= 3)) {
      custom(3, j, tm, ci, ti, csp, css, csb, bval);													// calls custom
      bc[b1][n][f] = bval;
      bc[b2][n][f] = bval;
      if (f == 1) {
        bc[b1][2][1] = bval/lr;
        bc[b2][2][1] = bval/lr;
      }
      else if (f == 3) {
        bc[b1][2][3] = (bval -bc[b1][2][2]*tmin)/fabs(tmax -tmin);
        bc[b2][2][3] = (bval -bc[b2][2][2]*tmin)/fabs(tmax -tmin);
      }
    }
  }

  // Loop over boundary surfaces
  iband = m*(band3 +m -3)/3;
  isemi = (iband +1)/2;
  p1 = 0;
  p2 = 0;
  p3 = 0;
  for (s = 1; s <= nsrf; ++s) {
    e = bel[s];
    b1 = int(2.0*(s -1) +1);
    b2 = b1 +1;
    n1 = p +m*(bgn[s][1] -1);
    n2 = p +m*(bgn[s][2] -1);
    dx = 0.5*(x[bgn[s][2]][1] -x[bgn[s][1]][1]);
    dy = 0.5*(x[bgn[s][2]][2] -x[bgn[s][1]][2]);

    // Loop over element nodes
    for (i = 1; i <= 2; ++i) {
      if (i == 1) {
        ni = n1;
        bi = b1;
      }
      else {
        ni = n2;
        bi = b2;
      }

      // Continuity equation
      if (p == 3) {
        mbc = 1; // mbc = 1: gradient spacified condition
        // mbc = 2: pressure specified at single boudnary
        // mbc = 3: pressure specified at multiple boundaries

        // Boundary mass equation closure
        if (mbc == 1) {
          n1 = 3 +m*(bgn[s][1] -1);
          n2 = 3 +m*(bgn[s][2] -1);
          n3 = m*(bgn[s][2] -bgn[s][1]);

          // Gradient specified condition
          if ((bc[b1][3][2] == 0.0) || (bc[b1][4][2] == 0.0)) {
            c[n1][isemi-2] = c[n1][isemi-2] -dy*0.75;
            c[n1][isemi-2+n3] = c[n1][isemi-2+n3] -dy*0.25;
            c[n1][isemi-1] = c[n1][isemi-1] +dx*0.75;
            c[n1][isemi-1+n3] = c[n1][isemi-1+n3] +dx*0.25;
            c[n2][isemi-2] = c[n2][isemi-2] -dy*0.75;
            c[n2][isemi-2-n3] = c[n2][isemi-2-n3] -dy*0.25;
            c[n2][isemi-1] = c[n2][isemi-1] +dx*0.75;
            c[n2][isemi-1-n3] = c[n2][isemi-1-n3] +dx*0.25;
          }
 
          // Velocity specified condition
          else {
            r[n1] = r[n1] +dy*bc[b1][3][3]/(bc[b1][3][2]
             +1.0E-16) -dx*bc[b1][4][3]/(bc[b1][4][2] +1.0E-16);
            r[n2] = r[n2] +dy*bc[b2][3][3]/(bc[b2][3][2]
             +1.0E-16) -dx*bc[b2][4][3]/(bc[b2][4][2] +1.0E-16);
          }
        }

        // Pressure based conditions
        else if (mbc >= 2) {
          np = 0;
          for (f = 3; f <= 4; ++f) {
            for (j = 0; j <= 1; ++j) {
              bi = int(2.0*(s -1) +1) +j;
              if ((bc[bi][f][1] == 0) && (bc[bi][f][3] != 0))
                {np = 1;}
              if ((bc[bi][f][2] == 0) && (bc[bi][f][3] == 0))
                {np = 2;}
            }
          }

          // Reference pressure specification
          if (np == 1) {p1 = p1++;}
          if (np != 1) {p2 = p2++;}
          if ((p1 > 1) && (p2 > 1) && (np != 1)) {p3 = p3++;}
          if (((mbc == 2) && (np == 1) && (p3 <= 1))
            || ((mbc == 3) && (np == 1))) {
            for (j = 1; j <= iband; ++j) {
              c[ni][j] = 0.0;
            }
            c[ni][isemi] = 1.0;
            r[ni] = 0.0;
          }

          // Zero gradient condition
          if (np == 2) {
            oi = int((bc[bi][3][1] -0.99999)*1.0E+03);
            f = isemi +m*(ie[e][o[oi+1]] -ie[e][o[oi+2]]);
            if (f > 0) {
              for (j = 1; j <= iband; ++j) {
                c[ni][j] = 0.0;
              }
              c[ni][isemi] = 1.0;
              c[ni][f] = -1.0;
              r[ni] = 0.0;
            }
          }
        }
      }

      // Dirichlet boundary condition
      else if (bc[bi][v][1] == 0.0) {
        for (j = 1; j <= iband; ++j) {
          c[ni][j] = 0.0;
        }
        c[ni][isemi] = 1.0;
        r[ni] = bc[bi][v][3]/(bc[bi][v][2] +1.0E-32);
      }

      // Neumann boundary condition
      else if ((bc[bi][v][1] > 1.0) && (bc[bi][v][1] < 1.01)
       && (bc[bi][v][2] == 0.0)) {
        oi = int((bc[bi][v][1] -0.99999)*1.0E+03);
        for (j = 1; j <= iband; ++j) {
          c[ni][j] = 0.0;
        }
        c[ni][isemi] = 1.0;
        c[ni][isemi+m*(ie[e][o[oi+1]] -ie[e][o[oi+2]])] = -1.0;
        r[ni] = sqrt(dx*dx +dy*dy)*bc[bi][v][3]/(bc[bi][v][1]
         +1.0E-32);
      }
    
      // Robin boundary condition
      else if ((bc[bi][v][1] != 0) && (bc[bi][v][2] != 0)) {
        tp = r[ni]/(c[ni][isemi] + 1.0E-32);
        for (j = 1; j <= iband; ++j) {
          if (j != isemi) {
            tp = tp-c[ni][j]*tn[j][1]/(c[ni][isemi]+1.0E-32);
          }
        }
        bc[bi][v][3] = bc[bi][v][3]/(bc[bi][v][1] +1.0E-32);
        bc[bi][v][2] = bc[bi][v][2]/(bc[bi][v][1] +1.0E-32);
        bc[bi][v][1] = 1.0;
        r[ni] = r[ni] +sqrt(dx*dx +dy*dy)*(bc[bi][v][3]
         -bc[bi][v][2]*tp);
      }
    }
  }

  // Common edge of doubly connected region
  for (i = 1; i <= nym; ++i) {
    if (dbc[i][1] > 0) {
      for (n = 1; n <= 4; ++n) {
        e = dbc[i][1];
        if (n <= 2) {s = dbc[i][n+1];}
        if (n >= 3) {s = dbc[i][n+2];}
        ni = p +m*(s -1);
        for (j = 1; j <= iband; ++j) {
          if (j != isemi) {c[ni][j] = 0.0;}
        }
        c[ni][isemi] = 1.0;
        if (v == 1) {r[ni] = cn[s][1];}
        if (v == 2) {r[ni] = tn[s][1];}
        if (v == 4) {r[ni] = vps[s][1];}
        if (v == 5) {r[ni] = vps[s][2];}
      }
    }
  }
}