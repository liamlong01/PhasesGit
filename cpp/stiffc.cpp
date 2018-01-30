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
//-- Local stiffness matrices: C equation (solid / liquid) --//
void stiffc(int e, double df, double fo, double xip[npe1][3],
 double flux[4][npe1][npe1], double sc[npe1],
 double aq[4][npe1][npe1], double rq[npe1])
{
  int i, j, k, n;
  double dfsc, advc, ci, ti, tm, le, jac, dfs, scc;
  double dx, dy, dqdx, dqdy, mflow, scadv, scdfs;
  double clm[npe1][4], rd[npe1][5], cph[npe1][3];

  le = pkl/(prl*phl*pdl0 +1.0E-32);
  dfs = 1.0/le;
  coeffc(e, cph, rd);																							// calls coeffc
  ipoint(2, e, dfs, xip, ic);																					// calls ipoint
  for (i = 1; i <= nnpe; ++i) {
    clm[i][1] = 0.0;
    clm[i][2] = 0.0;
    clm[i][3] = 0.0;
    for (n = 1; n <= nnpe; ++n) {
      clm[i][1] = clm[i][1] +nn[n][i]*cn[ie[e][n]][1];
      clm[i][2] = clm[i][2] +nn[n][i]*cph[i][2];
      clm[i][3] = clm[i][3] +nn[n][i]*fl[ie[e][n]][1];
    }
  }

  // Outer loop: local node dependence
  for (i = 1; i <= nnpe; ++i) {
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    sc[i] = 0.0;
    mflow = vi[e][i][1]*ds[i][2] -vi[e][i][2]*ds[i][1];
    shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);														// calls shape

    // Diffusion, advection, source terms
    for (j = 1; j <= nnpe; ++j) {
      n = ie[e][j];
      // dfsc = rd[i][1]*dn[j][3] -rd[i][2]*dn[j][4]
      // +rd[i][3]*dn[j][3] -rd[i][4]*dn[j][4];
      dfsc = (ds[i][2]*dn[j][3] -ds[i][1]*dn[j][4])
       *(fl[i][1]*pdl0 +(1.0 -fl[i][1])*pds0)/(pdl0 +1.0E-32);
      advc = -le*ic[e][1][1][i][4+j]*mflow;
      scdfs = (rd[i][1]*dn[j][3] -rd[i][2]*dn[j][4])
       *(cph[j][1] -cn[n][1]) -(rd[i][3]*dn[j][3]
       -rd[i][4]*dn[j][4])*(cph[j][2] -cn[n][1]);
      scadv = le*mflow*(clm[j][3]*clm[j][2] -clm[j][1]);
      sc[i] = sc[i] +vr*scadv +scdfs;
      flux[1][i][j] = vr*advc +dfsc;
    }
  }

  // Forming elemental stiffness matrix
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    for (j = 1; j <= nnpe; ++j) {
      if (i == j) {aq[1][i][j] = vr*le*ds[i][3]/df
       -flux[1][i][j] +flux[1][o[i]][j];}
      if (i != j) {aq[1][i][j] = -flux[1][i][j]
       +flux[1][o[i]][j];}
    }

    // Source terms (customized)
    scc = 0.0;
    ci = cn[n][1];
    ti = tmin +tn[n][1]*fabs(tmax-tmin);
    tm = fo*lr*lr*prl*phl/pkl;
    for (j = 1; j <= css[1][1]; ++j) {
      k = int(css[j+1][1]);
      if (k == 11) {custom(2, j, tm, ci, ti, csp, css,
       csb, scc);}																							// calls custom
    }

    // Transient and source terms (phase change + customized)
    rq[i] = vr*le*cn[n][2]*ds[i][3]/df +(sc[o[i]] -sc[i])
     +scc*ds[i][3]*(lr*lr*prl*phl/(pkl*vr));

	// cout << rq[i] << endl; 
  }
}