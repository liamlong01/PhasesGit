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
// Local stiffness matrices: u-V-p (y-momentum)
void stiffv(int e, double df, double flux[4][npe1][npe1],
 double sc[npe1], double aq[4][npe1][npe1], double rq[npe1])
{
  int i, j, n, m;
  double rr, tbar, cbar, pm, flip, dqdx, dqdy, dx, dy, pp;
  double jac, scf, advu, advv, advp, dfsu, dfsv, ppr, mflow;

  pm = (vsc*prl*phl)/(pkl +1.0E-16);
  cbar = 0.5*(pcs+ pcl);
  tbar = 0.5;

  // Start outer loop: local node dependence
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    m = ie[e][o[i]];
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    sc[i] = 0.0;
    flip = 2.0*fl[n][1]*fl[m][1]/(fl[n][1] +fl[m][1] +1.0E-16);
    rr = ((1.0 -flip)*prs +flip*prl)/(prl +1.0E-16);
    mflow = rr*vi[e][i][1]*ds[i][2] -rr*vi[e][i][2]*ds[i][1];
    shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);													// calls shape 

    // Prandtl number: molecular, turbulent + two-phase parts
    pp = k0*(1.0 -flip)*(1.0 -flip)/(flip*flip*flip+1.0E-16);
    //ppr = pr[e][i] +pp;
	ppr = mu_e;
    for (j = 1; j <= nnpe; ++j) {

      // Diffusion, advection flows
      advu = 0.0;
      dfsu = ppr*dn[j][4]*ds[i][2];
      advv = -ic[e][2][1][i][4+j]*mflow;
      dfsv = ppr*dn[j][3]*ds[i][2] -2.0*ppr*dn[j][4]*ds[i][1];
      advp = -ic[e][2][2][i][4+j]*mflow +nn[j][i]*ds[i][1];
      flux[1][i][j] = vr*advu +dfsu;
      flux[2][i][j] = vr*advv +dfsv;
      flux[3][i][j] = vr*advp;
    }
  }

  //  Forming elemental stiffness matrix
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    rr = ((1.0-fl[n][1])*prs +fl[n][1]*prl)/(prl +1.0E-16);
    for (j = 1; j <= nnpe; ++j) {
      if (i == j) {aq[2][i][j] = vr*rr*ds[i][3]/df
       -flux[2][i][j] +flux[2][o[i]][j];}
      if (i != j) {aq[2][i][j] = -flux[2][i][j]
       +flux[2][o[i]][j];}
      aq[1][i][j] = -flux[1][i][j] +flux[1][o[i]][j];
      aq[3][i][j] = -flux[3][i][j] +flux[3][o[i]][j];
    }

    // Sources (body forces and phase change), transient term
    scf = rat*pm*ds[i][3]*(tn[n][1] -tbar) +ras*pm*ds[i][3]
     *(cn[n][1] -cbar);
    rq[i] = vr*rr*vps[n][5]*ds[i][3]/df +(sc[o[i]]
     -sc[i]) +scf/vr;
  }
}