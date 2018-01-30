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
//-- Local stiffness matrices: temperature equation --//
void stifft(int e, double df, double fo, double cip[npe1],
 double tip[npe1], double flux[4][npe1][npe1], double sc[npe1],
 double aq[4][npe1][npe1], double rq[npe1])
{
  int i, j, k, n, m;
  double sct, ti, flip, kx, ky, ste, jac, rr, kkx, kky;
  double dx, dy, dqdx, dqdy, mflow, dfsq, advq, cval, hbar;
  double tphs, tphl, ec1ip, ec2ip, ec3ip, scadv, ci, tm;

  // Update customized properties
  ci = cip[1];
  kx = pks;
  ky = pks;
  ti = tmin +tip[1]*fabs(tmax-tmin);
  tm = fo*lr*lr*prl*phl/pkl;
  for (j = 1; j <= csp[1][1]; ++j) {
    k = int(csp[j+1][2]);
    custom(1, j, tm, ci, ti, csp, css, csb, cval);													// calls custom
    if (csp[j+1][1] == 21) {
      if (k == 2) {kx = cval;}
      if (k == 3) {ky = cval;}
      if ((k == 1) || (k > 3)) {
        kx = cval;
        ky = cval;
      }
    }
    if (csp[j+1][1] == 22) {phs = cval;}
  }

  ste = phl*fabs(tmax-tmin)/(pl + 1.0E-32);
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    m = ie[e][o[i]];
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    sc[i] = 0.0;
    shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);											// calls shape 

    // Integration point reference states
    tphs = fmax(tmlt -cip[i]*(tmlt -eps -tsol)/pcs -eps, tsol);
    tphl = tmlt -cip[i]*(tmlt +eps -tsol)/(pcl +1.0E-32) +eps;
    hbar = 0.5*(phs+phl);

    // Solid phase reference states
    if (tip[i] <= tphs) {
      ec1ip = 0.0;
      ec2ip = phs/(phl +1.0E-32);
      ec3ip = 0.0;
    }

    // Liquid phase reference states
    else if (tip[i] >= tphl) {
      ec1ip = hbar*(tphl-tphs)/(phl +1.0E-32) 
       +(phs*tphs)/(phl +1.0E-32) +1.0/(ste +1.0E-32);
      ec2ip = phl/(phl +1.0E-32);
      ec3ip = tphl;
    }

    // Melt phase reference states
    else {
      ec1ip = (phs*tphs)/(phl + 1.0E-32);
      ec2ip = hbar/(phl +1.0E-32) +1.0/(ste*(tphl -tphs)
       +1.0E-32);
      ec3ip = tphs;
    }
    flip = 2.0*fl[n][1]*fl[m][1]/(fl[n][1] +fl[m][1] +1.0E-32);
    rr = ((1.0-flip)*prs +flip*prl)/(prl +1.0E-32);
    kkx = ((1.0-flip)*kx +flip*pkl)/(pkl +1.0E-32);
    kky = ((1.0-flip)*ky +flip*pkl)/(pkl +1.0E-32);

    // Local entropy based corrections
    if (ao == 1) {kkx = vps[n][11];}
    if (ao == 1) {kky = vps[n][11];}
    mflow = rr*vi[e][i][1]*ds[i][2] -rr*vi[e][i][2]*ds[i][1];

    // Diffusive, advective and source terms
    for (j = 1; j <= nnpe; ++j) {
      n = ie[e][j];
      dfsq = kkx*ds[i][2]*dn[j][3] -kky*ds[i][1]*dn[j][4];
      advq = -ic[e][1][1][i][4+j]*ec2ip*mflow;
      flux[1][i][j] = vr*advq +dfsq;
      scadv = mflow*ic[e][1][1][i][4+j]*(ec1ip -ec2ip*ec3ip);
      sc[i] = sc[i] +vr*scadv;
    }
  }

  // Forming elemental stiffness matrix
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    rr = ((1.0-fl[n][1])*prs +fl[n][1]*prl)/(prl +1.0E-32);
    for (j = 1; j <= nnpe; ++j) {
      if (i == j) {aq[1][i][j] = vr*rr*ec[n][2]*ds[i][3]/df
       -flux[1][i][j] +flux[1][o[i]][j];}
      if (i != j) {aq[1][i][j] = -flux[1][i][j]
       +flux[1][o[i]][j];}
    }
	
    // Source terms (customized)
    sct = 0.0;
    ci = cip[i];
    ti = tmin +tip[i]*fabs(tmax-tmin);
    tm = fo*lr*lr*prl*phl/pkl;
    for (j = 1; j <= css[1][1]; ++j) {
      k = int(css[j+1][1]);
      if (k == 12) {custom(2, j, tm, ci, ti, csp, css,
       csb, sct);}																					// calls custom
    }

    // Transient and source terms (phase change + customized)
    rq[i] = -vr*rr*((ec[n][1] -ec[n][2]*ec[n][3])*ds[i][3]/df
     -(ec[n][4] +ec[n][5]*(tn[n][2]-ec[n][6]))*ds[i][3]/df)
     +sct*ds[i][3]*(lr*lr/(pkl*dtr*vr)) +(sc[o[i]] -sc[i]);
  }
}