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
//-- Control Volume entropy production rate --//
void entropy(double df, int nel, int nnp, int nsrf,
 double vps[nnpm][14])
{
  int e, i, j, k, n, jb, m, n1, n2, n3;
  double rhon, rhoip, kip, dx, dy, jac, tdim, ss, tt;
  double dtdx, dtdy, dudx, dudy, dvdx, dvdy, sf, ur, dtm;
  double tphs, tphl, teut, sr1, sr2, sr3, cpn;
  double tpip[npe1], sip[npe1], flux[4][npe1][npe1];

  ur = vr*pkl/(prl*phl*lr +1.0E-32);
  dtm = df*lr*lr*prl*phl/pkl;
  for (n = 1; n <= nnp; ++n) {
    for (i = 8; i <= 11; ++i) {
      vps[n][i] = 0.0;
    }
  }

  // Entropy equation of state
  teut = tmin +dtr*tsol;
  sf = pl/(tmin +dtr*tmlt +1.0E-32);
  for (n = 1; n <= nnp; ++n) {
    tdim = tmin +dtr*tn[n][1];
    tphs = tmin +dtr*(tmlt -eps -cn[n][1]*(tmlt -eps -tsol)
     /(pcs -eps));
    tphl = tmin +dtr*(tmlt +eps -cn[n][1]*(tmlt +eps -tsol)
     /(pcl +eps));

    // Solid phase reference states
    if (tn[n][1] <= tphs) {
      sr1 = 0.0;
      sr2 = phs;
      sr3 = teut;
    }

    // Liquid phase reference states
    else if (tn[n][1] >= tphl) {
      sr1 = log(tphl/tphs)*(phs*tphl -phl*tphs)/(tphl -tphs)
       +phs*log(tphs/teut) +(phl -phs +sf);
      sr2 = phl;
      sr3 = tphl;
    }

    // Melt phase reference states
    else {
      sr1 = phs*log(tphs/teut);
      sr2 = (phs*tphl -phl*tphs)/(tphl -tphs) +(phl -phs +sf)
       /log(tphl/tphs);
      sr3 = tphs;
    }
    vps[n][8] = sr1 +sr2*log(tdim/(sr3 +1.0E-16));
  }

  // Assembly of sub-control-volumes
  for (e = 1; e <= nel; ++e) {
    for (j = 1; j <= nnpe; ++j) {
      n = ie[e][j];
      ivo[j] = tn[n][1];
    }

    // Local element geometry
    for (i = 1; i <= nnpe; ++i) {
      m = ie[e][i];
      s[1][1] = lc[12+i][1];
      s[1][2] = lc[12+i][2];
      s[2][1] = lc[4+i][1];
      s[2][2] = lc[4+i][2];
      s[3][1] = lc[8+i][1];
      s[3][2] = lc[8+i][2];
      shape(e, nnpm, dx, dy, s, dn, ivo, jac, dtdx, dtdy);									// calls shape
      ds[i][3] = fabs(jac);
      ds[i][4] = fabs(jac);
      s[1][1] = lc[i][1];
      s[1][2] = lc[i][2];
      vps[m][9] = vps[m][9] +ds[i][3]*lr*lr;
      shape(e, nnpm, dx, dy, s, dn, ivo, jac, dtdx, dtdy);									// calls shape
      ds[i][1] = dx;
      ds[i][2] = dy;

      // Integration point interpolation
      ss = lc[i][1];
      tt = lc[i][2];
      nn[1][i] = 0.25*(1.0 +ss)*(1.0 +tt);
      nn[2][i] = 0.25*(1.0 -ss)*(1.0 +tt);
      nn[3][i] = 0.25*(1.0 -ss)*(1.0 -tt);
      nn[4][i] = 0.25*(1.0 +ss)*(1.0 -tt);
      rhoip = 0.0;
      kip = 0.0;
      tpip[i] = 0.0;
      sip[i] = 0.0;
      for (j = 1; j <= nnpe; ++j) {
        n = ie[e][j];
        rhoip = rhoip +nn[j][i]*((1.0 -fl[n][1])*prs
         +fl[n][1]*prl);
        kip = kip +nn[j][i]*((1.0 -fl[n][1])*pks +fl[n][1]*pkl);
        tpip[i] = tpip[i] +nn[j][i]*(tmin +dtr*tn[n][1]);
        sip[i] = sip[i] +nn[j][i]*vps[n][8];
      }

      // Entropy flux calculation
      flux[1][i][1] = kip*dtr*(dtdx*ds[i][2] -dtdy*ds[i][1])
       /tpip[i] +vr*lr*rhoip*sip[i]*(vi[e][i][1]*ds[i][2]
       -vi[e][i][2]*ds[i][1]);
    }

    // SCV entropy production rate
    for (i = 1; i <= nnpe; ++i) {
      n = ie[e][i];
	  rhon = (1.0 - fl[n][1])*prs + fl[n][1] * prl;														
      cpn = (1.0 -fl[n][1])*phs +fl[n][1]*phl;
      tdim = tmin +dtr*tn[n][1];
      vps[n][10] = vps[n][10] +rhon*cpn*dtr*(tn[n][1]
       -tn[n][2])*ds[i][3]*lr*lr/(dtm*tdim) +rhon*sf
       *(fl[n][1] -fl[n][2])*ds[i][3]*lr*lr/dtm
       -flux[1][i][1] +flux[1][o[i]][1];
      vps[n][11] = tdim*tdim/(dtr*dtdx*dtdx/lr
       +dtr*dtdy*dtdy/lr +1.0E-16);
	 // cout << rhon << endl;
	}
  }

  // Entropy based diffusivity
  for (n = 1; n <= nnp; ++n) {
    vps[n][10] = vps[n][10]/(vps[n][9] +1.0E-32);
    if (vps[n][10] < 0.0) {
      vps[n][11] = fabs(vps[n][10])*vps[n][11];
			if (vps[n][11] > 1.0) {vps[n][11] = 1.0;}
    }
    else {
      vps[n][11] = 0.0;
    }
  }

  // Boundary entropy production rate
  for (n = 1; n <= nsrf; ++n) {
    e = bel[n];
    n1 = bgn[n][1];
    n2 = bgn[n][2];
    for (i = 1; i <= 2; ++i) {
      for (k = 1; k <= 3; ++k) {
        for (j = 1; j <= nnpe; ++j) {
          n3 = ie[e][j];
          if (n1 == n3) jb = j;
          if (k == 1) ivo[j] = tn[n3][1];
          if (k == 2) ivo[j] = vps[n3][1];
          if (k == 3) ivo[j] = vps[n3][2];
        }
        m = 2*(jb -1) +i;
        s[1][1] = nc[m][1];
        s[1][2] = nc[m][2];
        s[2][1] = nc[8+m][1];
        s[2][2] = nc[8+m][2];
        s[3][1] = nc[16+m][1];
        s[3][2] = nc[16+m][2];
        if (k == 1) {shape(e, nnpm, dx, dy, s, dn, ivo,
         jac, dtdx, dtdy);}																		// calls shape
        if (k == 2) {shape(e, nnpm, dx, dy, s, dn, ivo,
         jac, dudx, dudy);}																		// calls shape
        if (k == 3) {shape(e, nnpm, dx, dy, s, dn, ivo,
         jac, dvdx, dvdy);}																		// calls shape
      }
      j = bgn[n][i];
      tdim = tmin +dtr*tn[j][1];
      vps[j][10] = pks*dtr*dtr*(dtdx*dtdx +dtdy*dtdy)
       /(tdim*tdim*lr*lr) +(prl*vsc)*vr*vr*(2.0*dudx*dudx
       +2.0*dvdy*dvdy +(dudy +dvdx)*(dudy +dvdx)
       -0.667*(dudx +dvdy)*(dudx +dvdy))/(tdim*lr*lr);
      vps[j][11] = 0.0;
    }
  }
}