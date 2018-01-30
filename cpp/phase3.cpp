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
//-- Surface: three-phase heat, mass / momentum balances --//
void phase3(int nel, int nnp, int &cmf, int &cmx, int &cr,
 int nsrf, double df, int ph[nnpm][3], double fl[nnpm][3],
 double cn[nnpm][6], double ec[nnpm][7])
{
  int e, i, j, k, m, n, phj, ns, nf, rn, rg;
  double hconv, u0, g0, hf, dx, dy, vx1, vx2, cl;
  double ds, dm, di, t0, ts, tf, tm, ti, tp, tw;
  double dudx, dudy, chi, e0, beta, cw, at, bt;
  double dt, tauf, hs, taus, mr, flr, dtm, ar2;

  // Initialized problem parameters
  rg = 0; // rime (0) or rime / glaze (1) ice
  g0 = 1.0E+08; // initialized liquid water content
  u0 = 0.0; // initialized reference velocity
  dtm = df*lr*lr*prl*phl/pkl; // dimensional time step
  cr = 0; // phase convergence parameter
  cmx = 2; // number of iterations parameter
  hconv = 500.0; // Comparison with Myers (IJHMT, 1999)
  beta = 0.5; // empirical collection efficiency

  // Reference properties: air / water / ice system
  ur = vr*pkl/(prl*phl*lr +1.0E-32); // dimensional velocity
  tp = 273.0; // phase change temperature of water
  dt = 1.0; // characteristic temperature difference
  chi = 8.52; // evaporation coefficient
  e0 = 44.4; // vapor pressure constant
  cw = 4220.0; // specific heat of unfrozen water
  tf = (tp -tmin)/(tmax -tmin +1.0E-32);
  tw = tp; // reference wall temperature
  for (i = 1; i <= nsrf; ++i) {
    for (k = 1; k <= 2; ++k) {
      j = 2*(i -1) +k;
      if ((bc[j][3][1] == 0) && (bc[j][3][3] == 0) &&
       (bc[j][4][1] == 0) && (bc[j][4][3] == 0) &&
       (bc[j][2][2] != 0)) {
        tw = fmin(tw, tmin +(tmax -tmin)*bc[j][2][3]
         /(bc[j][2][2] +1.0E-32));
      }
    }
  }

  // Mass conservation: incoming droplets
  for (e = 1; e <= nel; ++e) {
    for (j = 1; j <= nnpe; ++j) {
      i = ie[e][j];
      tn[i][3] = 0.0;
      cn[i][1] = fmax(cn[i][1], 0.0);

      // Change of phase and mass in filled volumes
      if (cn[i][4] > 0.0) {cn[i][1] = cn[i][4];}
      cn[i][4] = 0.0;
      if (cmf == 1) {cmx = fmax(cmx, int(cn[i][1]));}
      cn[i][1] = fmin(cn[i][1], 2.0);
      if (cn[i][1] >= 1.0) {
        ph[i][1] = 1;
        for (k = 1; k <= nnpe; ++k) {
          n = ie[e][k];

          // Mass re-adjustment in over-filled volumes
          if (ph[n][1] != 2) {fl[n][1] = 0.0;}
          if ((cn[n][1] > 0.0) && (cn[n][1] < 0.999)) {
            ph[n][1] = 2;
            cn[n][4] = fabs(cn[i][1] -1.0);
          }
        }
        cn[i][2] = 1.0;
      }
      if (ph[i][1] != ph[i][2]) {cr = 1;}

      // Reference parameters
      g0 = fmax(fmin(1000.0*cn[i][1], g0), 1.0E-08);
      u0 = fmax(ur*sqrt(vps[i][1]*vps[i][1] 
       +vps[i][2]*vps[i][2]), u0);

      // Initial solid layer temperature
      if (ph[i][1] == 1) {
        for (k = 1; k <= nnpe; ++k) {
          n = ie[e][k];
          tn[i][2] = fmin(tn[n][2], tn[i][2]);
          if (ph[n][1] != 1) {tn[n][2] = tn[i][2];}
        }
      }
    }
  }

  // Neighbour filled volumes
  for (i = 1; i <= nnp; ++i) {
    rn = 0;
    if ((cn[i][1] > 0.999) && (ph[i][1] == 1))
     {cn[i][1] = 1.0;}
    if (nnb[i][1] == 9) {
      for (j = 2; j <= nnb[i][1] +1; ++j) {
        n = nnb[i][j];
        if ((cn[i][1] < 1) && (fl[n][1] == 0)) {rn = rn++;}
      }
      if (rn >= 8) {
        cn[i][1] = 1.0;
        fl[i][1] = 0.0;
        ph[i][1] = 1;
      }
    }
  }

  // Start heat / momentum balances
  for (e = 1; e <= nel; ++e) {
    phj = 0;
    for (j = 1; j <= nnpe; ++j) {
      i = ie[e][j];
      if (ph[i][1] == 1) {phj = 1;}
      if (ph[i][1] == 3) {phj = 3;}
      if ((ph[i][1] == 2) && (ph[ie[e][o[j]]][1] == 2)
       && (ph[ie[e][o[j+2]]][1] == 2)) {phj = 1;}
    }

    // Loop over sub-control volumes
    for (j = 1; j <= nnpe; ++j) {
      i = ie[e][j];
      n = ie[e][o[j]];
      m = ie[e][o[j+2]];
      ns = n;
      nf = m;

      // Tangential / normal directions: phase interface
      vx1 = vps[i][1]*(x[i][1] -x[n][1]) +vps[i][2]
       *(x[i][2] -x[n][2]);
      vx2 = vps[i][1]*(x[i][1] -x[m][1]) +vps[i][2]
       *(x[i][2] -x[m][2]);
      if (fabs(vx1) > fabs(vx2)) {
        ns = m;
        nf = n;
      }
      t0 = tmin +(tmax -tmin)*tn[nf][1];
      tm = fmax(tmin +(tmax -tmin)*tn[ns][1], tw);
      if ((vx1 < 0) || (vx2 < 0)) {
        t0 = tmin +(tmax -tmin)*tn[i][1];
        tm = fmax(tmin +(tmax -tmin)*tn[nf][1], tw);
      }

      // Rime ice
      if (rg == 0) {
        cn[i][5] = 1.0;
      }

      // Inside the solid
      else if (ph[i][1] == 1) {
        ec[i][2] = phs/(phl +1.0E-32);
      }

      // Solid side of phase interface with unfrozen water
      else if ((phj == 1) && (ph[i][1] == 2)) {
        if (cn[m][1] > cn[n][1]) {n = m;}
        dm = lr*sqrt(pow(x[nf][1] -x[i][1], 2.0)
         +pow(x[nf][2] -x[i][2], 2.0));
        at = fmin(fmax((hconv*(tp -t0) +0.333*pks*(tp -tm)/dm)
         /(chi*e0*(tp -t0) +0.5*pks*(tp -tm)/dm +beta*u0*g0
         *cw*dt +hconv*dt +1.0E-32), 0.0), 1.0);

        // Add conduction to heat balance (rime; cn[i][5] = 1)
        ds = lr*sqrt(pow(x[i][1] -x[n][1], 2.0)
         +pow(x[i][2] -x[n][2], 2.0));
        di = 0.5*ds*(1.0 +2.0*fmin(0.5, cn[n][1]));
        ti = at*tf +(1.0 -at)*tn[i][1];
        cn[i][5] = fmin(fmax(cn[i][5] +fabs((tmax -tmin)*pks
         *(tf -ti)/(di*u0*g0*pl +1.0E-32)), 0.0), 1.0);
      }

      // Air side of phase interface (incoming droplets)
      else if ((phj == 3) && (ph[i][1] == 2)) {
        ds = lr*sqrt(pow(x[ns][1] -x[i][1],2.0)
         +pow(x[ns][2] -x[i][2],2.0));
        s[1][1] = lc[12+j][1];
        s[1][2] = lc[12+j][2];
        s[2][1] = lc[4+j][1];
        s[2][2] = lc[4+j][2];
        s[3][1] = lc[8+j][1];
        s[3][2] = lc[8+j][2];
        shape(e, nnpm, dx, dy, s, dn, ivo, ar2, dudx, dudy);												// calls shape 

        // Film momentum balance for shear-driven runback
        hf = fl[nf][1]*lr*lr*ar2/(0.5*ds +1.0E-32);
        tauf = 0.5*(0.008 +2.0E-05*(u0*hf/vsc))*prl*u0*u0;
        hs = fl[ns][1]*lr*lr*ar2/(0.5*ds +1.0E-32);
        taus = 0.5*(0.008 +2.0E-05*(u0*hs/vsc))*prl*u0*u0;
        mr = 0.5*fabs(tauf*hf*hf -taus*hs*hs)*dtm/vsc;
        flr = mr/(1000.0*lr*lr*ar2 +1.0E-32);
        fl[i][1] = fmax(fmin(fl[i][1] +flr, 1.0 -cn[i][1]), 0.0);

        // Laminar flow, smooth surface of arbitrary shape
        //hconv = 0.2926*pkl*pow(u0/(vsc*ds), 0.5);

        // Turbulent flow, smooth surface of arbitrary shape
        //hconv = 0.0287*prl*phl*u0*pow(vsc*prl*phl/pkl, -0.4)
        // *pow(vsc/(u0*ds), 0.2);

        // Temperature rise due to latent heat released
        ts = t0 +(u0*g0*pl +0.5*hconv*u0*u0/phl)
         /(hconv +u0*g0*phl);
        tn[i][3] = fmin((ts -tmin)/(tmax -tmin), tf);

        // Criterion for unfrozen liquid layer (glaze ice)
        if (ts >= tp) {
          bt = fmin(fmax((u0*g0*cw*dt -0.5*g0*u0*u0*u0
           +chi*e0*dt)/(u0*g0*pl +1.0E-32), 0.0), 1.0);
          cn[i][5] = fmin(fmax((cn[i][5] +hconv*(tp -t0))
           /(u0*g0*pl +1.0E-32) +bt, 0.0), 1.0);
        }
        ec[i][1] = 0.0;
        ec[i][2] = 1.0E+08;
        ec[i][3] = 0.0;
      }

      // Droplets in the freestream
      else if (ph[i][1] == 3) {
        ec[i][2] = 1.0;
      }
      ec[i][4] = ec[i][1];
      ec[i][5] = ec[i][2];
      ec[i][6] = ec[i][3];
    }
  }

  // Volume fraction change: surface runoff
  for (i = 1; i <= nnp; ++i) {
    if ((ph[i][1] == 2) && (tn[i][3] >= tf)) {
      cl = (1.0 -cn[i][5])*(cn[i][1] -cn[i][2]);
      fl[i][1] = fmin(fmax(cl, 0.0), 1.0);
      cn[i][1] = fmin(fmax(cn[i][1] -cl, 0.0), 1.0);
    }
    cn[i][5] = 0.0;
  }
}