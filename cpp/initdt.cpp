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
//#include "functInit.h"  //Declares all functions

///////////////////////////////////////////////////////////////////////////////
//-- Initial conditions and problem parameters --/
void initdt(int nnp, int nel, int nsrf, int &band1,
 int nnb[nnpm][11], int ph[nnpm][3], double ec[nnpm][7],
 double fl[nnpm][3], double tn[nnpm][4],
 double vi[nelm][npe1][7], double pr[nelm][npe1])
{
  int e, i, j, k, n, m, iw, i1, i2;
  int dif, maxdif, j1, j2, j3, j4;
  double cphs, cphl, tphs, tphl, um, vm;
  double ste, hbar, u1, v1, vdn;

  // Start from initial conditions
  for (e = 1; e <= nel; ++e) {
    for (n = 1; n <= nnpe; ++n) {
      i = ie[e][n];
      vi[e][n][1] = vps[i][1];
      vi[e][n][2] = vps[i][2];
      vi[e][n][5] = vps[i][1];
      vi[e][n][6] = vps[i][2];
      pr[e][n] = (vsc*prl*phl)/(pkl +1.0E-32);
    }
  }

  // Initial phase distributions
  ste = phl*fabs(tmax-tmin)/(pl + 1.0E-32);
  hbar = 0.5*(phs+phl);
  for (i = 1; i <= nnp; ++i) {
    tphs = tmlt -eps -cn[i][1]*(tmlt -eps -tsol)/(pcs+1.0E-32);
    tphl = tmlt +eps -cn[i][1]*(tmlt +eps -tsol)/(pcl+1.0E-32);
    cphs = pcs*(tmlt -eps -tn[i][1])/(tmlt -eps -tsol);
    cphl = pcl*(tmlt +eps -tn[i][1])/(tmlt +eps -tsol);

    // Solid phase reference states
    if (tn[i][1] <= tphs) {
      ph[i][1] = 1;
      fl[i][1] = 0.0;
      ec[i][1] = 0.0;
      ec[i][2] = phs/(phl +1.0E-32);
      ec[i][3] = 0.0;
    }

    // Liquid phase reference states
    else if (tn[i][1] >= tphl) {
      ph[i][1] = 3;
      fl[i][1] = 1.0;
      ec[i][1] = hbar*(tphl-tphs)/(phl +1.0E-32)
       +(phs*tphs)/(phl +1.0E-32) +1.0/(ste +1.0E-32);
      ec[i][2] = phl/(phl +1.0E-32);
      ec[i][3] = tphl;
    }

    // Melt phase reference states
    else {
      ph[i][1] = 2;
      ec[i][1] = (phs*tphs)/(phl +1.0E-32);
      ec[i][2] = hbar/(phl +1.0E-32) +1.0/(ste*(tphl -tphs)
       +1.0E-32);
      ec[i][3] = tphs;
      fl[i][1] = (tn[i][1] -tphs)/(tphl -tphs);
    }

	//cout << fl[i][1] << endl;
  }

  // Droplet velocities
  if (ao == 2) {
    um = 0.0;
    vm = 0.0;
    for (i = 1; i <= nsrf; ++i) {
      for (k = 1; k <= 2; ++k) {
        j = 2*(i -1) +k;
        u1 = sqrt(um*um +vm*vm);
        v1 = sqrt(pow(bc[j][3][3]/(bc[j][3][2] +1.0E-32), 2.0)
         +pow(bc[j][4][3]/(bc[j][4][2] +1.0E-32), 2.0));
        if ((bc[j][3][1] == 0) && (bc[j][4][1] == 0)
         && (v1 > u1)) {
          um = bc[j][3][3]/(bc[j][3][2] +1.0E-32);
          vm = bc[j][4][3]/(bc[j][4][2] +1.0E-32);
        }
      }
    }

    for (i = 1; i <= nel; ++i) {
      for (k = 1; k <= nnpe; ++k) {
        vi[i][k][5] = um;
        vi[i][k][6] = vm;
      }
    }

    // Phase distribution
    for (i = 1; i <= nnp; ++i) {
      fl[i][1] = 1.0;
      ph[i][1] = 3;
    }

    for (i = 1; i <= nsrf; ++i) {
      n = bgn[i][1];
      m = bgn[i][2];
      vdn = um*(x[m][2] -x[n][2]) +vm*(x[m][1] -x[n][1]);
      for (k = 1; k <= 2; ++k) {
        j = 2*(i -1) +k;
        if ((bc[j][3][1] == 0) && (bc[j][3][3] == 0) &&
         (bc[j][4][1] == 0) && (bc[j][4][3] == 0) &&
         (vdn > 0.0)) {
          for (n = 1; n <= nnpe; ++n) {
            m = ie[bel[i]][n];
            fl[m][1] = 0.0;
            ph[m][1] = 2;
          }
        }
      }
    }
  }

  // Find matrix bandwidth
  maxdif = 0;
  for (n = 1; n <= nel; ++n) {
    for (i = 1; i <= nnpe; ++i) {
      for (k = 1; k <= nnpe; ++k) {
        dif = abs(ie[n][i] -ie[n][k]);
        if (dif > maxdif) {maxdif = dif;}
      }
    }
  }
  band1 = 2*(maxdif +1) -1;

  // Neighbour node array
  for (i = 1; i <= nnp; ++i) {
    nnb[i][1] = 1;
    nnb[i][2] = i;
  }
  for (e = 1; e <= nel; ++e) {
    for (j1 = 1; j1 <= nnpe; ++j1) {
      i1 = ie[e][j1];
      for (j2 = 1; j2 <= nnpe; ++j2) {
        i2 = ie[e][j2];
        iw = 0;
        j4 = 2;

        // Global node already in list
        for (j3 = 2; j3 <= nnb[i1][1]+1; ++j3) {
          if (nnb[i1][j3] == i2) {
            iw = 1;
            j4 = j3 -1;
          }
        }

        // Node not in list and not diagonal node
        if ((iw == 0) && (nnb[i1][2] != i2)) {
          nnb[i1][1] = nnb[i1][1] +1;
          nnb[i1][nnb[i1][1]+1] = i2;
          j4 = nnb[i1][1];
        }
      }
    }
  }
}



																									// calls no functions