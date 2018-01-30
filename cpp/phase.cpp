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
//-- Iteration rules for solid / liquid phase change --//
void phase(int nnp, int ph[nnpm][3], double fl[nnpm][3],
 double ec[nnpm][7], int &tr)
{
  int i, j, r2;
  double ste, tphs, tphl, cphs, cphl, hbar;

  hbar = 0.5*(phl +phs);
  ste = phl*dtr/(pl +1.0E-32);
  for (i = 1; i <= nnp; ++i) {
    tphs = fmax(tmlt -eps -cn[i][1]*(tmlt -eps -tsol)/pcs, tsol);
    tphl = fmax(tmlt +eps -cn[i][1]*(tmlt +eps -tsol)/pcl, tsol);
    cphs = pcs*(tmlt -eps -tn[i][1])/(tmlt -eps -tsol);
    cphl = pcl*(tmlt +eps -tn[i][1])/(tmlt +eps -tsol);
    if (tn[i][1] <= tphs) {
      ph[i][1] = 1;
      fl[i][1] = 0.0;
    }
    else if (tn[i][1] >= tphl) {
      ph[i][1] = 3;
      fl[i][1] = 1.0;
    }
    else {
      ph[i][1] = 2;
      fl[i][1] = (tn[i][1] -tphs)/(tphl -tphs);
    }
  }

  // Compare tentative vs computed phases
  tr = 0;
  for (i = 1; i <= nnp; ++i) {
    if (ph[i][1] != ph[i][2]) {tr = 1;}
  }

  // Rule 1: Pass through melt phase
  if (tr == 1) {
    for (i = 1; i <= nnp; ++i) {
      if (abs(ph[i][1] -ph[i][2]) > 1) {ph[i][1] = 2;}//Changed from "fabs" to "abs"
    }
  }

  // Rule 2: Adjacent CV phase change criterion
  if (tr == 1) {
    for (i = 1; i <= nnp; ++i) {
      r2 = 0;
      if (nnb[i][1] == 9) {
        for (j = 2; j <= nnb[i][1]+1; ++j) {
          if (ph[nnb[i][j]][2] != ph[i][2]) {r2 = 1;}
        }
        if ((r2 == 0) && (ph[i][1] != ph[i][2]))
          {ph[i][1] = ph[i][2];}
      }
    }
  }

  // Energy equation of state
  for (i = 1; i <= nnp; ++i) {
    tphs = tmlt -cn[i][1]*(tmlt-eps -tsol)/(pcs +1.0E-32) -eps;
    tphl = tmlt -cn[i][1]*(tmlt+eps -tsol)/(pcl +1.0E-32) +eps;
    cphs = pcs*(tmlt -eps -tn[i][1])/(tmlt -eps -tsol);
    cphl = pcl*(tmlt +eps -tn[i][1])/(tmlt +eps -tsol);

    // Solid phase reference states
    if (ph[i][1] == 1) {
      ec[i][1] = 0.0;
      ec[i][2] = phs/(phl +1.0E-32);
      ec[i][3] = 0.0;
    }

    // Liquid phase reference states
    else if (ph[i][1] == 3) {
      ec[i][1] = hbar*(tphl -tphs)/(phl +1.0E-32)
       +(phs*tphs)/(phl +1.0E-32) +1.0/(ste +1.0E-32);
      ec[i][2] = phl/phl;
      ec[i][3] = tphl;
    }

    // Melt phase reference states
    else {
      ec[i][1] = (phs*tphs)/(phl +1.0E-32);
      ec[i][2] = hbar/(phl +1.0E-32) +1.0/(ste*(tphl -tphs)
       +1.0E-32);
      ec[i][3] = tphs;
    }
  }
}																											// calls no functions 