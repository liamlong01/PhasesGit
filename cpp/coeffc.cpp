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
//-- Diffusion coefficients: two-phase concentration equation --//
void coeffc(int e, double cph[npe1][3], double rd[npe1][5])
{
  int i, j;
  double dds, ddl, tphs, tphl, cphs, cphl;
  double res[npe1][3][3];

  dds = pds0/(pdl0 +1.0E-32);
  ddl = pdl0/(pdl0 +1.0E-32);
  for (j = 1; j <= nnpe; ++j) {
    i = ie[e][j];

    // Binary phase diagram: solidus and liquidus lines
    tphs = fmax(tmlt -cn[i][1]*(tmlt -eps -tsol)/pcs -eps, tsol);
    tphl = fmax(tmlt -cn[i][1]*(tmlt +eps -tsol)/pcl +eps, tsol);
    cphs = pcs*(tmlt -eps -tn[i][1])/(tmlt -eps -tsol);
    cphl = pcl*(tmlt +eps -tn[i][1])/(tmlt +eps -tsol);

    // Phase compositions and resistances
    if (fl[i][1] == 0.0) {
      cph[j][1] = cn[i][1];
      cph[j][2] = 0.0;
      res[j][1][1] = 1.0/(dds +1.0E-32);
      res[j][1][2] = 1.0/(dds +1.0E-32);
      res[j][2][1] = 1.0E+08;
      res[j][2][2] = 1.0E+08;
    }
    else if (fl[i][1] == 1.0) {
      cph[j][1] = 0.0;
      cph[j][2] = cn[i][1];
      res[j][1][1] = 1.0E+08;
      res[j][1][2] = 1.0E+08;
      res[j][2][1] = 1.0/(ddl +1.0E-32);
      res[j][2][2] = 1.0/(ddl +1.0E-32);
    }
    else {
      cph[j][1] = cphs;
      cph[j][2] = cphl;
      res[j][1][1] = 1.0/(dds*(1.0-fl[i][1]) +1.0E-32);
      res[j][1][2] = 1.0/(dds*(1.0-fl[i][1]) +1.0E-32);
      res[j][2][1] = 1.0/(ddl*fl[i][1] +1.0E-32);
      res[j][2][2] = 1.0/(ddl*fl[i][1] +1.0E-32);
    }
  }

  // Diffusion coeffients: solid / liquid phases
  for (j = 1; j <= nnpe; ++j) {
    rd[j][1] = 2.0*ds[j][2]/(fabs(res[j][1][1])
     +fabs(res[o[j+2]][1][1]));
    rd[j][2] = 2.0*ds[j][1]/(fabs(res[j][1][2])
     +fabs(res[o[j+2]][1][2]));
    rd[j][3] = 2.0*ds[j][2]/(fabs(res[j][2][1])
     +fabs(res[o[j+2]][2][1]));
    rd[j][4] = 2.0*ds[j][1]/(fabs(res[j][2][2])
     +fabs(res[o[j+2]][2][2]));
  }
}																											// calls no functions 