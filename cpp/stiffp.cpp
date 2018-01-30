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
//-- Local stiffness matrices: u-v-P equation (continuity) --//
void stiffp(int e, double flux[4][npe1][npe1], double sc[npe1],
 double aq[4][npe1][npe1], double rq[npe1])
{
  int i, j, n, m;
  double flip, rr, advu, advv, advp;
  double dqdx, dqdy, jac, dx, dy;

  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i]; 
    m = ie[e][o[i]];
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    sc[i] = 0.0;
    flip = 2.0*fl[n][1]*fl[m][1]/(fl[n][1] +fl[m][1] +1.0E-32);
    rr = ((1.0-flip)*prs +flip*prl)/(prl +1.0E-32);
    shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);								// calls shape 

    // Mass flux across 2-D sub-surface
    for (j = 1; j <= nnpe; ++j) {
      advu = ic[e][1][1][i][4+j]*ds[i][2];
      advv = -ic[e][2][1][i][4+j]*ds[i][1];
      advp = ic[e][1][2][i][4+j]*ds[i][2] -ic[e][2][2][i][4+j]
       *ds[i][1];
	  flux[1][i][j] = advu; 
      flux[2][i][j] = advv;
      flux[3][i][j] = advp; 
    }
  }

  // Forming elemental stiffness matrix
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    for (j = 1; j <= nnpe; ++j) {
      aq[1][i][j] = -flux[1][i][j] +flux[1][o[i]][j];
      aq[2][i][j] = -flux[2][i][j] +flux[2][o[i]][j];
      aq[3][i][j] = -flux[3][i][j] +flux[3][o[i]][j];
    }
    rq[i] = 0.0;
  }
}