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
//-- Updated integration point velocities --//
void uvcor(int nel, double vi[nelm][npe1][7])
{
  int e, n, i, j;
  double uic, vic, pxc, pyc;

  // Re-construct momentum equations
  for (e = 1; e <= nel; ++e) {
    for (j = 1; j <= nnpe; ++j) {
      n = ie[e][j];
      uic = 0.0;
      vic = 0.0;
      pxc = 0.0;
      pyc = 0.0;
      for (i = 1; i <= nnpe; ++i) {
        n = ie[e][i];
        uic = uic +ic[e][1][1][j][4+i]*vps[n][1];
        vic = vic +ic[e][2][1][j][4+i]*vps[n][2];
        pxc = pxc +ic[e][1][2][j][4+i]*vps[n][3];
        pyc = pyc +ic[e][2][2][j][4+i]*vps[n][3];
      }
      vi[e][j][1] = uic +pxc;
      vi[e][j][2] = vic +pyc;
    }
  }
}																													// calls no functions 