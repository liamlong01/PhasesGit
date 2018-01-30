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
//-- Re-initialize matrices --//
void nullmx(double r[3*nnpm], double z[3*nnpm],
 double c[3*nnpm][6*nym+11], double aq[4][npe1][npe1])
{
  int i, j, v;
  for (i = 0; i < 3*nnpm; ++i) {
    r[i] = 0.0;
    z[i] = 0.0;
  }

  for (i = 0; i < 3*nnpm; ++i) {
    for (j = 0; j < 6*nym+11; ++j) {
      c[i][j] = 0.0;
    }
  }

  for (v = 0; v < 4; ++v) {
    for (i = 0; i < npe1; ++i) {
      for (j = 0; j < npe1; ++j) {
        aq[v][i][j] = 0.0;
      }
    }
  }
}																											// calls no functions