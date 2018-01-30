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
//-- Customized properties, source terms and boundary conditions --//
void custom(int n, int m, double &tm, double &ci, double &ti,
 double csp[11][7], double css[11][7], double csb[4*nym][8],
 double &cval)
{
  int g, j, k;
  double cs[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  j = 0;
  for (g = 1; g <= 5; ++g) {
    if (n == 1) {cs[g] = csp[m+1][g];}
    if (n == 2) {cs[g] = css[m+1][g];}
    if (n == 3) {cs[g] = csb[m+1][g+2];}
  }
  if (n <= 3) k = int(cs[2]);
  if (n == 1) {
    if ((k == 2) || (k == 3)) {k = int(cs[3]);}
    j = 1;
  }
  if (k == 1) {cval = cs[3+j] +tm*cs[4+j] +tm*tm*cs[5+j];}
  if (k == 11) {cval = cs[3+j] +ci*cs[4+j] +ci*ci*cs[5+j];}
  if (k == 12) {cval = cs[3+j] +ti*cs[4+j] +ti*ti*cs[5+j];}
}																											// calls no functions 