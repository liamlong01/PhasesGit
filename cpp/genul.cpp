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
//-- LU decomposition into lower - upper triangular matrices --//
void genul(int m, double c[3*nnpm][6*nym+11], int nym,
 int &neq, int nnpm, int band3)
{
  int isemi, iband, hirow, lowrow, neqm1;
  int irow, iloc, ieq, jcol, locol, hicol;
  double piv, factor;

  // Begin row reduction
  iband = int(m*(band3 +m -3)/3);
  isemi = (iband +1)/2;
  neqm1 = neq -1;
  for (ieq = 1; ieq <= neqm1; ++ieq) {
    lowrow = ieq +1;
    hirow = ieq +isemi -1;
    if (hirow > neq) {hirow = neq;}
    piv = c[ieq][isemi];
    for (irow = lowrow; irow <= hirow; ++irow) {
      iloc = irow -ieq;
      factor = c[irow][isemi-iloc]/(piv +1.0E-16);
      c[irow][isemi-iloc] = factor;
      locol = isemi -iloc +1;
      hicol = locol +isemi -2;
      for (jcol = locol; jcol <= hicol; ++jcol) {
        c[irow][jcol] = c[irow][jcol] -factor*c[ieq][jcol+iloc];
      }
    }
  }
}																										// calls no functions 