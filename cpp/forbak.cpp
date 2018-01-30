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
//-- Banded matrix solver by forward - backward substitution --//
void forbak(int m, double c[3*nnpm][6*nym+11], double z[3*nnpm],
 double r[3*nnpm], int nym, int &neq, int nnpm, int band3)
{
  int jw, isemi, iband, jpiv, lim1, lim2, iloc, jcol;
  double sum;

  // Forward substitution
  iband = int(m*(band3 +m -3)/3);
  isemi = (iband +1)/2;
  z[1] = r[1];
  sum = 0.0;
  for (jpiv = 2; jpiv <= neq; ++jpiv) {
    lim1 = 1;
    if (jpiv < isemi) {lim1 = isemi -jpiv +1;}
    lim2 = isemi -1;
    sum = r[jpiv];
    iloc = jpiv -isemi;
    for (jcol = lim1; jcol <= lim2; ++jcol) {
      sum = sum -c[jpiv][jcol]*z[iloc+jcol];
    }
    z[jpiv] = sum;
  }

  // Back substitution
  z[neq] = z[neq]/(c[neq][isemi] +1.0E-32);
  for (jw = 2; jw <= neq; ++jw) {
    jpiv = neq -jw +1;
    lim1 = isemi +1;
    lim2 = iband;
    if (jw < isemi) {lim2 = isemi +jw -1;}
    sum = z[jpiv];
    iloc = jpiv -isemi;
    for (jcol = lim1; jcol <= lim2; ++jcol) {
      sum = sum -c[jpiv][jcol]*z[iloc+jcol];
    }
    z[jpiv] = sum/(c[jpiv][isemi] +1.0E-32);
  }
}																												// calls no functions 