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
//-- Shape functions: bilinear quadrilateral elements --//
void shape(int e, int nnpm, double &dx, double &dy,
 double s[4][3], double dn[npe1][5], double ivo[npe1],
 double &jac, double &dqdx, double &dqdy)
{
  int j;
  double dxds, dxdt, dyds, dydt;

  // Derivatives dn/ds and dn/dt
  dn[1][1] = 0.25*(1.0 +s[1][2]);
  dn[2][1] = -0.25*(1.0 +s[1][2]);
  dn[3][1] = -0.25*(1.0 -s[1][2]);
  dn[4][1] = 0.25*(1.0 -s[1][2]);
  dn[1][2] = 0.25*(1.0 +s[1][1]);
  dn[2][2] = 0.25*(1.0 -s[1][1]);
  dn[3][2] = -0.25*(1.0 -s[1][1]);
  dn[4][2] = -0.25*(1.0 +s[1][1]);
  dxds = 0.0;
  dxdt = 0.0;
  dyds = 0.0;
  dydt = 0.0;

  // Global - local coordinate transformation
  for (j = 1; j <= nnpe; ++j) {
    dxds = dxds +x[ie[e][j]][1]*dn[j][1];
    dxdt = dxdt +x[ie[e][j]][1]*dn[j][2];
    dyds = dyds +x[ie[e][j]][2]*dn[j][1];
    dydt = dydt +x[ie[e][j]][2]*dn[j][2];
  }
  dx = dxds*(s[2][1] -s[3][1]) +dxdt*(s[2][2] -s[3][2]);
  dy = dyds*(s[2][1] -s[3][1]) +dydt*(s[2][2] -s[3][2]);
  jac = dxds*dydt -dyds*dxdt;

  // Derivatives dn/dx and dn/dy
  for (j = 1; j <= nnpe; ++j) {
    dn[j][3] = (dydt*dn[j][1] -dyds*dn[j][2])/(jac +1.0E-32);
    dn[j][4] = (-dxdt*dn[j][1] +dxds*dn[j][2])/(jac +1.0E-32);
  }
  dqdx = 0.0;
  dqdy = 0.0;

  // Scalar derivatives within element
  for (j = 1; j <= nnpe; ++j) {
    dqdx = dqdx +dn[j][3]*ivo[j];
    dqdy = dqdy +dn[j][4]*ivo[j];
  }
}																												// calls no functions