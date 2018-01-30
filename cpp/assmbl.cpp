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
//-- Assembly of elemental stiffness equations --//
void assmbl(int v, int p, int m, int nel, int nelm,
 int band3, double df, double fo, double nn[npe1][npe1],
 double aq[4][npe1][npe1], double c[3*nnpm][6*nym+11],
 double rq[npe1], double r[3*nnpm])
{
  int e, i, j, n, iw, jw, im, iband, isemi;
  double jac, dx, dy, ss, tt, dqdx, dqdy;
  double xip[npe1][3], cip[npe1], tip[npe1];
  static double sc[npe1], flux[4][npe1][npe1];

  iband = int(m*(band3 +m -3)/3);
  isemi = (iband +1)/2;

  // Start outer loop over all elements
  for (e = 1; e <= nel; ++e) {
    for (j = 1; j <= nnpe; ++j) {
      s[1][1] = lc[12+j][1];
      s[1][2] = lc[12+j][2];
      s[2][1] = lc[4+j][1];
      s[2][2] = lc[4+j][2];
      s[3][1] = lc[8+j][1];
      s[3][2] = lc[8+j][2];
      shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);												// calls shape
      ds[j][3] = fabs(jac);
      ds[j][4] = fabs(jac);
      s[1][1] = lc[j][1];
      s[1][2] = lc[j][2];
      shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);												// calls shape
      ds[j][1] = dx;
      ds[j][2] = dy;

      // Integration point coordinates
      ss = lc[j][1];
      tt = lc[j][2];
      nn[1][j] = 0.25*(1.0 +ss)*(1.0 +tt);
      nn[2][j] = 0.25*(1.0 -ss)*(1.0 +tt);
      nn[3][j] = 0.25*(1.0 -ss)*(1.0 -tt);
      nn[4][j] = 0.25*(1.0 +ss)*(1.0 -tt);
      xip[j][1] = 0.0;
      xip[j][2] = 0.0;
      cip[j] = 0.0;
      tip[j] = 0.0;
      for (i = 1; i <= nnpe; ++i) {
        n = ie[e][i];
        xip[j][1] = xip[j][1] +nn[i][j]*x[n][1];
        xip[j][2] = xip[j][2] +nn[i][j]*x[n][2];
        cip[j] = cip[j] +nn[i][j]*cn[n][1];
        tip[j] = tip[j] +nn[i][j]*tn[n][1];
      }
    }

    // Forming stiffness matrices
    if (v == 1) {stiffc(e, df, fo, xip, flux, sc, aq, rq);}												// calls stiffc
    if (v == 2) {stifft(e, df, fo, cip, tip, flux, sc, aq, rq);}										// calls stifft
    if (v == 3) {stiffu(e, df, xip, flux, sc, aq, rq);}													// calls stiffu
    if (v == 4) {stiffv(e, df, flux, sc, aq, rq);}														// calls stiffv
    if (v == 5) {stiffp(e, flux, sc, aq, rq);}															// calls stiffp
    if (v == 6) {stiffdv(e, vi);}																		// calls stiffdv
    if (v == 7) {stiffdc(e, df, xip, flux, aq, rq);}													// calls stiffdc

    // Local stiffness matrix -> global matrix
    for (iw = 1; iw <= nnpe; ++iw) {
      i  = ie[e][iw];
      im = p +m*(ie[e][iw] -1);
      r[im] = r[im] +rq[iw];
      for (jw = 1; jw <= nnpe; ++jw) {
        j = ie[e][jw];
        n = isemi +m*(j -i);
        c[im][n+1-p] = c[im][n+1-p] +aq[1][iw][jw];
        c[im][n+2-p] = c[im][n+2-p] +aq[2][iw][jw];
        c[im][n+3-p] = c[im][n+3-p] +aq[3][iw][jw];
      }
    }
  }
}
