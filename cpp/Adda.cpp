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
//-- Adda --//
double Adda(int w1, int w2, int w3, int w4, double w5, int nnp,
 int nel, int nx, int ny, int nsrf)

// Exported DLL -> Visual Basic GUI
// 1. New Project -> Win32 Project -> DLL
// 2. use: '#include <windows.h>' and 'double _declspec ...'
// 3. use: phasesb.def source file: LIBRARY phasesa
//                                  EXPORTS Adda
// 4. remove: 'void main()' and previous 'double Adda ...'
// 5. Debug -> Build -> c:\ ... \debug\phasesa.dll created
//#include <windows.h>
//double _declspec (dllexport) _stdcall Adda (int w1,
// int w2, int w3, int w4, double w5, int nnp, int nel,
// int nx, int ny, int nsrf)
{
  int i, n, band1;
  double cd, cl, um, addv;
  static double ivo[5] = {0.0,0.0,0.0,0.0};
  band3 = 6*ny +11;

  // Import problem data parameters
  if ((w1 > 0) && (w1 < 30)) {
    if (w1 == 1) {se = int(w5);}
    if (w1 == 2) {ao = int(w5);}
    if (w1 == 3) {tk = int(w5);}
    if (w1 == 4) {tvk = int(w5);}
    if (w1 == 5) {tstp = int(w5);}
	if (w1 == 6) {df = w5; } //cout << w5 << endl;}
	if (w1 == 7) {lr = w5; } // cout << w5 << endl;}
	if (w1 == 8) {wr = w5; } // cout << w5 << endl;}
	if (w1 == 9) {tol = w5;} //cout << w5 << endl;}
	if (w1 == 10) {ur = w5;} // cout << w5 << endl;}
    if (w1 == 11) {pcs = w5;}
    if (w1 == 12) {pcl = w5;}
	if (w1 == 13) {tmlt = w5;}
    if (w1 == 14) {tsol = w5;}
    if (w1 == 15) {bt = w5;}
    if (w1 == 16) {bs = w5;}
    if (w1 == 17) {tmin = w5;}
    if (w1 == 18) {tmax = w5;}
    if (w1 == 19) {vsc = w5;}
    if (w1 == 20) {prs = w5;}
    if (w1 == 21) {prl = w5;}
    if (w1 == 22) {pds0 = w5;}
    if (w1 == 23) {pdl0 = w5;}
    if (w1 == 24) {pks = w5;}
    if (w1 == 25) {pkl = w5;}
    if (w1 == 26) {phs = w5;}
    if (w1 == 27) {phl = w5;}
    if (w1 == 28) {pl = w5;}
    if (w1 == 29) {
		tmlt = (tmlt - tmin) / (tmax - tmin + 1.0E-32); //cout << tmlt << endl;
      tsol = (tsol -tmin)/(tmax -tmin +1.0E-32);
	  rat = 9.8*bt*(tmax - tmin)*lr*lr*lr*prl*phl / (vsc*pkl + 1.0E-32); //cout << rat << endl;
	  ras = 9.8*bs*(pcl - pcs)*lr*lr*lr*prl*phl / (vsc*pkl + 1.0E-32); //cout << ras << endl;
      pds0 = 1.0E-08*pds0;
      pdl0 = 1.0E-08*pdl0;
      dtr = fabs(tmax -tmin);

      // Set reference velocity based on boundary conditions
      // (boundary specified: ur, fixed wall: pkl/(prl*phl*lr))
      vr = ur/(pkl/(prl*phl*lr));
      if ((se <= 3) || (ao == 2)) {tvk = 1;}
    }
  }
  
  // Import problem mesh parameters
  if (w1 == 31) {ie[w2][w3] = w4;}
  if (w1 == 32) {x[w2][1] = w5/lr;}
  if (w1 == 33) {x[w2][2] = w5/lr;}
  if (w1 == 34) {dbc[w2][w3] = w4;}

  // Import boundary conditions
  if (w1 == 41) {bgn[w2][w3] = w4;}
  if (w1 == 42) {bel[w2] = w3;}
  if (w1 == 43) {bc[w4][w2][w3] = w5;}

  // Import initial conditions
  if ((w1 > 50) && (w1 < 60)) {
    for (i = 1; i <=4; ++i) {
		n = ie[w2][i];// cout << n << endl;
      if (w1 == 51) {cn[n][1] = w5;}
	  if (w1 == 52) {tn[n][1] = w5;} //cout << w5 << endl;
      if (w1 == 52) {tn[n][2] = w5;}
      if (w1 == 53) {vps[n][1] = w5;}
      if (w1 == 54) {
        vps[n][2] = w5;
        vps[n][3] = 0.0;
      }
    }
  }
																											    	//calling initdt
  if (w1 == 60) {
    initdt(nnp, nel, nsrf, band1, nnb, ph, ec, fl, tn, vi, pr);
    band3 = 3*band1;

	//cout << initdt << endl;

  }

  // Variable property data
  if (w1 == 71) {
    csp[1][1] = w2;
    csp[w2+1][w3] = w5;
    }
  else if (w1 == 72) {
    css[1][1] = w2;
    css[w2+1][w3] = w5;
    }
  else if (w1 == 73) {
    csb[1][1] = w2;
    csb[w2+1][w3] = w5;
  }

  // Loop over time                                                                                                 //calling cntrl
  if (w1 == 0) {
	  fo = w1*df; //cout << fo << endl;
	cntrl(se, nnp, nel, nsrf, band3, df, fo, ph, ec, fl, cn, tn, z, vps);
  }

  // Output data
  addv = 0.0;
  if (w1 == 100) {
	if (w2 == 1) addv = tn[w3][1]; //cout << tn[w3][1] << endl;
    if (w2 == 2) addv = cn[w3][1];
    if (w2 == 3) addv = (1.0F -fl[w3][1])*prs +fl[w3][1]*prl;
    if (w2 == 4) addv = vps[w3][1];
    if (w2 == 5) addv = vps[w3][2];
    if (w2 == 6) addv = 0.0F;
    if (w2 == 7) addv = vps[w3][3];
    if (w2 == 8) addv = vps[w3][10];
    if (w2 == 9) addv = fl[w3][1];
    if (w2 == 10) {
      cdcl(nnp, nsrf, um, cd, cl);                                                                               // calling cdcl
      if (w3 == 1) addv = cd;
      if (w3 == 2) addv = cl;
      if (w3 == 3) addv = um;
    }
  }

  // EXE output test data
  if (w1 == 200) addv = tn[w2][1];
  if (w1 == 201) addv = tn[w2][1];
  return (addv);
}