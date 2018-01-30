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

//const int nnpe = 4; // number of nodal points per element
//const int npe1 = nnpe + 1; // C++ arrays start at 0, add one
//const int nnpm = 2601; // maximum number of nodal points
//const int nelm = 2500; // maximum number of elements
//const int nsrfm = 320; // maximum number of boundary surfaces
//const int nym = 51; // maximum number of rows
//const double k0 = 1.0E+06; // permeability parameter
//const double eps = 0.001; // energy equation of state
//const int tph = 0; 
//const int fn = 2; 
//const int stol = 99;

//extern int se, ao, tk, tvk, tstp, band3;
//extern double df, prs, prl, pks, pkl, pcs, pcl;
//extern double pds0, pdl0, phl, phs, pl, rat, ras;
//extern double tol, lr, vr, ur, wr, vsc, tsol, tmlt;
//extern double tmin, tmax, dtr, bt, bs, fo;
//extern int ie[nelm][npe1]; // local to global node mapping
//extern int bgn[nsrfm][3]; // boundary global node numbers
//extern int bel[nsrfm]; // boundary element numbers
//extern int nnb[nnpm][11]; // neighbouring nodes array
//extern int ph[nnpm][3]; // phase number (time: n+1, n)
//extern int dbc[nym][7]; // doubly connected domain nodes
//extern double x[nnpm][3]; // global x and y coordinates
//extern double s[4][3]; // local coordinates within element
//extern double ds[npe1][5]; // elemental lengths, areas
//extern double dn[npe1][5]; // shape functions, derivatives
//extern double ic[nelm][3][3][npe1][2*npe1];
//extern double bc[2*nsrfm][5][4]; // boundary conditions
//extern double css[11][7]; // customized source terms
//extern double csb[4*nym][8]; // customiuzed boundary conditions
//extern double csp[11][7]; // customized properties
//extern double pr[nelm][npe1]; // Pr: molecular + turbulent parts
//extern double ec[nnpm][7]; // equation of state (time: n+1, n)
//extern double ivo[npe1]; // miscellaneous storage array
//extern double nn[npe1][npe1]; // shape function interpolation
//extern double cn[nnpm][6]; // C (t: n+1, n, iteration m, droplets)
//extern double tn[nnpm][4]; // T (time: n+1, n, iteration m)
//extern double fl[nnpm][3]; // liquid fraction (time: n+1, n)
//extern double vi[nelm][npe1][7]; // main flow, droplet velocities
//extern double c[3*nnpm][6*nym+11]; // global stiffness matrix
//extern double r[3*nnpm]; // global right hand side
//extern double z[3*nnpm]; // global solution array
//extern double vps[nnpm][14]; // 1 - 3: U-V-P (time: n+1)

//static const int o[12] = {0,4,1,2,3,4,1,2,2,1,3,4};
/*static const double lc[17][3] = {0.0,0.0,0.0, 0.0,0.0,0.5,
 0.0,-0.5,0.0, 0.0,0.0,-0.5, 0.0,0.5,0.0, 0.0,0.0,0.0,
 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,1.0,
 0.0,-1.0,0.0, 0.0,0.0,-1.0, 0.0,1.0,0.0, 0.0,0.5,0.5,
 0.0,-0.5,0.5, 0.0,-0.5,-0.5, 0.0,0.5,-0.5};*/
/*static const double nc[25][3] = {0.0,0.0,0.0, 0.0,0.5,1.0,
 0.0,-0.5,1.0, 0.0,-1.0,0.5, 0.0,-1.0,-0.5, 0.0,-0.5,-1.0,
 0.0,0.5,-1.0, 0.0,1.0,-0.5, 0.0,1.0,0.5, 0.0,1.0,1.0,
 0.0,0.0,1.0, 0.0,-1.0,1.0, 0.0,-1.0,0.0, 0.0,-1.0,-1.0,
 0.0,0.0,-1.0, 0.0,1.0,-1.0, 0.0,1.0,0.0, 0.0,0.0,1.0,
 0.0,-1.0,1.0, 0.0,-1.0,0.0, 0.0,-1.0,-1.0, 0.0,0.0,-1.0,
 0.0,1.0,-1.0, 0.0,1.0,0.0, 0.0,1.0,1.0};*/

///////////////////////////////////////////////////////////////////////////////      calls no functions 
