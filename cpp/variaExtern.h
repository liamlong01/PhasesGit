// -------------------------------------------------
// C++ PHASES - Copyright by G. F. Naterer
// Coding format: (i) // comments, (ii) column < 65,
// (ii) line continued -> one column indentation,
// (iv) start block -> two column indentation,
// (v) code optimization: set maximum speed options
// -------------------------------------------------

// This Header file contains all global variables external declarations.
// Any cpp file, besides the main phases.cpp file, needs
// to put -- #include "variaExtern.h" -- at the top of the code.

extern int se;
extern int ao;
extern int tk;
extern int tvk;
extern int tstp;
extern int band3;
extern double df;
extern double prs;
extern double prl;
extern double pks;
extern double pkl;
extern double pcs;
extern double pcl;
extern double pds0;
extern double pdl0;
extern double phl;
extern double phs;
extern double pl;
extern double rat;
extern double ras;
extern double tol;
extern double lr;
extern double vr;
extern double ur;
extern double wr;
extern double vsc;
extern double tsol;
extern double tmlt;
extern double mu_e;
extern double tmin;
extern double tmax;
extern double dtr;
extern double bt;
extern double bs;
extern double fo;
extern int calcppr;
extern int ie[nelm][npe1]; // local to global node mapping
extern int bgn[nsrfm][3]; // boundary global node numbers
extern int bel[nsrfm]; // boundary element numbers
extern int nnb[nnpm][11]; // neighbouring nodes array
extern int ph[nnpm][3]; // phase number (time: n+1, n)
extern int dbc[nym][7]; // doubly connected domain nodes
extern double x[nnpm][3]; // global x and y coordinates
extern double s[4][3]; // local coordinates within element
extern double ds[npe1][5]; // elemental lengths, areas
extern double dn[npe1][5]; // shape functions, derivatives
extern double ic[nelm][3][3][npe1][2*npe1];
extern double bc[2*nsrfm][5][4]; // boundary conditions
extern double css[11][7]; // customized source terms
extern double csb[4*nym][8]; // customiuzed boundary conditions
extern double csp[11][7]; // customized properties
extern double pr[nelm][npe1]; // Pr: molecular + turbulent parts
extern double ec[nnpm][7]; // equation of state (time: n+1, n)
extern double ivo[npe1]; // miscellaneous storage array
extern double nn[npe1][npe1]; // shape function interpolation
extern double cn[nnpm][6]; // C (t: n+1, n, iteration m, droplets)
extern double tn[nnpm][4]; // T (time: n+1, n, iteration m)
extern double fl[nnpm][3]; // liquid fraction (time: n+1, n)
extern double vi[nelm][npe1][7]; // main flow, droplet velocities
extern double c[3*nnpm][6*nym+11]; // global stiffness matrix
extern double r[3*nnpm]; // global right hand side
extern double z[3*nnpm]; // global solution array
extern double vps[nnpm][14]; // 1 - 3: U-V-P (time: n+1)

// 4 - 9: U-V (time: n, iteration m), s (entropy), area (CV)
// 10 - 13: entropy based diffusivity, smoothing coeffieicnts

// NON-DIMENSIONAL VARIABLE DEFINITIONS
//
// x[i][1] = x/lr, x[i][2] = y/lr (lr = reference length [m])
// vps[i][1] = u/v0 (v0 = ur = vr*al/lr is reference velocity)
// vps[i][2] = v/v0 (where diffusivity is al = pkl/(prl*phl))
// df = al*dt/(lr*lr) (where dt = time step size [s])
// tn[i][1] = (T-tmin)/(tmax-tmin) (where T = temperature [K])
// vps[i][3] = p/(prl*v0*v0) (where p is pressure [Pa])
