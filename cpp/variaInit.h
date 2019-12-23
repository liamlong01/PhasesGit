// -------------------------------------------------
// C++ PHASES - Copyright by G. F. Naterer
// Coding format: (i) // comments, (ii) column < 65,
// (ii) line continued -> one column indentation,
// (iv) start block -> two column indentation,
// (v) code optimization: set maximum speed options
// -------------------------------------------------

// This Header file contains all global variable declarations and
// any initializations if needed.
// Only the main file, phases.cpp needs
// to put -- #include "variaInit.h" -- at the top of the code.

int se; // there are 3 main solution blocks in cntrl, C(concentration equation), T (temperature) and U/V/P (velocities/ pressures),
		// se indicates which block is solved, numbered 1,2,3 or 4 (all together), Options: 
		// se = 1(C), 2(T), 3(UVP), 12(CT), 13(CUVP), 23(TUVP), 4(CTUVP)
int ao; // an addition capability of droplet flows was added later, ao = 2 means that the droplet flow solver is added, use ao = 1 (without droplets)
int tk; // number of times cycling through the temperature equation loop until a convergence tolerance is reached 
int tvk; // same as tk, but instead the number of time through the temperature/velocity loop 
int tstp; // number of time steps
int band3; // bandwidth if the matrix storing coefficients 
double df; // size of non-dimensional timestep, Fourier number
double prs; // density of a solid
double prl; // density of a liquid 
double pks; // conductivity of a solid  
double pkl; // conductivity of a liquid 
double pcs; // concentration of the soluidus 
double pcl; // concentration of the liquidus 
double pds0; // diffusivity of the solid 
double pdl0; // diffusivity of the liquid  
double phl; // specific heat of the liquid 
double phs; // specific heat of the solid 
double pl; // latent heat of fusion 
double rat; // thermal Rayleigh number 
double ras; // solutal Rayleigh number 
double tol; // convergence tolerance 
double lr; // reference length scale 
double vr; // reference v-velocity
double ur; // reference u-velocity 
double wr; // reference w-velocity 
double vsc; // viscosity 
double tsol; // temperature of the solidus 
double tmlt; // melting point temperature 
double tmin; // minimum temperature used to non-dimensionalize the temperature 
double tmax; // maximum temperature used to non-dimensionalize the temperature 
double dtr; // reference temperature, difference between max and min temperature 
double bt; // thermal expansion coefficent 
double bs; // solutal expansion coefficient 
double fo;// Fourier number, timestep the number of total steps 
double mu_e;
int ie[nelm][npe1]; // local to global node mapping
int bgn[nsrfm][3]; // boundary global node numbers
int bel[nsrfm]; // boundary element numbers
int nnb[nnpm][11]; // neighbouring nodes array
int ph[nnpm][3]; // phase number (time: n+1, n)
int dbc[nym][7]; // doubly connected domain nodes
double x[nnpm][3]; // global x and y coordinates
double s[4][3]; // local coordinates within element
double ds[npe1][5]; // elemental lengths, areas
double dn[npe1][5]; // shape functions, derivatives
double ic[nelm][3][3][npe1][2*npe1];
double bc[2*nsrfm][5][4]; // boundary conditions
double css[11][7]; // customized source terms
double csb[4*nym][8]; // customized boundary conditions
double csp[11][7]; // customized properties
double pr[nelm][npe1]; // Pr: molecular + turbulent parts
double ec[nnpm][7]; // equation of state (time: n+1, n)
double ivo[npe1]; // miscellaneous storage array
double nn[npe1][npe1]; // shape function interpolation
double cn[nnpm][6]; // C (t: n+1, n, iteration m, droplets)
double tn[nnpm][4]; // T (time: n+1, n, iteration m)
double fl[nnpm][3]; // liquid fraction (time: n+1, n)
double vi[nelm][npe1][7]; // main flow, droplet velocities
double c[3*nnpm][6*nym+11]; // global stiffness matrix
double r[3*nnpm]; // global right hand side
double z[3*nnpm]; // global solution array
double vps[nnpm][14]; // 1 - 3: U-V-P (time: n+1)
// 4 - 9: U-V (time: n, iteration m), s (entropy), area (CV)
// 10 - 13: entropy based diffusivity, smoothing coeffieicnts

// NON-DIMENSIONAL VARIABLE DEFINITIONS
//
// x[i][1] = x/lr, x[i][2] = y/lr (lr = reference length [m])
//vps[i][1] = u/v0 (v0 = ur = vr*al/lr is reference velocity)
//vps[i][2] = v/v0 (where diffusivity is al = pkl/(prl*phl))
// df = al*dt/(lr*lr) (where dt = time step size [s])
// tn[i][1] = (T-tmin)/(tmax-tmin) (where T = temperature [K])
// vps[i][3] = p/(prl*v0*v0) (where p is pressure [Pa])

