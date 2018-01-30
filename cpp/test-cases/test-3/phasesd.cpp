
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

const int nnpe = 4; // number of nodal points per element
const int npe1 = nnpe + 1; // C++ arrays start at 0, add one
const int nnpm = 2601; // maximum number of nodal points
const int nelm = 2500; // maximum number of elements
const int nsrfm = 320; // maximum number of boundary surfaces
const int nym = 51; // maximum number of rows
const double k0 = 1.0E+06; // permeability parameter
const double eps = 0.001; // energy equation of state
const int tph = 0; const int fn = 100; const int stol = 99;

static int se, ao, tk, tvk, tstp, band3;
static double df, prs, prl, pks, pkl, pcs, pcl;
static double pds0, pdl0, phl, phs, pl, rat, ras;
static double tol, lr, vr, ur, wr, vsc, tsol, tmlt;
static double tmin, tmax, dtr, bt, bs, fo;
static int ie[nelm][npe1]; // local to global node mapping
static int bgn[nsrfm][3]; // boundary global node numbers
static int bel[nsrfm]; // boundary element numbers
static int nnb[nnpm][11]; // neighbouring nodes array
static int ph[nnpm][3]; // phase number (time: n+1, n)
static int dbc[nym][7]; // doubly connected domain nodes
static double x[nnpm][3]; // global x and y coordinates
static double s[4][3]; // local coordinates within element
static double ds[npe1][5]; // elemental lengths, areas
static double dn[npe1][5]; // shape functions, derivatives
static double ic[nelm][3][3][npe1][2*npe1];
static double bc[2*nsrfm][5][4]; // boundary conditions
static double css[11][7]; // customized source terms
static double csb[4*nym][8]; // customiuzed boundary conditions
static double csp[11][7]; // customized properties
static double pr[nelm][npe1]; // Pr: molecular + turbulent parts
static double ec[nnpm][7]; // equation of state (time: n+1, n)
static double ivo[npe1]; // miscellaneous storage array
static double nn[npe1][npe1]; // shape function interpolation
static double cn[nnpm][6]; // C (t: n+1, n, iteration m, droplets)
static double tn[nnpm][4]; // T (time: n+1, n, iteration m)
static double fl[nnpm][3]; // liquid fraction (time: n+1, n)
static double vi[nelm][npe1][7]; // main flow, droplet velocities
static double c[3*nnpm][6*nym+11]; // global stiffness matrix
static double r[3*nnpm]; // global right hand side
static double z[3*nnpm]; // global solution array
static double vps[nnpm][14]; // 1 - 3: U-V-P (time: n+1)
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

static const int o[12] = {0,4,1,2,3,4,1,2,2,1,3,4};
static const double lc[17][3] = {0.0,0.0,0.0, 0.0,0.0,0.5,
 0.0,-0.5,0.0, 0.0,0.0,-0.5, 0.0,0.5,0.0, 0.0,0.0,0.0,
 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,1.0,
 0.0,-1.0,0.0, 0.0,0.0,-1.0, 0.0,1.0,0.0, 0.0,0.5,0.5,
 0.0,-0.5,0.5, 0.0,-0.5,-0.5, 0.0,0.5,-0.5};
static const double nc[25][3] = {0.0,0.0,0.0, 0.0,0.5,1.0,
 0.0,-0.5,1.0, 0.0,-1.0,0.5, 0.0,-1.0,-0.5, 0.0,-0.5,-1.0,
 0.0,0.5,-1.0, 0.0,1.0,-0.5, 0.0,1.0,0.5, 0.0,1.0,1.0,
 0.0,0.0,1.0, 0.0,-1.0,1.0, 0.0,-1.0,0.0, 0.0,-1.0,-1.0,
 0.0,0.0,-1.0, 0.0,1.0,-1.0, 0.0,1.0,0.0, 0.0,0.0,1.0,
 0.0,-1.0,1.0, 0.0,-1.0,0.0, 0.0,-1.0,-1.0, 0.0,0.0,-1.0,
 0.0,1.0,-1.0, 0.0,1.0,0.0, 0.0,1.0,1.0};

void assmbl(int v, int p, int m, int nel, int nelm,
 int band3, double df, double fo, double nn[npe1][npe1],
 double c[3*nnpm][6*nym+11], double r[3*nnpm]);
void bndry(int v, int p, int m, int nsrf, double fo,
 double bc[2*nsrfm][5][4], double c[3*nnpm][6*nym+11],
 double r[3*nnpm]);
void cdcl(int nnp, int nsrf, double &um, double &cd, double &cl);
void cntrl(int se, int nnp, int nel, int nsrf,
 int band3, double df, double fo, int ph[nnpm][3],
 double ec[nnpm][7], double fl[nnpm][3], double cn[nnpm][6],
 double tn[nnpm][4], double z[3*nnpm], double vps[nnpm][14]);
void coeffc(int e, double cph[npe1][3], double rd[npe1][5]);
void custom(int n, int m, double &tm, double &ci, double &ti,
 double csp[11][7], double css[11][7], double csb[4*nym][8],
 double &cval);
void dsmooth(int nel, int nnp, double vps[nnpm][14]);
void entropy(double df, int nel, int nnp, int nsrf,
 double vps[nnpm][14]);
void filedata(int fn);
void forbak(int m, double c[3*nnpm][6*nym+11], double z[3*nnpm],
 double r[3*nnpm], int nym, int &neq, int nnpm, int band3);
void genul(int m, double c[3*nnpm][6*nym+11], int nym,
 int &neq, int nnpm, int band3);
void initdt(int nnp, int nel, int nsrf, int &band1,
 int nnb[nnpm][11], int ph[nnpm][3], double ec[nnpm][7],
 double fl[nnpm][3], double tn[nnpm][4],
 double vi[nelm][npe1][7], double pr[nelm][npe1]);
void ipoint(int eqn, int e, double dfs, double xip[npe1][3],
 double ic[nelm][3][3][npe1][2*npe1]);
void nullmx(double r[3*nnpm], double z[3*nnpm],
 double c[3*nnpm][6*nym+11]);
void phase(int nnp, int ph[nnpm][3], double fl[nnpm][3],
 double ec[nnpm][7], int &tr);
void phase3(int nel, int nnp, int &cmf, int &cmx, int &cr,
 int nsrf, double df, int ph[nnpm][3], double fl[nnpm][3],
 double cn[nnpm][6], double ec[nnpm][7]);
void shape(int e, int nnpm, double &dx, double &dy,
 double s[4][3], double dn[npe1][5], double ivo[npe1],
 double &jac, double &dqdx, double &dqdy);
void stiffc(int e, double df, double fo, double xip[npe1][3],
 double flux[4][npe1][npe1], double sc[npe1],
 double aq[4][npe1][npe1], double rq[npe1]);
void stiffdc(int e, double df, double xip[npe1][3],
 double flux[4][npe1][npe1], double aq[4][npe1][npe1],
 double rq[npe1]);
void stiffdv(int e, double vi[nelm][npe1][7]);
void stiffp(int e, double flux[4][npe1][npe1], double sc[npe1],
 double aq[4][npe1][npe1], double rq[npe1]);
void stifft(int e, double df, double fo, double cip[npe1],
 double tip[npe1], double flux[4][npe1][npe1], double sc[npe1],
 double aq[4][npe1][npe1], double rq[npe1]);
void stiffu(int e, double df, double xip[npe1][3],
 double flux[4][npe1][npe1], double sc[npe1],
 double aq[4][npe1][npe1], double rq[npe1]);
void stiffv(int e, double df, double flux[4][npe1][npe1],
 double sc[npe1], double aq[4][npe1][npe1], double rq[npe1]);
void uvcor(int nel, double vi[nelm][npe1][7]);

// Exported EXE -> console
// 1. New Project -> Win32 Project -> console application
// 2. use: 'void main()', both 'double Adda()' declarations
// 4. remove: '#include <windows.h>', 'double _declspec ...'
void main()
{
  filedata(fn);
}
double Adda(int w1, int w2, int w3, int w4, double w5, int nnp,
 int nel, int nx, int ny, int nsrf);
double Adda(int w1, int w2, int w3, int w4, double w5, int nnp,
 int nel, int nx, int ny, int nsrf)

// Exported DLL -> Visual Basic GUI
// 1. New Project -> Win32 Project -> DLL
// 2. use: '#include <windows.h>', 'double _declspec ...'
// 3. use: phasesb.def source file: LIBRARY phasesb
//                                  EXPORTS Adda
// 1. remove: 'void main()', previous 'double Adda ...'
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
    if (w1 == 6) {df = w5;}
    if (w1 == 7) {lr = w5;}
    if (w1 == 8) {wr = w5;}
    if (w1 == 9) {tol = w5;}
    if (w1 == 10) {ur = w5;}
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
      tmlt = (tmlt -tmin)/(tmax -tmin +1.0E-32);
      tsol = (tsol -tmin)/(tmax -tmin +1.0E-32);
      rat = 9.8*bt*(tmax -tmin)*lr*lr*lr*prl*phl/(vsc*pkl
       +1.0E-32);
      ras = 9.8*bs*(pcl -pcs)*lr*lr*lr*prl*phl/(vsc*pkl
       +1.0E-32);
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
      n = ie[w2][i];
      if (w1 == 51) {cn[n][1] = w5;}
      if (w1 == 52) {tn[n][1] = w5;}
      if (w1 == 52) {tn[n][2] = w5;}
      if (w1 == 53) {vps[n][1] = w5;}
      if (w1 == 54) {
        vps[n][2] = w5;
        vps[n][3] = 0.0;
      }
    }
  }
  if (w1 == 60) {
    initdt(nnp, nel, nsrf, band1, nnb, ph, ec, fl, tn, vi, pr);
    band3 = 3*band1;
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

  // Loop over time
  if (w1 == 0) {
		fo = w2*df;
    cntrl(se, nnp, nel, nsrf, band3, df, fo, ph, ec, fl,
     cn, tn, z, vps);
	}

  // Output data
	addv = 0.0;
  if (w1 == 100) {
    if (w2 == 1) addv = tn[w3][1];
    if (w2 == 2) addv = cn[w3][1];
    if (w2 == 3) addv = (1.0F -fl[w3][1])*prs +fl[w3][1]*prl;
    if (w2 == 4) addv = vps[w3][1];
    if (w2 == 5) addv = vps[w3][2];
    if (w2 == 6) addv = 0.0F;
    if (w2 == 7) addv = vps[w3][3];
    if (w2 == 8) addv = vps[w3][10];
    if (w2 == 9) addv = fl[w3][1];
    if (w2 == 10) {
      cdcl(nnp, nsrf, um, cd, cl);
      if (w3 == 1) addv = cd;
      if (w3 == 2) addv = cl;
      if (w3 == 3) addv = um;
    }
  }

  // EXE output test data
  if (w1 == 200) addv = cn[w2][1];
  if (w1 == 201) addv = tn[w2][1];
  return (addv);
}

// Assembly of elemental stiffness equations
void assmbl(int v, int p, int m, int nel, int nelm,
 int band3, double df, double fo, double nn[npe1][npe1],
 double c[3*nnpm][6*nym+11], double r[3*nnpm])
{
  int e, i, j, n, iw, jw, im, iband, isemi;
  double jac, dx, dy, ss, tt, dqdx, dqdy;
  double xip[npe1][3], cip[npe1], tip[npe1];
  static double aq[4][npe1][npe1], rq[npe1];
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
      shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);
      ds[j][3] = fabs(jac);
      ds[j][4] = fabs(jac);
      s[1][1] = lc[j][1];
      s[1][2] = lc[j][2];
      shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);
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
    if (v == 1) {stiffc(e, df, fo, xip, flux, sc, aq, rq);}
    if (v == 2) {stifft(e, df, fo, cip, tip, flux, sc, aq, rq);}
    if (v == 3) {stiffu(e, df, xip, flux, sc, aq, rq);}
    if (v == 4) {stiffv(e, df, flux, sc, aq, rq);}
    if (v == 5) {stiffp(e, flux, sc, aq, rq);}
    if (v == 6) {stiffdv(e, vi);}
    if (v == 7) {stiffdc(e, df, xip, flux, aq, rq);}

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

// Boundary conditions
void bndry(int v, int p, int m, int nsrf, double fo,
 double bc[2*nsrfm][5][4], double c[3*nnpm][6*nym+11],
 double r[3*nnpm])
{
  int e, f, i, j, n, np, s, n1, n2, ni, p1, p2, p3;
  int n3, mbc, b1, b2, bi, isemi, iband, oi;
  double dx, dy, ci, ti, tm, tp, bval;

  // Boundary coefficients (customized)
  bval = 0.0;
  tm = fo*lr*lr*prl*phl/pkl;
  for (j = 1; j <= csb[1][1]; ++j) {
    i = int(csb[j+1][1]);
    ci = cn[bgn[i][1]][1];
    ti = tmin +tn[bgn[i][1]][1]*fabs(tmax -tmin);
    n = int(csb[j+1][2]) -10;
    f = int(csb[j+1][3]) -30;
    b1 = int(2.0*(i -1) +1);
    b2 = b1 +1;
    if ((f >= 1) && (f <= 3)) {
      custom(3, j, tm, ci, ti, csp, css, csb, bval);
      bc[b1][n][f] = bval;
      bc[b2][n][f] = bval;
      if (f == 1) {
        bc[b1][2][1] = bval/lr;
        bc[b2][2][1] = bval/lr;
      }
      else if (f == 3) {
        bc[b1][2][3] = (bval -bc[b1][2][2]*tmin)/fabs(tmax -tmin);
        bc[b2][2][3] = (bval -bc[b2][2][2]*tmin)/fabs(tmax -tmin);
      }
    }
  }

  // Loop over boundary surfaces
  iband = m*(band3 +m -3)/3;
  isemi = (iband +1)/2;
  p1 = 0;
  p2 = 0;
  p3 = 0;
  for (s = 1; s <= nsrf; ++s) {
    e = bel[s];
    b1 = int(2.0*(s -1) +1);
    b2 = b1 +1;
    n1 = p +m*(bgn[s][1] -1);
    n2 = p +m*(bgn[s][2] -1);
    dx = 0.5*(x[bgn[s][2]][1] -x[bgn[s][1]][1]);
    dy = 0.5*(x[bgn[s][2]][2] -x[bgn[s][1]][2]);

    // Loop over element nodes
    for (i = 1; i <= 2; ++i) {
      if (i == 1) {
        ni = n1;
        bi = b1;
      }
      else {
        ni = n2;
        bi = b2;
      }

      // Continuity equation
      if (p == 3) {
        mbc = 3; // mbc = 1: gradient spacified condition
        // mbc = 2: pressure specified at single boudnary
        // mbc = 3: pressure specified at multiple boundaries

        // Boundary mass equation closure
        if (mbc == 1) {
          n1 = 3 +m*(bgn[s][1] -1);
          n2 = 3 +m*(bgn[s][2] -1);
          n3 = m*(bgn[s][2] -bgn[s][1]);

          // Gradient specified condition
          if ((bc[b1][3][2] == 0.0) || (bc[b1][4][2] == 0.0)) {
            c[n1][isemi-2] = c[n1][isemi-2] -dy*0.75;
            c[n1][isemi-2+n3] = c[n1][isemi-2+n3] -dy*0.25;
            c[n1][isemi-1] = c[n1][isemi-1] +dx*0.75;
            c[n1][isemi-1+n3] = c[n1][isemi-1+n3] +dx*0.25;
            c[n2][isemi-2] = c[n2][isemi-2] -dy*0.75;
            c[n2][isemi-2-n3] = c[n2][isemi-2-n3] -dy*0.25;
            c[n2][isemi-1] = c[n2][isemi-1] +dx*0.75;
            c[n2][isemi-1-n3] = c[n2][isemi-1-n3] +dx*0.25;
          }
 
          // Velocity specified condition
          else {
            r[n1] = r[n1] +dy*bc[b1][3][3]/(bc[b1][3][2]
             +1.0E-16) -dx*bc[b1][4][3]/(bc[b1][4][2] +1.0E-16);
            r[n2] = r[n2] +dy*bc[b2][3][3]/(bc[b2][3][2]
             +1.0E-16) -dx*bc[b2][4][3]/(bc[b2][4][2] +1.0E-16);
          }
        }

        // Pressure based conditions
        else if (mbc >= 2) {
          np = 0;
          for (f = 3; f <= 4; ++f) {
            for (j = 0; j <= 1; ++j) {
              bi = int(2.0*(s -1) +1) +j;
              if ((bc[bi][f][1] == 0) && (bc[bi][f][3] != 0))
                {np = 1;}
              if ((bc[bi][f][2] == 0) && (bc[bi][f][3] == 0))
                {np = 2;}
            }
          }

          // Reference pressure specification
          if (np == 1) {p1 = p1++;}
          if (np != 1) {p2 = p2++;}
          if ((p1 > 1) && (p2 > 1) && (np != 1)) {p3 = p3++;}
          if (((mbc == 2) && (np == 1) && (p3 <= 1))
            || ((mbc == 3) && (np == 1))) {
            for (j = 1; j <= iband; ++j) {
              c[ni][j] = 0.0;
            }
            c[ni][isemi] = 1.0;
            r[ni] = 0.0;
          }

          // Zero gradient condition
          if (np == 2) {
            oi = int((bc[bi][3][1] -0.99999)*1.0E+03);
            f = isemi +m*(ie[e][o[oi+1]] -ie[e][o[oi+2]]);
            if (f > 0) {
              for (j = 1; j <= iband; ++j) {
                c[ni][j] = 0.0;
              }
              c[ni][isemi] = 1.0;
              c[ni][f] = -1.0;
              r[ni] = 0.0;
            }
          }
        }
      }

      // Dirichlet boundary condition
      else if (bc[bi][v][1] == 0.0) {
        for (j = 1; j <= iband; ++j) {
          c[ni][j] = 0.0;
        }
        c[ni][isemi] = 1.0;
        r[ni] = bc[bi][v][3]/(bc[bi][v][2] +1.0E-32);
      }

      // Neumann boundary condition
      else if ((bc[bi][v][1] > 1.0) && (bc[bi][v][1] < 1.01)
       && (bc[bi][v][2] == 0.0)) {
        oi = int((bc[bi][v][1] -0.99999)*1.0E+03);
        for (j = 1; j <= iband; ++j) {
          c[ni][j] = 0.0;
        }
        c[ni][isemi] = 1.0;
        c[ni][isemi+m*(ie[e][o[oi+1]] -ie[e][o[oi+2]])] = -1.0;
        r[ni] = sqrt(dx*dx +dy*dy)*bc[bi][v][3]/(bc[bi][v][1]
         +1.0E-32);
      }
    
      // Robin boundary condition
      else if ((bc[bi][v][1] != 0) && (bc[bi][v][2] != 0)) {
        tp = r[ni]/(c[ni][isemi] + 1.0E-32);
        for (j = 1; j <= iband; ++j) {
          if (j != isemi) {
            tp = tp-c[ni][j]*tn[j][1]/(c[ni][isemi]+1.0E-32);
          }
        }
        bc[bi][v][3] = bc[bi][v][3]/(bc[bi][v][1] +1.0E-32);
        bc[bi][v][2] = bc[bi][v][2]/(bc[bi][v][1] +1.0E-32);
        bc[bi][v][1] = 1.0;
        r[ni] = r[ni] +sqrt(dx*dx +dy*dy)*(bc[bi][v][3]
         -bc[bi][v][2]*tp);
      }
    }
  }

  // Common edge of doubly connected region
  for (i = 1; i <= nym; ++i) {
    if (dbc[i][1] > 0) {
      for (n = 1; n <= 4; ++n) {
        e = dbc[i][1];
        if (n <= 2) {s = dbc[i][n+1];}
        if (n >= 3) {s = dbc[i][n+2];}
        ni = p +m*(s -1);
        for (j = 1; j <= iband; ++j) {
          if (j != isemi) {c[ni][j] = 0.0;}
        }
        c[ni][isemi] = 1.0;
        if (v == 1) {r[ni] = cn[s][1];}
        if (v == 2) {r[ni] = tn[s][1];}
        if (v == 4) {r[ni] = vps[s][1];}
        if (v == 5) {r[ni] = vps[s][2];}
      }
    }
  }
}

// Main control sequence of inter-equation iterations
void cntrl(int se, int nnp, int nel, int nsrf,
 int band3, double df, double fo, int ph[nnpm][3],
 double ec[nnpm][7], double fl[nnpm][3], double cn[nnpm][6],
 double tn[nnpm][4], double z[3*nnpm], double vps[nnpm][14])
{
  // Options: se = 1(C), 2(T), 3(UVP), 12(CT)
  // Options: se = 13(CUVP), 23(TUVP), 4(CTUVP)
  int i, j, n, cr, tr, tph, cmf, cmx, vk, np3;
  double cdif, cphs, cphl, tphs, tphl;
  double vdif, tdif, ste, rlx, hbar;

  // Store field variables
  ste = phl*dtr/(pl +1.0E-32);
  np3 = nnp*3;
	tdif = 0.0;
  for (i = 1; i <= nnp; ++i) {
    cn[i][2] = cn[i][1];
    cn[i][3] = cn[i][1];
    cn[i][4] = cn[i][1];
    tn[i][2] = tn[i][1];
    tn[i][3] = tn[i][1];
    fl[i][2] = fl[i][1];
    ph[i][2] = ph[i][1];
    vps[i][4] = vps[i][1];
    vps[i][5] = vps[i][2];
    vps[i][6] = vps[i][1];
    vps[i][7] = vps[i][2];
    ec[i][4] = ec[i][1];
    ec[i][5] = ec[i][2];
    ec[i][6] = ec[i][3];
  }

  // Start C - T - V iteration sequence
  for (vk = 1; vk <= tvk; ++vk) {

    // Start C equation iterations
    if ((se == 1) || (se == 12) || (se == 13) || (se == 4)) {
      for (i = 1; i <= nnp; ++i) {
        cn[i][2] = cn[i][3];
        cn[i][4] = 0.0;
        fl[i][1] = fl[i][2];
      }
      cmx = 2;
      cr = 0;
      for (cmf = 1; cmf <= tk; ++cmf) {
        cdif = 0.0;

        // C equation: solid / liquid phase change
        if (ao != 2) {
          nullmx(r, z, c);
          assmbl(1, 1, 1, nel, nelm, band3, df, fo, nn, c, r);
          bndry(1, 1, 1, nsrf, fo, bc, c, r);
          genul(1, c, nym, nnp, nnpm, band3);
          forbak(1, c, z, r, nym, nnp, nnpm, band3);

          // Solution residual
          for (i = 1; i <= nnp; ++i) {
            cn[i][1] = z[i];
            cdif = cdif +fabs(cn[i][1] -cn[i][3])/nnp;
            cn[i][3] = cn[i][1];
          }
        }

        // C equation: droplet flow
        else if (ao == 2) {
          for (i = 1; i <= nnp; ++i) {
            ph[i][2] = ph[i][1];
            cn[i][5] = 0.0;
          }

          // External flow: droplet momentum / mass equations
          assmbl(6, 1, 1, nel, nelm, band3, df, fo, nn, c, r);
          nullmx(r, z, c);
          assmbl(7, 1, 1, nel, nelm, band3, df, fo, nn, c, r);
          bndry(1, 1, 1, nsrf, fo, bc, c, r);
          genul(1, c, nym, nnp, nnpm, band3);
          forbak(1, c, z, r, nym, nnp, nnpm, band3);
          for (j = 1; j <= nnp; ++j) {
            cn[j][1] = z[j];
          }

          // Surface: three-phase heat, mass / momentum balances
          phase3(nel, nnp, cmf, cmx, cr, nsrf, df, ph, fl, cn, ec);
          if (cmf >= cmx) {cr = 0;}
        }

        // Quit early for small residual
        if ((cr == 0) && (cdif < tol)) {goto cnt3000;}
      }
    }

    // Start temperature - phase iteration algorithm
    cnt3000: if ((se == 2) || (se == 12) || (se == 23)
     || (se == 4)) {
      for (tph = 1; tph <= tk; ++tph) {
        for (i = 1; i <= nnp; ++i) {
          tn[i][3] = tn[i][1];
          ph[i][2] = ph[i][1];
        }

        // Energy equation: single / multiphase flows
        nullmx(r, z, c);
        assmbl(2, 1, 1, nel, nelm, band3, df, fo, nn, c, r);
        bndry(2, 1, 1, nsrf, fo, bc, c, r);
				genul(1, c, nym, nnp, nnpm, band3);
        forbak(1, c, z, r, nym, nnp, nnpm, band3);
      
        // Temperature solution residual
        tdif = 0.0;
        tr = 0;
        for (i = 1; i <= nnp; ++i) {
          tn[i][1] = z[i];
          tdif = tdif +fabs(tn[i][1] -tn[i][3])/nnp;
        }

        // Apply phase interface motion rules
        if (ao != 2) {phase(nnp, ph, fl, ec, tr);}
        
        // Quit early if tentative phase matches new phase
        if ((tr == 0) || (tdif < tol)) {goto cnt6000;}
      }
    }

    // Phase - temperature consistency
    cnt6000: if (ao != 2) {
      for (i = 1; i <= nnp; ++i) {
        tphs = tmlt -eps -cn[i][1]*(tmlt -eps -tsol)
         /(pcs +1.0E-32);
        tphl = tmlt +eps -cn[i][1]*(tmlt +eps -tsol)
         /(pcl +1.0E-32);
        cphs = pcs*(tmlt -eps -tn[i][1])/(tmlt -eps -tsol);
        cphl = pcl*(tmlt +eps -tn[i][1])/(tmlt +eps -tsol);
        hbar = 0.5*(phs+phl);

        // Solid phase reference states
        if (tn[i][1] <= tphs) {
          ph[i][1] = 1;
          fl[i][1] = 0.0;
          ec[i][1] = 0.0;
          ec[i][2] = phs/(phl +1.0E-32);
          ec[i][3] = 0.0;
        }

        // Liquid phase reference states
        else if (tn[i][1] >= tphl) {
          ph[i][1] = 3;
          fl[i][1] = 1.0;
          ec[i][1] = hbar*(tphl-tphs)/(phl +1.0E-32) 
           +(phs*tphs)/(phl +1.0E-32) +1.0/(ste +1.0E-32);
          ec[i][2] = phl/(phl +1.0E-32);
          ec[i][3] = tphl;
        }

        // Melt phase reference states
        else {
          ph[i][1] = 2;
          ec[i][1] = (phs*tphs)/(phl +1.0E-32);
          ec[i][2] = hbar/(phl +1.0E-32) +1.0/(ste*(tphl -tphs)
           +1.0E-32);
          ec[i][3] = tphs;
          fl[i][1] = (tn[i][1] -tphs)/(tphl -tphs);
        }
      }
    }

    // Start U - V - P solution iterations
    if ((se == 3) || (se == 13) || (se == 23) || (se == 4)) {
      for (n = 1; n <= tk; ++n) {

        // Element assembly: x momentum, y momentum, mass equations
        nullmx(r, z, c);
        assmbl(3, 1, 3, nel, nelm, band3, df, fo, nn, c, r);
        assmbl(4, 2, 3, nel, nelm, band3, df, fo, nn, c, r);
        assmbl(5, 3, 3, nel, nelm, band3, df, fo, nn, c, r);

        // Boundary conditions, U-V-P direct simultaneous solution
        bndry(3, 1, 3, nsrf, fo, bc, c, r);
        bndry(4, 2, 3, nsrf, fo, bc, c, r);
        bndry(4, 3, 3, nsrf, fo, bc, c, r);
				genul(3, c, nym, np3, nnpm, band3);
        forbak(3, c, z, r, nym, np3, nnpm, band3);
        uvcor(nel, vi);
        for (i = 1; i <= nnp; ++i) {
          vps[i][1] = z[3*(i -1) +1];
          vps[i][2] = z[3*(i -1) +2];
          vps[i][3] = z[3*(i -1) +3];
        }

        // Velocity residual
        vdif = 0.0;
        rlx = 0.5;
        for (i = 1; i <= nnp; ++i) {
          vdif = vdif +fabs(vps[i][6] -vps[i][1])/nnp;
          vps[i][6] = vps[i][1] +rlx*(vps[i][6] -vps[i][1]);
          vps[i][7] = vps[i][2] +rlx*(vps[i][7] -vps[i][2]);
        } 
        if (vdif <= tol) {goto cnt8200;}
      }
      
      // Quit early for small residual
      cnt8200: if ((tdif < tol) && (vdif < tol)) {goto cnt9600;}
    }
  }

  // Entropy extensions
  cnt9600: if (ao == 1) {
    entropy(df, nel, nnp, nsrf, vps);
    dsmooth(nel, nnp, vps);
  }
}

// LU decomposition into lower - upper triangular matrices
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
}

// Banded matrix solver by forward - backward substitution
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
}

// Initial conditions and problem parameters
void initdt(int nnp, int nel, int nsrf, int &band1,
 int nnb[nnpm][11], int ph[nnpm][3], double ec[nnpm][7],
 double fl[nnpm][3], double tn[nnpm][4],
 double vi[nelm][npe1][7], double pr[nelm][npe1])
{
  int e, i, j, k, n, m, iw, i1, i2;
	int dif, maxdif, j1, j2, j3, j4;
  double cphs, cphl, tphs, tphl, um, vm;
  double ste, hbar, u1, v1, vdn;

  // Start from initial conditions
  for (e = 1; e <= nel; ++e) {
    for (n = 1; n <= nnpe; ++n) {
      i = ie[e][n];
      vi[e][n][1] = vps[i][1];
      vi[e][n][2] = vps[i][2];
      vi[e][n][5] = vps[i][1];
      vi[e][n][6] = vps[i][2];
      pr[e][n] = (vsc*prl*phl)/(pkl +1.0E-32);
    }
  }

  // Initial phase distributions
  ste = phl*fabs(tmax-tmin)/(pl + 1.0E-32);
  hbar = 0.5*(phs+phl);
  for (i = 1; i <= nnp; ++i) {
    tphs = tmlt -eps -cn[i][1]*(tmlt -eps -tsol)/(pcs
     +1.0E-32);
    tphl = tmlt +eps -cn[i][1]*(tmlt +eps -tsol)/(pcl
     +1.0E-32);
    cphs = pcs*(tmlt -eps -tn[i][1])/(tmlt -eps -tsol);
    cphl = pcl*(tmlt +eps -tn[i][1])/(tmlt +eps -tsol);

    // Solid phase reference states
    if (tn[i][1] <= tphs) {
      ph[i][1] = 1;
      fl[i][1] = 0.0;
      ec[i][1] = 0.0;
      ec[i][2] = phs/(phl +1.0E-32);
      ec[i][3] = 0.0;
    }

    // Liquid phase reference states
    else if (tn[i][1] >= tphl) {
      ph[i][1] = 3;
      fl[i][1] = 1.0;
      ec[i][1] = hbar*(tphl-tphs)/(phl +1.0E-32)
       +(phs*tphs)/(phl +1.0E-32) +1.0/(ste +1.0E-32);
      ec[i][2] = phl/(phl +1.0E-32);
      ec[i][3] = tphl;
    }

    // Melt phase reference states
    else {
      ph[i][1] = 2;
      ec[i][1] = (phs*tphs)/(phl +1.0E-32);
      ec[i][2] = hbar/(phl +1.0E-32) +1.0/(ste*(tphl -tphs)
       +1.0E-32);
      ec[i][3] = tphs;
      fl[i][1] = (tn[i][1] -tphs)/(tphl -tphs);
    }
  }

  // Droplet velocities
  if (ao == 2) {
    um = 0.0;
    vm = 0.0;
    for (i = 1; i <= nsrf; ++i) {
      for (k = 1; k <= 2; ++k) {
        j = 2*(i -1) +k;
        u1 = sqrt(um*um +vm*vm);
        v1 = sqrt(pow(bc[j][3][3]/(bc[j][3][2] +1.0E-32), 2.0)
         +pow(bc[j][4][3]/(bc[j][4][2] +1.0E-32), 2.0));
        if ((bc[j][3][1] == 0) && (bc[j][4][1] == 0)
         && (v1 > u1)) {
          um = bc[j][3][3]/(bc[j][3][2] +1.0E-32);
          vm = bc[j][4][3]/(bc[j][4][2] +1.0E-32);
        }
      }
    }
    for (i = 1; i <= nel; ++i) {
      for (k = 1; k <= nnpe; ++k) {
        vi[i][k][5] = um;
        vi[i][k][6] = vm;
      }
    }

    // Phase distribution
    for (i = 1; i <= nnp; ++i) {
      fl[i][1] = 1.0;
      ph[i][1] = 3;
    }
    for (i = 1; i <= nsrf; ++i) {
      n = bgn[i][1];
      m = bgn[i][2];
      vdn = um*(x[m][2] -x[n][2]) +vm*(x[m][1] -x[n][1]);
      for (k = 1; k <= 2; ++k) {
        j = 2*(i -1) +k;
        if ((bc[j][3][1] == 0) && (bc[j][3][3] == 0) &&
         (bc[j][4][1] == 0) && (bc[j][4][3] == 0) &&
         (vdn > 0.0)) {
          for (n = 1; n <= nnpe; ++n) {
            m = ie[bel[i]][n];
            fl[m][1] = 0.0;
            ph[m][1] = 2;
          }
        }
      }
    }
  }

  // Find matrix bandwidth
  maxdif = 0;
  for (n = 1; n <= nel; ++n) {
    for (i = 1; i <= nnpe; ++i) {
      for (k = 1; k <= nnpe; ++k) {
        dif = abs(ie[n][i] -ie[n][k]);
        if (dif > maxdif) {maxdif = dif;}
      }
    }
  }
  band1 = 2*(maxdif +1) -1;

  // Neighbour node array
  for (i = 1; i <= nnp; ++i) {
    nnb[i][1] = 1;
    nnb[i][2] = i;
  }
  for (e = 1; e <= nel; ++e) {
    for (j1 = 1; j1 <= nnpe; ++j1) {
      i1 = ie[e][j1];
      for (j2 = 1; j2 <= nnpe; ++j2) {
        i2 = ie[e][j2];
        iw = 0;
        j4 = 2;

        // Global node already in list
        for (j3 = 2; j3 <= nnb[i1][1]+1; ++j3) {
          if (nnb[i1][j3] == i2) {
            iw = 1;
            j4 = j3 -1;
          }
        }

        // Node not in list and not diagonal node
        if ((iw == 0) && (nnb[i1][2] != i2)) {
          nnb[i1][1] = nnb[i1][1] +1;
          nnb[i1][nnb[i1][1]+1] = i2;
          j4 = nnb[i1][1];
        }
      }
    }
  }
}

// Re-initialize matrices
void nullmx(double r[3*nnpm], double z[3*nnpm],
 double c[3*nnpm][6*nym+11])
{
  int i, j;
  for (i = 0; i < 3*nnpm; ++i) {
    r[i] = 0.0;
    z[i] = 0.0;
  }

  for (i = 0; i < 3*nnpm; ++i) {
    for (j = 0; j < 6*nym+11; ++j) {
      c[i][j] = 0.0;
    }
  }
}

// Shape functions: bilinear quadrilateral elements
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
}

// Local stiffness matrices: temperature equation
void stifft(int e, double df, double fo, double cip[npe1],
 double tip[npe1], double flux[4][npe1][npe1], double sc[npe1],
 double aq[4][npe1][npe1], double rq[npe1])
{
  int i, j, k, n, m;
  double sct, ti, flip, kx, ky, ste, jac, rr, kkx, kky;
  double dx, dy, dqdx, dqdy, mflow, dfsq, advq, cval, hbar;
  double tphs, tphl, ec1ip, ec2ip, ec3ip, scadv, ci, tm;

  // Update customized properties
  ci = cip[1];
  kx = pks;
  ky = pks;
  ti = tmin +tip[1]*fabs(tmax-tmin);
  tm = fo*lr*lr*prl*phl/pkl;
	for (j = 1; j <= csp[1][1]; ++j) {
    k = int(csp[j+1][2]);
    custom(1, j, tm, ci, ti, csp, css, csb, cval);
    if (csp[j+1][1] == 21) {
      if (k == 2) {kx = cval;}
      if (k == 3) {ky = cval;}
      if ((k == 1) || (k > 3)) {
        kx = cval;
        ky = cval;
      }
    }
    if (csp[j+1][1] == 22) {phs = cval;}
  }

  ste = phl*fabs(tmax-tmin)/(pl + 1.0E-32);
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    m = ie[e][o[i]];
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    sc[i] = 0.0;
    shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);

    // Integration point reference states
    tphs = max(tmlt -cip[i]*(tmlt -eps -tsol)/pcs -eps, tsol);
    tphl = tmlt -cip[i]*(tmlt +eps -tsol)/(pcl +1.0E-32) +eps;
    hbar = 0.5*(phs+phl);

    // Solid phase reference states
    if (tip[i] <= tphs) {
      ec1ip = 0.0;
      ec2ip = phs/(phl +1.0E-32);
      ec3ip = 0.0;
    }

    // Liquid phase reference states
    else if (tip[i] >= tphl) {
      ec1ip = hbar*(tphl-tphs)/(phl +1.0E-32) 
       +(phs*tphs)/(phl +1.0E-32) +1.0/(ste +1.0E-32);
      ec2ip = phl/(phl +1.0E-32);
      ec3ip = tphl;
    }

    // Melt phase reference states
    else {
      ec1ip = (phs*tphs)/(phl + 1.0E-32);
      ec2ip = hbar/(phl +1.0E-32) +1.0/(ste*(tphl -tphs)
       +1.0E-32);
      ec3ip = tphs;
    }
    flip = 2.0*fl[n][1]*fl[m][1]/(fl[n][1] +fl[m][1] +1.0E-32);
    rr = ((1.0-flip)*prs +flip*prl)/(prl +1.0E-32);
    kkx = ((1.0-flip)*kx +flip*pkl)/(pkl +1.0E-32);
    kky = ((1.0-flip)*ky +flip*pkl)/(pkl +1.0E-32);

    // Local entropy based corrections
    if (ao == 1) {kkx = vps[n][11];}
    if (ao == 1) {kky = vps[n][11];}
    mflow = rr*vi[e][i][1]*ds[i][2] -rr*vi[e][i][2]*ds[i][1];

    // Diffusive, advective and source terms
    for (j = 1; j <= nnpe; ++j) {
      n = ie[e][j];
      dfsq = kkx*ds[i][2]*dn[j][3] -kky*ds[i][1]*dn[j][4];
      advq = -ic[e][1][1][i][4+j]*ec2ip*mflow;
      flux[1][i][j] = vr*advq +dfsq;
      scadv = mflow*ic[e][1][1][i][4+j]*(ec1ip -ec2ip*ec3ip);
      sc[i] = sc[i] +vr*scadv;
    }
  }

  // Forming elemental stiffness matrix
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    rr = ((1.0-fl[n][1])*prs +fl[n][1]*prl)/(prl +1.0E-32);
    for (j = 1; j <= nnpe; ++j) {
      if (i == j) {aq[1][i][j] = vr*rr*ec[n][2]*ds[i][3]/df
       -flux[1][i][j] +flux[1][o[i]][j];}
      if (i != j) {aq[1][i][j] = -flux[1][i][j]
       +flux[1][o[i]][j];}
    }
	
    // Source terms (customized)
    sct = 0.0;
    ci = cip[i];
    ti = tmin +tip[i]*fabs(tmax-tmin);
    tm = fo*lr*lr*prl*phl/pkl;
    for (j = 1; j <= css[1][1]; ++j) {
      k = int(css[j+1][1]);
      if (k == 12) {custom(2, j, tm, ci, ti, csp, css,
       csb, sct);}
    }

    // Transient and source terms (phase change + customized)
    rq[i] = -vr*rr*((ec[n][1] -ec[n][2]*ec[n][3])*ds[i][3]/df
     -(ec[n][4] +ec[n][5]*(tn[n][2]-ec[n][6]))*ds[i][3]/df)
     +sct*ds[i][3]*(lr*lr/(pkl*dtr*vr)) +(sc[o[i]] -sc[i]);
  }
}

// Diffusion coefficients: two-phase concentration equation
void coeffc(int e, double cph[npe1][3], double rd[npe1][5])
{
  int i, j;
  double dds, ddl, tphs, tphl, cphs, cphl;
  double res[npe1][3][3];

  dds = pds0/(pdl0 +1.0E-32);
  ddl = pdl0/(pdl0 +1.0E-32);
  for (j = 1; j <= nnpe; ++j) {
    i = ie[e][j];

    // Binary phase diagram: solidus and liquidus lines
    tphs = max(tmlt -cn[i][1]*(tmlt -eps -tsol)/pcs -eps, tsol);
    tphl = max(tmlt -cn[i][1]*(tmlt +eps -tsol)/pcl +eps, tsol);
    cphs = pcs*(tmlt -eps -tn[i][1])/(tmlt -eps -tsol);
    cphl = pcl*(tmlt +eps -tn[i][1])/(tmlt +eps -tsol);

    // Phase compositions and resistances
    if (fl[i][1] == 0.0) {
      cph[j][1] = cn[i][1];
      cph[j][2] = 0.0;
      res[j][1][1] = 1.0/(dds +1.0E-32);
      res[j][1][2] = 1.0/(dds +1.0E-32);
      res[j][2][1] = 1.0E+08;
      res[j][2][2] = 1.0E+08;
    }
    else if (fl[i][1] == 1.0) {
      cph[j][1] = 0.0;
      cph[j][2] = cn[i][1];
      res[j][1][1] = 1.0E+08;
      res[j][1][2] = 1.0E+08;
      res[j][2][1] = 1.0/(ddl +1.0E-32);
      res[j][2][2] = 1.0/(ddl +1.0E-32);
    }
    else {
      cph[j][1] = cphs;
      cph[j][2] = cphl;
      res[j][1][1] = 1.0/(dds*(1.0-fl[i][1]) +1.0E-32);
      res[j][1][2] = 1.0/(dds*(1.0-fl[i][1]) +1.0E-32);
      res[j][2][1] = 1.0/(ddl*fl[i][1] +1.0E-32);
      res[j][2][2] = 1.0/(ddl*fl[i][1] +1.0E-32);
    }
  }

  // Diffusion coeffients: solid / liquid phases
  for (j = 1; j <= nnpe; ++j) {
    rd[j][1] = 2.0*ds[j][2]/(fabs(res[j][1][1])
     +fabs(res[o[j+2]][1][1]));
    rd[j][2] = 2.0*ds[j][1]/(fabs(res[j][1][2])
     +fabs(res[o[j+2]][1][2]));
    rd[j][3] = 2.0*ds[j][2]/(fabs(res[j][2][1])
     +fabs(res[o[j+2]][2][1]));
    rd[j][4] = 2.0*ds[j][1]/(fabs(res[j][2][2])
     +fabs(res[o[j+2]][2][2]));
  }
}

// Integration point convection operators
void ipoint(int eqn, int e, double dfs, double xip[npe1][3],
 double ic[nelm][3][3][npe1][2*npe1])
{
  int i, j, k, m, n;
  int iip1, iip2, np1, np2, np3, np4;
  double rr, ai, bi, x1, y1, x2, y2, xm, ym, piv;
  double xc, yc, sm, mflow, vmag, xi, yi, p, flip;
  double l, u1, u2, v1, v2, lt2, ld2, dx, dy, msum;
  double pm, dqdx, dqdy, jac, a, b, m1, m2, fct[npe1];
  double sl[nelm][npe1][4], cc[npe1][2*npe1];

  // Re-initialize condensing coefficients
  for (i = 1; i <= nnpe; ++i) {
    fct[i] = 0.0;
    for (n = 1; n <= 2*nnpe; ++n) {
      cc[i][n] = 0.0;
      for (k = 1; k <= 2; ++k) {
        for (m = 1; m <= 2; ++m) {
          ic[e][m][k][i][n] = 0.0;
        }
      }
    }
  }

  // Streamwise direction
  for (j = 1; j <= nnpe; ++j) {
    n = ie[e][j];
    u1 = 0.0;
    v1 = 0.0;
    u2 = 0.0;
    v2 = 0.0;
    s[1][1] = lc[j][1];
    s[1][2] = lc[j][2];
    rr = ((1.0 -fl[n][1])*prs +fl[n][1]*prl)/(prl +1.0E-16);
    shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);

    // Droplet flow velocities
    if (eqn == 4) {
      m = o[j+2];
      u1 = vi[e][j][3];
      v1 = vi[e][j][4];
      u2 = vi[e][m][3];
      v2 = vi[e][m][4];
    }

    // Main flow velocities
    else {
      for (i = 1; i <= nnpe; ++i) {
        n = ie[e][i];
        m = o[j+2];
        u1 = u1 +nn[i][j]*vps[n][6];
        v1 = v1 +nn[i][j]*vps[n][7];
        u2 = u2 +nn[i][m]*vps[n][6];
        v2 = v2 +nn[i][m]*vps[n][7];
      }
    }
    vmag = sqrt(u1*u1 +v1*v1);
    ai = u1/(vmag +1.0E-16);
    bi = v1/(vmag +1.0E-16);

    // Sub-control-volume coordinates
    mflow = u1*ds[j][2] -v1*ds[j][1];
    if (mflow < 0.0) {
      np1 = o[j+1];
      np2 = o[j+2];
      np3 = o[j+3];
      np4 = o[j];
      iip1 = o[j+2];
      iip2 = o[j];
    }
    else if (mflow >= 0.0) {
      np1 = o[j+2];
      np2 = o[j+1];
      np3 = o[j];
      np4 = o[j+3];
      iip1 = o[j];
      iip2 = o[j+2];
    }
    x1 = x[ie[e][np1]][1];
    y1 = x[ie[e][np1]][2];
    x2 = x[ie[e][np2]][1];
    y2 = x[ie[e][np2]][2];
    xm = 0.5*(x2 +x[ie[e][np3]][1]);
    ym = 0.5*(y2 +x[ie[e][np3]][2]);
    xc = 0.25*(x1 +x2 +x[ie[e][o[j+3]]][1] +x[ie[e][o[j]]][1]);
    yc = 0.25*(y1 +y2 +x[ie[e][o[j+3]]][2] +x[ie[e][o[j]]][2]);

    // Equate lines at streamline / element intersection
    p = (bi*(xip[j][1] -xm) -ai*(xip[j][2] -ym))/((x2 -xm)*bi
     -(y2 -ym)*ai +1.0E-16);

    // Adjacent external surface upwinding
    if (p > 1.0) {
      p = (bi*(xip[j][1] -x2) -ai*(xip[j][2] -y2))/(((x1 -x2)
       *bi -(y1 -y2)*ai) +1.0E-16);
      xi = x2 +p*(x1 -x2);
      yi = y2 +p*(y1 -y2);
      b = sqrt((x1 -x2)*(x1 -x2) +(y1 -y2)*(y1 -y2));
      a = b*(1.0 -fabs(p));
      l = sqrt((xip[j][1] -xi)*(xip[j][1] -xi) +(xip[j][2]
       -yi)*(xip[j][2] -yi));
      ic[e][1][1][j][np2] = (a/b)*rr*vmag/(l +1.0E-16);
      ic[e][1][1][j][np1] = (1.0 -a/b)*rr*vmag/(l +1.0E-16);
      ic[e][2][1][j][np2] = (a/b)*rr*vmag/(l +1.0E-16);
      ic[e][2][1][j][np1] = (1.0 -a/b)*rr*vmag/(l +1.0E-16);
    }

    // Opposite external surface upwinding
    else if ((p >= 0.0) && (p <= 1.0)) {
      xi = xm +p*(x2 -xm);
      yi = ym +p*(y2 -ym);
      b = 2.0*sqrt((x2 -xm)*(x2 -xm) +(y2 -ym)*(y2 -ym));
      a = 0.5*b*(1.0 +fabs(p));
      l = sqrt(pow(xip[j][1] -xi, 2.0) +pow(xip[j][2] -yi, 2.0));
      ic[e][1][1][j][np2] = (a/b)*rr*vmag/(l +1.0E-16);
      ic[e][1][1][j][np3] = (1.0 -a/b)*rr*vmag/(l +1.0E-16);
      ic[e][2][1][j][np2] = (a/b)*rr*vmag/(l +1.0E-16);
      ic[e][2][1][j][np3] = (1.0 -a/b)*rr*vmag/(l +1.0E-16);
    }

    // Adjacent interior surface upwinding
    else {
      p = (bi*(xip[j][1] -xm) -ai*(xip[j][2] -ym))/(((xc -xm)*bi
       -(yc -ym)*ai) +1.0E-16);

      // Outer interior section upwinding
      if (p <= 0.5) {
        xi = xm +p*(xc -xm);
        yi = ym +p*(yc -ym);
        b = 0.5*sqrt((xc -xm)*(xc -xm) +(yc -ym)*(yc -ym));
        a = b*(1.0 -2.0*fabs(p));
        l = sqrt(pow(xip[j][1] -xi, 2.0) +pow(xip[j][2]
         -yi, 2.0));
        ic[e][1][1][j][np3] = 0.5*(a/b)*rr*vmag/(l +1.0E-16);
        ic[e][1][1][j][np2] = 0.5*(a/b)*rr*vmag/(l +1.0E-16);
        ic[e][2][1][j][np3] = 0.5*(a/b)*rr*vmag/(l +1.0E-16);
        ic[e][2][1][j][np2] = 0.5*(a/b)*rr*vmag/(l +1.0E-16);
        cc[j][iip1] = -(1.0 -a/b)*rr*vmag/(l +1.0E-16);
      }

      // Inner interior section upwinding
      else {
        xi = xm +p*(xc-xm);
        yi = ym +p*(yc-ym);
        b = sqrt(pow(xip[iip1][1] -xip[iip2][1], 2.0)
         +pow(xip[iip1][2] -xip[iip2][2], 2.0) +1.0E-16);
        a = sqrt(pow(xip[iip2][1] -xi, 2.0) +pow(xip[iip2][2]
         -yi, 2.0));
        l = sqrt(pow(xip[j][1] -xi, 2.0) +pow(xip[j][2]
         -yi, 2.0));
        cc[j][iip1] = -(a/b)*rr*vmag/(l +1.0E-16);
        cc[j][iip2] = -(1.0 -a/b)*rr*vmag/(l +1.0E-16);
      }
    }

    // Multiphase region permeability
    m = ie[e][o[j]];
    flip = 2.0*fl[n][1]*fl[m][1]/(fl[n][1] +fl[m][1] +1.0E-16);
    pm = 1.0/(1.0 +k0*(1.0 -flip)*(1.0 -flip)/(flip*flip*flip
     +1.E-16));
    if ((eqn <= 2) || (eqn == 4)) {pm = 1.0;}

    // Mass weighted (positive coefficient) upwinding
    lt2 = ds[j][1]*ds[j][1] +ds[j][2]*ds[j][2];
    ld2 = 0.5*ds[j][3]*ds[j][3]/lt2 +0.375*lt2;
    m1 = u1*ds[j][2] -v1*ds[j][1];
    m2 = u2*ds[o[j+2]][2] -v2*ds[o[j+2]][1];
    if (m1 >= 0.0) {sm = max(min(-m2/(m1 +1.0E-16), 1.0), 0.0);}
    if (m1 < 0.0) {sm = max(min(m2/(m1 +1.0E-16), 1.0), 0.0);}
    if (eqn <= 1) {
      for (n = 1; n <= 2*nnpe; ++n) {
        cc[j][n] = 0.0;
        for (k = 1; k <= 2; ++k) {
          for (m = 1; m <= 2; ++m) {
            ic[e][m][k][j][n] = 0.0;
          }
        }
      }
      cc[j][j] = pm*rr*vmag/(l +1.0E-16) +dfs/(ld2 +1.0E-16);
      cc[j][iip1] = -sm*rr*vmag/l;
      cc[j][4+j] = 1.0;
      ic[e][1][1][j][np2] = (1.0 -sm)*rr*vmag/(l +1.0E-16);
      ic[e][2][1][j][np2] = (1.0 -sm)*rr*vmag/(l +1.0E-16);
    }

    // Skew upwind convection operator
    if (eqn > 1) {
      cc[j][j] = pm*rr*vmag/(l +1.0E-16) +dfs/(ld2 +1.0E-16);
      cc[j][4+j] = 1.0;
      sl[e][j][1] = sm;
      sl[e][j][2] = l;
      sl[e][j][3] = vmag;

      // Nodal velocity influence coefficients
      for (i = 1; i <= nnpe; ++i) {
        ic[e][1][1][j][i] = pm*ic[e][1][1][j][i] +dfs*nn[i][j]
         /(ld2 +1.0E-16);
        ic[e][2][1][j][i] = pm*ic[e][2][1][j][i] +dfs*nn[i][j]
         /(ld2 +1.0E-16);
      }
    }

    // Nodal pressure influence coefficients
    for (i = 1; i <= nnpe; ++i) {
      ic[e][1][2][j][i] = -pm*dn[i][3];
      ic[e][2][2][j][i] = -pm*dn[i][4];
    }
  }

  // Integration point -> nodal dependence (matrix inversion)
  for (i = 1; i <= nnpe; ++i) {
    piv = cc[i][i];
    for (m = 1; m <= 2*nnpe; ++m) {
      cc[i][m] = cc[i][m]/piv;
      for (k = 1; k <= nnpe; ++k) {
        if ((k != i) && (m == i)) {fct[k] = cc[k][i]/cc[i][i];}
        if (k != i) {cc[k][m] = cc[k][m] -fct[k]*cc[i][m];}
      }
    }
  }

  // Isolate integration point variable dependence
  for (m = 1; m <= 2; ++m) {
    for (n = 1; n <= 2; ++n) {
      for (i = 1; i <= nnpe; ++i) {
        for (k = 1; k <= nnpe; ++k) {
          msum = 0.0;
          for (j = 1; j <= nnpe; ++j) {
            msum = msum +cc[i][4+j]*ic[e][m][n][j][k];
          }
          ic[e][m][n][i][4+k] = msum;
        }
      }
    }
  }
}

// Iteration rules for solid / liquid phase change
void phase(int nnp, int ph[nnpm][3], double fl[nnpm][3],
 double ec[nnpm][7], int &tr)
{
  int i, j, r2;
  double ste, tphs, tphl, cphs, cphl, hbar;

  hbar = 0.5*(phl +phs);
  ste = phl*dtr/(pl +1.0E-32);
  for (i = 1; i <= nnp; ++i) {
    tphs = max(tmlt -eps -cn[i][1]*(tmlt -eps -tsol)/pcs, tsol);
    tphl = max(tmlt +eps -cn[i][1]*(tmlt +eps -tsol)/pcl, tsol);
    cphs = pcs*(tmlt -eps -tn[i][1])/(tmlt -eps -tsol);
    cphl = pcl*(tmlt +eps -tn[i][1])/(tmlt +eps -tsol);
    if (tn[i][1] <= tphs) {
      ph[i][1] = 1;
      fl[i][1] = 0.0;
    }
    else if (tn[i][1] >= tphl) {
      ph[i][1] = 3;
      fl[i][1] = 1.0;
    }
    else {
      ph[i][1] = 2;
      fl[i][1] = (tn[i][1] -tphs)/(tphl -tphs);
    }
  }

  // Compare tentative vs computed phases
  tr = 0;
  for (i = 1; i <= nnp; ++i) {
    if (ph[i][1] != ph[i][2]) {tr = 1;}
  }

  // Rule 1: Pass through melt phase
  if (tr == 1) {
    for (i = 1; i <= nnp; ++i) {
      if (fabs(ph[i][1] -ph[i][2]) > 1) {ph[i][1] = 2;}
    }
  }

  // Rule 2: Adjacent CV phase change criterion
  if (tr == 1) {
    for (i = 1; i <= nnp; ++i) {
      r2 = 0;
      if (nnb[i][1] == 9) {
        for (j = 2; j <= nnb[i][1]+1; ++j) {
          if (ph[nnb[i][j]][2] != ph[i][2]) {r2 = 1;}
        }
        if ((r2 == 0) && (ph[i][1] != ph[i][2]))
          {ph[i][1] = ph[i][2];}
      }
    }
  }

  // Energy equation of state
  for (i = 1; i <= nnp; ++i) {
    tphs = tmlt -cn[i][1]*(tmlt-eps -tsol)/(pcs +1.0E-32) -eps;
    tphl = tmlt -cn[i][1]*(tmlt+eps -tsol)/(pcl +1.0E-32) +eps;
    cphs = pcs*(tmlt -eps -tn[i][1])/(tmlt -eps -tsol);
    cphl = pcl*(tmlt +eps -tn[i][1])/(tmlt +eps -tsol);

    // Solid phase reference states
    if (ph[i][1] == 1) {
      ec[i][1] = 0.0;
      ec[i][2] = phs/(phl +1.0E-32);
      ec[i][3] = 0.0;
    }

    // Liquid phase reference states
    else if (ph[i][1] == 3) {
      ec[i][1] = hbar*(tphl -tphs)/(phl +1.0E-32)
       +(phs*tphs)/(phl +1.0E-32) +1.0/(ste +1.0E-32);
      ec[i][2] = phl/phl;
      ec[i][3] = tphl;
    }

    // Melt phase reference states
    else {
      ec[i][1] = (phs*tphs)/(phl +1.0E-32);
      ec[i][2] = hbar/(phl +1.0E-32) +1.0/(ste*(tphl -tphs)
       +1.0E-32);
      ec[i][3] = tphs;
    }
  }
}

// Local stiffness matrices: C equation (solid / liquid)
void stiffc(int e, double df, double fo, double xip[npe1][3],
 double flux[4][npe1][npe1], double sc[npe1],
 double aq[4][npe1][npe1], double rq[npe1])
{
  int i, j, k, n;
  double dfsc, advc, ci, ti, tm, le, jac, dfs, scc;
  double dx, dy, dqdx, dqdy, mflow, scadv, scdfs;
  double clm[npe1][4], rd[npe1][5], cph[npe1][3];

  le = pkl/(prl*phl*pdl0 +1.0E-32);
  dfs = 1.0/le;
  coeffc(e, cph, rd);
  ipoint(2, e, dfs, xip, ic);
  for (i = 1; i <= nnpe; ++i) {
    clm[i][1] = 0.0;
    clm[i][2] = 0.0;
    clm[i][3] = 0.0;
    for (n = 1; n <= nnpe; ++n) {
      clm[i][1] = clm[i][1] +nn[n][i]*cn[ie[e][n]][1];
      clm[i][2] = clm[i][2] +nn[n][i]*cph[i][2];
      clm[i][3] = clm[i][3] +nn[n][i]*fl[ie[e][n]][1];
    }
  }

  // Outer loop: local node dependence
  for (i = 1; i <= nnpe; ++i) {
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    sc[i] = 0.0;
    mflow = vi[e][i][1]*ds[i][2] -vi[e][i][2]*ds[i][1];
    shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);

    // Diffusion, advection, source terms
    for (j = 1; j <= nnpe; ++j) {
      n = ie[e][j];
      // dfsc = rd[i][1]*dn[j][3] -rd[i][2]*dn[j][4]
      // +rd[i][3]*dn[j][3] -rd[i][4]*dn[j][4];
      dfsc = (ds[i][2]*dn[j][3] -ds[i][1]*dn[j][4])
       *(fl[i][1]*pdl0 +(1.0 -fl[i][1])*pds0)/(pdl0 +1.0E-32);
      advc = -le*ic[e][1][1][i][4+j]*mflow;
      scdfs = (rd[i][1]*dn[j][3] -rd[i][2]*dn[j][4])
       *(cph[j][1] -cn[n][1]) -(rd[i][3]*dn[j][3]
       -rd[i][4]*dn[j][4])*(cph[j][2] -cn[n][1]);
      scadv = le*mflow*(clm[j][3]*clm[j][2] -clm[j][1]);
      sc[i] = sc[i] +vr*scadv +scdfs;
      flux[1][i][j] = vr*advc +dfsc;
    }
  }

  // Forming elemental stiffness matrix
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    for (j = 1; j <= nnpe; ++j) {
      if (i == j) {aq[1][i][j] = vr*le*ds[i][3]/df
       -flux[1][i][j] +flux[1][o[i]][j];}
      if (i != j) {aq[1][i][j] = -flux[1][i][j]
       +flux[1][o[i]][j];}
    }

    // Source terms (customized)
    scc = 0.0;
    ci = cn[n][1];
    ti = tmin +tn[n][1]*fabs(tmax-tmin);
    tm = fo*lr*lr*prl*phl/pkl;
    for (j = 1; j <= css[1][1]; ++j) {
      k = int(css[j+1][1]);
      if (k == 11) {custom(2, j, tm, ci, ti, csp, css,
       csb, scc);}
    }

    // Transient and source terms (phase change + customized)
    rq[i] = vr*le*cn[n][2]*ds[i][3]/df +(sc[o[i]] -sc[i])
     +scc*ds[i][3]*(lr*lr*prl*phl/(pkl*vr));
  }
}

// Local stiffness matrices: u-v-P equation (continuity)
void stiffp(int e, double flux[4][npe1][npe1], double sc[npe1],
 double aq[4][npe1][npe1], double rq[npe1])
{
  int i, j, n, m;
  double flip, rr, advu, advv, advp;
  double dqdx, dqdy, jac, dx, dy;

  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i]; 
    m = ie[e][o[i]];
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    sc[i] = 0.0;
    flip = 2.0*fl[n][1]*fl[m][1]/(fl[n][1] +fl[m][1] +1.0E-32);
    rr = ((1.0-flip)*prs +flip*prl)/(prl +1.0E-32);
    shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);

    // Mass flux across 2-D sub-surface
    for (j = 1; j <= nnpe; ++j) {
      advu = ic[e][1][1][i][4+j]*ds[i][2];
      advv = -ic[e][2][1][i][4+j]*ds[i][1];
      advp = ic[e][1][2][i][4+j]*ds[i][2] -ic[e][2][2][i][4+j]
       *ds[i][1];
      flux[1][i][j] = advu;
      flux[2][i][j] = advv;
      flux[3][i][j] = advp;
    }
  }

  // Forming elemental stiffness matrix
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    for (j = 1; j <= nnpe; ++j) {
      aq[1][i][j] = -flux[1][i][j] +flux[1][o[i]][j];
      aq[2][i][j] = -flux[2][i][j] +flux[2][o[i]][j];
      aq[3][i][j] = -flux[3][i][j] +flux[3][o[i]][j];
    }
    rq[i] = 0.0;
  }
}

// Local stiffness matrices: U-v-p (x-momentum)
void stiffu(int e, double df, double xip[npe1][3],
 double flux[4][npe1][npe1], double sc[npe1],
 double aq[4][npe1][npe1], double rq[npe1])
{
  int i, j, n, m;
  double mflow, dfs, advu, advv, advp, dx, dy, jac;
  double dqdx, dqdy, flip, dfsu, dfsv, pp, ppr, rr;

  // Integration point influence coefficients
  dfs = 0.25*(pr[e][1] +pr[e][2] +pr[e][3] +pr[e][4]);
  ipoint(3, e, dfs, xip, ic);

  // Start outer loop: local node dependence
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    m = ie[e][o[i]];
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    sc[i] = 0.0;
    flip = 2.0*fl[n][1]*fl[m][1]/(fl[n][1] +fl[m][1] +1.0E-16);
    rr = ((1.0 -flip)*prs +flip*prl)/(prl +1.0E-16);
    mflow = rr*vi[e][i][1]*ds[i][2] -rr*vi[e][i][2]*ds[i][1];
    shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);

    // Prandtl number: molecular, turbulent + two-phase parts
    pp = k0*(1.0 -flip)*(1.0 -flip)/(flip*flip*flip +1.0E-16);
    ppr = pr[e][i] +pp;
    for (j = 1; j <= nnpe; ++j) {

      // Diffusion, advection flux
      advu = -ic[e][1][1][i][4+j]*mflow;
      dfsu = 2.0*ppr*dn[j][3]*ds[i][2] -ppr*dn[j][4]*ds[i][1];
      advv = 0.0;
      dfsv = -ppr*dn[j][3]*ds[i][1];
      advp = -ic[e][1][2][i][4+j]*mflow -nn[j][i]*ds[i][2];
      flux[1][i][j] = vr*advu +dfsu;
      flux[2][i][j] = vr*advv +dfsv;
      flux[3][i][j] = vr*advp;
    }
  }


  // Forming elemental stiffness matrix
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    rr = ((1.0-fl[n][1])*prs +fl[n][1]*prl)/(prl +1.0E-16);
    for (j = 1; j <= nnpe; ++j) {
      if (i == j) {aq[1][i][j] = vr*rr*ds[i][3]/df
       -flux[1][i][j] +flux[1][o[i]][j];}
      if (i != j) {aq[1][i][j] = -flux[1][i][j]
       +flux[1][o[i]][j];}
      aq[2][i][j] = -flux[2][i][j] +flux[2][o[i]][j];
      aq[3][i][j] = -flux[3][i][j] +flux[3][o[i]][j];
    }

    // Transient and source terms (phase interactions)
    rq[i] = vr*rr*vps[n][4]*ds[i][3]/df +(sc[o[i]] -sc[i]);
  }
}

// Local stiffness matrices: u-V-p (y-momentum)
void stiffv(int e, double df, double flux[4][npe1][npe1],
 double sc[npe1], double aq[4][npe1][npe1], double rq[npe1])
{
  int i, j, n, m;
  double rr, tbar, cbar, pm, flip, dqdx, dqdy, dx, dy, pp;
  double jac, scf, advu, advv, advp, dfsu, dfsv, ppr, mflow;

  pm = (vsc*prl*phl)/(pkl +1.0E-16);
  cbar = 0.5*(pcs+ pcl);
  tbar = 0.5;

  // Start outer loop: local node dependence
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    m = ie[e][o[i]];
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    sc[i] = 0.0;
    flip = 2.0*fl[n][1]*fl[m][1]/(fl[n][1] +fl[m][1] +1.0E-16);
    rr = ((1.0 -flip)*prs +flip*prl)/(prl +1.0E-16);
    mflow = rr*vi[e][i][1]*ds[i][2] -rr*vi[e][i][2]*ds[i][1];
    shape(e, nnpm, dx, dy, s, dn, ivo, jac, dqdx, dqdy);

    // Prandtl number: molecular, turbulent + two-phase parts
    pp = k0*(1.0 -flip)*(1.0 -flip)/(flip*flip*flip+1.0E-16);
    ppr = pr[e][i] +pp;
    for (j = 1; j <= nnpe; ++j) {

      // Diffusion, advection flows
      advu = 0.0;
      dfsu = ppr*dn[j][4]*ds[i][2];
      advv = -ic[e][2][1][i][4+j]*mflow;
      dfsv = ppr*dn[j][3]*ds[i][2] -2.0*ppr*dn[j][4]*ds[i][1];
      advp = -ic[e][2][2][i][4+j]*mflow +nn[j][i]*ds[i][1];
      flux[1][i][j] = vr*advu +dfsu;
      flux[2][i][j] = vr*advv +dfsv;
      flux[3][i][j] = vr*advp;
    }
  }

  //  Forming elemental stiffness matrix
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    rr = ((1.0-fl[n][1])*prs +fl[n][1]*prl)/(prl +1.0E-16);
    for (j = 1; j <= nnpe; ++j) {
      if (i == j) {aq[2][i][j] = vr*rr*ds[i][3]/df
       -flux[2][i][j] +flux[2][o[i]][j];}
      if (i != j) {aq[2][i][j] = -flux[2][i][j]
       +flux[2][o[i]][j];}
      aq[1][i][j] = -flux[1][i][j] +flux[1][o[i]][j];
      aq[3][i][j] = -flux[3][i][j] +flux[3][o[i]][j];
    }

    // Sources (body forces and phase change), transient term
    scf = rat*pm*ds[i][3]*(tn[n][1] -tbar) +ras*pm*ds[i][3]
     *(cn[n][1] -cbar);
    rq[i] = vr*rr*vps[n][5]*ds[i][3]/df +(sc[o[i]]
     -sc[i]) +scf/vr;
  }
}

// Updated integration point velocities
void uvcor(int nel, double vi[nelm][npe1][7])
{
  int e, n, i, j;
  double uic, vic, pxc, pyc;

  // Re-construct momentum equations
  for (e = 1; e <= nel; ++e) {
    for (j = 1; j <= nnpe; ++j) {
      n = ie[e][j];
      uic = 0.0;
      vic = 0.0;
      pxc = 0.0;
      pyc = 0.0;
      for (i = 1; i <= nnpe; ++i) {
        n = ie[e][i];
        uic = uic +ic[e][1][1][j][4+i]*vps[n][1];
        vic = vic +ic[e][2][1][j][4+i]*vps[n][2];
        pxc = pxc +ic[e][1][2][j][4+i]*vps[n][3];
        pyc = pyc +ic[e][2][2][j][4+i]*vps[n][3];
      }
      vi[e][j][1] = uic +pxc;
      vi[e][j][2] = vic +pyc;
    }
  }
}

// Control Volume entropy production rate
void entropy(double df, int nel, int nnp, int nsrf,
 double vps[nnpm][14])
{
  int e, i, j, k, n, jb, m, n1, n2, n3;
  double rhon, rhoip, kip, dx, dy, jac, tdim, ss, tt;
  double dtdx, dtdy, dudx, dudy, dvdx, dvdy, sf, ur, dtm;
  double tphs, tphl, teut, sr1, sr2, sr3, cpn;
  double tpip[npe1], sip[npe1], flux[4][npe1][npe1];

  ur = vr*pkl/(prl*phl*lr +1.0E-32);
  dtm = df*lr*lr*prl*phl/pkl;
  for (n = 1; n <= nnp; ++n) {
    for (i = 8; i <= 11; ++i) {
      vps[n][i] = 0.0;
    }
  }

  // Entropy equation of state
  teut = tmin +dtr*tsol;
  sf = pl/(tmin +dtr*tmlt +1.0E-32);
  for (n = 1; n <= nnp; ++n) {
    tdim = tmin +dtr*tn[n][1];
    tphs = tmin +dtr*(tmlt -eps -cn[n][1]*(tmlt -eps -tsol)
     /(pcs -eps));
    tphl = tmin +dtr*(tmlt +eps -cn[n][1]*(tmlt +eps -tsol)
     /(pcl +eps));

    // Solid phase reference states
    if (tn[n][1] <= tphs) {
      sr1 = 0.0;
      sr2 = phs;
      sr3 = teut;
    }

    // Liquid phase reference states
    else if (tn[n][1] >= tphl) {
      sr1 = log(tphl/tphs)*(phs*tphl -phl*tphs)/(tphl -tphs)
       +phs*log(tphs/teut) +(phl -phs +sf);
      sr2 = phl;
      sr3 = tphl;
    }

    // Melt phase reference states
    else {
      sr1 = phs*log(tphs/teut);
      sr2 = (phs*tphl -phl*tphs)/(tphl -tphs) +(phl -phs +sf)
       /log(tphl/tphs);
      sr3 = tphs;
    }
    vps[n][8] = sr1 +sr2*log(tdim/(sr3 +1.0E-16));
  }

  // Assembly of sub-control-volumes
  for (e = 1; e <= nel; ++e) {
    for (j = 1; j <= nnpe; ++j) {
      n = ie[e][j];
      ivo[j] = tn[n][1];
    }

    // Local element geometry
    for (i = 1; i <= nnpe; ++i) {
      m = ie[e][i];
      s[1][1] = lc[12+i][1];
      s[1][2] = lc[12+i][2];
      s[2][1] = lc[4+i][1];
      s[2][2] = lc[4+i][2];
      s[3][1] = lc[8+i][1];
      s[3][2] = lc[8+i][2];
      shape(e, nnpm, dx, dy, s, dn, ivo, jac, dtdx, dtdy);
      ds[i][3] = fabs(jac);
      ds[i][4] = fabs(jac);
      s[1][1] = lc[i][1];
      s[1][2] = lc[i][2];
      vps[m][9] = vps[m][9] +ds[i][3]*lr*lr;
      shape(e, nnpm, dx, dy, s, dn, ivo, jac, dtdx, dtdy);
      ds[i][1] = dx;
      ds[i][2] = dy;

      // Integration point interpolation
      ss = lc[i][1];
      tt = lc[i][2];
      nn[1][i] = 0.25*(1.0 +ss)*(1.0 +tt);
      nn[2][i] = 0.25*(1.0 -ss)*(1.0 +tt);
      nn[3][i] = 0.25*(1.0 -ss)*(1.0 -tt);
      nn[4][i] = 0.25*(1.0 +ss)*(1.0 -tt);
      rhoip = 0.0;
      kip = 0.0;
      tpip[i] = 0.0;
      sip[i] = 0.0;
      for (j = 1; j <= nnpe; ++j) {
        n = ie[e][j];
        rhoip = rhoip +nn[j][i]*((1.0 -fl[n][1])*prs
         +fl[n][1]*prl);
        kip = kip +nn[j][i]*((1.0 -fl[n][1])*pks +fl[n][1]*pkl);
        tpip[i] = tpip[i] +nn[j][i]*(tmin +dtr*tn[n][1]);
        sip[i] = sip[i] +nn[j][i]*vps[n][8];
      }

      // Entropy flux calculation
      flux[1][i][1] = kip*dtr*(dtdx*ds[i][2] -dtdy*ds[i][1])
       /tpip[i] +vr*lr*rhoip*sip[i]*(vi[e][i][1]*ds[i][2]
       -vi[e][i][2]*ds[i][1]);
    }

    // SCV entropy production rate
    for (i = 1; i <= nnpe; ++i) {
      n = ie[e][i];
      rhon = (1.0 -fl[n][1])*prs +fl[n][1]*prl;
      cpn = (1.0 -fl[n][1])*phs +fl[n][1]*phl;
      tdim = tmin +dtr*tn[n][1];
      vps[n][10] = vps[n][10] +rhon*cpn*dtr*(tn[n][1]
       -tn[n][2])*ds[i][3]*lr*lr/(dtm*tdim) +rhon*sf
       *(fl[n][1] -fl[n][2])*ds[i][3]*lr*lr/dtm
       -flux[1][i][1] +flux[1][o[i]][1];
      vps[n][11] = tdim*tdim/(dtr*dtdx*dtdx/lr
       +dtr*dtdy*dtdy/lr +1.0E-16);
    }
  }

  // Entropy based diffusivity
  for (n = 1; n <= nnp; ++n) {
    vps[n][10] = vps[n][10]/(vps[n][9] +1.0E-32);
    if (vps[n][10] < 0.0) {
      vps[n][11] = fabs(vps[n][10])*vps[n][11];
			if (vps[n][11] > 1.0) {vps[n][11] = 1.0;}
    }
    else {
      vps[n][11] = 0.0;
    }
  }

  // Boundary entropy production rate
  for (n = 1; n <= nsrf; ++n) {
    e = bel[n];
    n1 = bgn[n][1];
    n2 = bgn[n][2];
    for (i = 1; i <= 2; ++i) {
      for (k = 1; k <= 3; ++k) {
        for (j = 1; j <= nnpe; ++j) {
          n3 = ie[e][j];
          if (n1 == n3) jb = j;
          if (k == 1) ivo[j] = tn[n3][1];
          if (k == 2) ivo[j] = vps[n3][1];
          if (k == 3) ivo[j] = vps[n3][2];
        }
        m = 2*(jb -1) +i;
        s[1][1] = nc[m][1];
        s[1][2] = nc[m][2];
        s[2][1] = nc[8+m][1];
        s[2][2] = nc[8+m][2];
        s[3][1] = nc[16+m][1];
        s[3][2] = nc[16+m][2];
        if (k == 1) {shape(e, nnpm, dx, dy, s, dn, ivo,
         jac, dtdx, dtdy);}
        if (k == 2) {shape(e, nnpm, dx, dy, s, dn, ivo,
         jac, dudx, dudy);}
        if (k == 3) {shape(e, nnpm, dx, dy, s, dn, ivo,
         jac, dvdx, dvdy);}
      }
      j = bgn[n][i];
      tdim = tmin +dtr*tn[j][1];
      vps[j][10] = pks*dtr*dtr*(dtdx*dtdx +dtdy*dtdy)
       /(tdim*tdim*lr*lr) +(prl*vsc)*vr*vr*(2.0*dudx*dudx
       +2.0*dvdy*dvdy +(dudy +dvdx)*(dudy +dvdx)
       -0.667*(dudx +dvdy)*(dudx +dvdy))/(tdim*lr*lr);
      vps[j][11] = 0.0;
    }
  }
}

// Smoothing of entropy based diffusivity
void dsmooth(int nel, int nnp, double vps[nnpm][14])
{
  int e, i, j, k, kn, n, m;
  double jac, dx, dy, dmdx, dmdy;
  double flux[4][npe1][npe1];

  for (n = 1; n <= nnp; ++n) {
    vps[n][12] = 0.0;
    vps[n][13] = 0.0;
  }
  kn = 5;

  // Assembly of diffusion terms
  for (k = 1; k <= kn; ++k) {
    for (e = 1; e <= nel; ++e) {
      for (i = 1; i <= nnpe; ++i) {
        m = ie[e][i];
        s[1][1] = lc[i][1];
        s[1][2] = lc[i][2];
        s[2][1] = lc[4+i][1];
        s[2][2] = lc[4+i][2];
        s[3][1] = lc[8+i][1];
        s[3][2] = lc[8+i][2];
        shape(e, nnpm, dx, dy, s, dn, ivo, jac, dmdx, dmdy);
        ds[i][1] = dx;
        ds[i][2] = dy;
        for (j = 1; j <= nnpe; ++j) {
          flux[1][i][j] = dn[j][3]*ds[i][2] -dn[j][4]*ds[i][1];
        }
      }

      // Iterative Jacobi point solver
      for (i = 1; i <= nnpe; ++i) {
        n = ie[e][i];
        for (j = 1; j <= nnpe; ++j) {
          m = ie[e][j];
          if (i == j) {
            vps[n][12] = vps[n][12] -flux[1][i][j]
             +flux[1][o[i]][j];
          }
          else {
            vps[n][13] = vps[n][13] -flux[1][i][j]*vps[m][11]
             +flux[1][o[i]][j]*vps[m][11];
          }
        }
      }
    }

    // Updated diffusivity field
    for (n = 1; n <= nnp; ++n) {
      if (vps[n][10] > 0.0) {vps[n][11] = fabs(vps[n][13]
       /(vps[n][12] +1.0E-16));}
      if (k == kn) {vps[n][11] = vps[n][11]/(pkl +1.0E-32)
       +((1.0 -fl[n][1])*pks +fl[n][1]*pkl)/(pkl +1.0E-32);}
      vps[n][12] = 0.0;
      vps[n][13] = 0.0;
    }
  }
}

// Customized properties, source terms and boundary conditions
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
}

// Drag and lift coefficients by pressure and viscous forces
void cdcl(int nnp, int nsrf, double &um, double &cd, double &cl)
{
  int e, i, j, k, m, n, p, r, n1, n2;
  double sa, dudx, dudy, dvdx, dvdy, dsx, dsy;
  double jac, dx, dy, u1, drag, lift, dels, pw;

  // Reference or maximum velocity
  um = 0.0;
  sa = 0.0;
  ur = vr*pkl/(prl*phl*lr +1.0E-32);
  for (i = 1; i <= nsrf; ++i) {
    for (k = 1; k <= 2; ++k) {
      j = 2*(i -1) +k;
      m = bgn[i][1];
      u1 = sqrt(pow(bc[j][3][3]/(bc[j][3][2] +1.0E-32), 2.0)
       +pow(bc[j][4][3]/(bc[j][4][2] +1.0E-32), 2.0));
      if ((bc[j][3][1] == 0) && (bc[j][4][1] == 0)
       && (u1 > um)) {um = u1;}
    }
  }
  if (um == 0.0) {
    for (i = 1; i <= nnp; ++i) {
      u1 = sqrt(vps[i][1]*vps[i][1] +vps[i][2]*vps[i][2]);
      if (u1 > um) {um = u1;}
    }
  }
  um = um*ur;

  // Boundary surface geometry
  drag = 0.0;
  lift = 0.0;
  for (i = 1; i <= nsrf; ++i) {
    n1 = bgn[i][1];
    n2 = bgn[i][2];
    e = bel[i];
    dels = sqrt(pow(x[n2][1] -x[n1][1], 2.0) +pow(x[n2][2]
     -x[n1][2], 2.0));
    dsx = (x[n2][1] -x[n1][1])/(dels +1.0E-32);
    dsy = (x[n2][2] -x[n1][2])/(dels +1.0E-32);
    for (k = 1; k <= 2; ++k) {
      j = 2*(i -1) +k;
      p = bgn[i][k];
      pw = vps[p][3];

      // Drag and lift coefficients on no-slip boundaries
      if ((bc[j][3][1] == 0) && (bc[j][3][3] == 0)
       && (bc[j][4][1] == 0) && (bc[j][4][3] == 0)) {
        for (r = 1; r <= 2; ++r) {
          for (n = 1; n <= nnpe; ++n) {
            m = ie[e][n];
            if (m == p) {
              s[1][1] = nc[8+2*(n-1)+1][1];
              s[1][2] = nc[8+2*(n-1)+1][2];
            }
            ivo[n] = vps[m][r];
          }
          if (r == 1) {shape(e, nnpm, dx, dy, s, dn, ivo,
           jac, dudx, dudy);}
          if (r == 2) {shape(e, nnpm, dx, dy, s, dn, ivo,
           jac, dvdx, dvdy);}
        }
        drag = drag +(0.5*lr*dels)*(dsy*prl*pow(um, 2.0)*pw
         +dsx*prl*vsc*um*(dudy +dvdx)/(lr +1.0E-32));
        lift = lift +(0.5*lr*dels)*(-dsx*prl*pow(um, 2.0)*pw
         +dsy*prl*vsc*um*(dudy +dvdx)/(lr +1.0E-32));
        sa = sa +0.5*dels;
      }
    }
  }
  cd = fabs(drag/(0.5*prl*um*um*lr*sa +1.0E-32));
  cl = lift/(0.5*prl*um*um*lr*sa +1.0E-32);
}

// Element momentum equations: droplet trajectories
void stiffdv(int e, double vi[nelm][npe1][7])
{
  int n, i, j, k;
  double du, dv, dd, dx, dy, rf, rd, cd;
  double uc, vc, jac, dudx, dudy, dvdx, dvdy;

  dd = 0.001*pcs +1.0E-32;
  ur = vr*pkl/(prl*phl*lr +1.0E-32);
  for (i = 1; i <= nnpe; ++i) {
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    s[2][1] = lc[4+i][1];
    s[2][2] = lc[4+i][2];
    s[3][1] = lc[8+i][1];
    s[3][2] = lc[8+i][2];
    for (j = 1; j <= 2; ++j) {
      for (k = 1; k <= nnpe; ++k) {
        n = ie[e][k];
        ivo[k] = vps[n][3+j];
      }
      if (j == 1) {shape(e, nnpm, dx, dy, s, dn, ivo,
       jac, dudx, dudy);}
      if (j == 1) {shape(e, nnpm, dx, dy, s, dn, ivo,
       jac, dvdx, dvdy);}
    }

    // Droplet velocities: momentum equations
    du = ur*(vi[e][i][3] -vi[e][i][1]);
    dv = ur*(vi[e][i][4] -vi[e][i][2]);
    rd = (fabs(du) +fabs(dv))*dd/(vsc +1.0E-32) +1.0E-32;
    cd = 24.0/rd +4.73/pow(rd, 0.37) +0.00624*pow(rd, 0.38);
    uc = ur*ur*(vi[e][i][3]*dudx +vi[e][i][4]*dudy)
     /(lr +1.0E-32);
    vc = -9.8 +ur*ur*(vi[e][i][3]*dvdx +vi[e][i][4]*dvdy)
     /(lr +1E-32);
    rf = (dd*dd*24/rd)/(18*vsc*cd +1.0E-32);
    vi[e][i][3] = vi[e][i][5] +uc*rf/(ur +1.0E-32);
    vi[e][i][4] = vi[e][i][6] +vc*rf/(ur +1.0E-32);
  }
}

// Stiffness matrices: C equation (continuity / droplet flow)
void stiffdc(int e, double df, double xip[npe1][3],
 double flux[4][npe1][npe1], double aq[4][npe1][npe1],
 double rq[npe1])
{
  int m, n, i, j;
  double advc, mflow, mw;

  // Convection coefficients
  ipoint(4, e, 0.0, xip, ic);
  for (i = 1; i <= nnpe; ++i) {
    s[1][1] = lc[i][1];
    s[1][2] = lc[i][2];
    m = ie[e][i];
    mflow = vi[e][i][3]*ds[i][2] -vi[e][i][4]*ds[i][1];

    // Mass flux of droplets
    for (j = 1; j <= nnpe; ++j) {
      n = ie[e][j];
      mw = 2.0 -2.0/(1.0 +pow(fl[n][1], 0.1));
      if (ph[n][1] == 2) {mw = 0.0;}
      advc = -ic[e][1][1][i][4+j]*mflow*mw;
      flux[1][i][j] = vr*advc;
    }
  }

  // Local stiffness matrix: phase volume fraction
  for (i = 1; i <= nnpe; ++i) {
    n = ie[e][i];
    for (j = 1; j <= nnpe; ++j) {
      if (i == j) {aq[1][i][j] = vr*ds[i][3]/df -flux[1][i][j]
       +flux[1][o[i]][j];}
      if (i != j) {aq[1][i][j] = -flux[1][i][j]
       +flux[1][o[i]][j];}
    }
    rq[i] = vr*cn[n][2]*ds[i][3]/df;
  }
}

// Surface: three-phase heat, mass / momentum balances
void phase3(int nel, int nnp, int &cmf, int &cmx, int &cr,
 int nsrf, double df, int ph[nnpm][3], double fl[nnpm][3],
 double cn[nnpm][6], double ec[nnpm][7])
{
  int e, i, j, k, m, n, phj, ns, nf, rn, rg;
  double hconv, u0, g0, hf, dx, dy, vx1, vx2, cl;
  double ds, dm, di, t0, ts, tf, tm, ti, tp, tw;
  double dudx, dudy, chi, e0, beta, cw, at, bt;
  double dt, tauf, hs, taus, mr, flr, dtm, ar2;

  // Initialized problem parameters
  rg = 1; // rime (0) or rime / glaze (1) ice
  g0 = 1.0E+08; // initialized liquid water content
  u0 = 0.0; // initialized reference velocity
  dtm = df*lr*lr*prl*phl/pkl; // dimensional time step
  cr = 0; // phase convergence parameter
  cmx = 2; // number of iterations parameter
  hconv = 500.0; // Comparison with Myers (IJHMT, 1999)
  beta = 0.5; // empirical collection efficiency

  // Reference properties: air / water / ice system
  ur = vr*pkl/(prl*phl*lr +1.0E-32); // dimensional velocity
  tp = 273.0; // phase change temperature of water
  dt = 1.0; // characteristic temperature difference
  chi = 8.52; // evaporation coefficient
  e0 = 44.4; // vapor pressure constant
  cw = 4220.0; // specific heat of unfrozen water
  tf = (tp -tmin)/(tmax -tmin +1.0E-32);
  tw = tp; // reference wall temperature
  for (i = 1; i <= nsrf; ++i) {
    for (k = 1; k <= 2; ++k) {
      j = 2*(i -1) +k;
      if ((bc[j][3][1] == 0) && (bc[j][3][3] == 0) &&
       (bc[j][4][1] == 0) && (bc[j][4][3] == 0) &&
       (bc[j][2][2] != 0)) {
        tw = min(tw, tmin +(tmax -tmin)*bc[j][2][3]
         /(bc[j][2][2] +1.0E-32));
      }
    }
  }

  // Mass conservation: incoming droplets
  for (e = 1; e <= nel; ++e) {
    for (j = 1; j <= nnpe; ++j) {
      i = ie[e][j];
      tn[i][3] = 0.0;
      cn[i][1] = max(cn[i][1], 0.0);

      // Change of phase and mass in filled volumes
      if (cn[i][4] > 0.0) {cn[i][1] = cn[i][4];}
      cn[i][4] = 0.0;
      if (cmf == 1) {cmx = max(cmx, int(cn[i][1]));}
      cn[i][1] = min(cn[i][1], 2.0);
      if (cn[i][1] >= 1.0) {
        ph[i][1] = 1;
        for (k = 1; k <= nnpe; ++k) {
          n = ie[e][k];

          // Mass re-adjustment in over-filled volumes
          if (ph[n][1] != 2) {fl[n][1] = 0.0;}
          if ((cn[n][1] > 0.0) && (cn[n][1] < 0.999)) {
            ph[n][1] = 2;
            cn[n][4] = fabs(cn[i][1] -1.0);
          }
        }
        cn[i][2] = 1.0;
      }
      if (ph[i][1] != ph[i][2]) {cr = 1;}

      // Reference parameters
      g0 = max(min(1000.0*cn[i][1], g0), 1.0E-08);
      u0 = max(ur*sqrt(vps[i][1]*vps[i][1] 
       +vps[i][2]*vps[i][2]), u0);

      // Initial solid layer temperature
      if (ph[i][1] == 1) {
        for (k = 1; k <= nnpe; ++k) {
          n = ie[e][k];
          tn[i][2] = min(tn[n][2], tn[i][2]);
          if (ph[n][1] != 1) {tn[n][2] = tn[i][2];}
        }
      }
    }
  }

  // Neighbour filled volumes
  for (i = 1; i <= nnp; ++i) {
    rn = 0;
    if ((cn[i][1] > 0.999) && (ph[i][1] == 1))
     {cn[i][1] = 1.0;}
    if (nnb[i][1] == 9) {
      for (j = 2; j <= nnb[i][1] +1; ++j) {
        n = nnb[i][j];
        if ((cn[i][1] < 1) && (fl[n][1] == 0)) {rn = rn++;}
      }
      if (rn >= 8) {
        cn[i][1] = 1.0;
        fl[i][1] = 0.0;
        ph[i][1] = 1;
      }
    }
  }

  // Start heat / momentum balances
  for (e = 1; e <= nel; ++e) {
    phj = 0;
    for (j = 1; j <= nnpe; ++j) {
      i = ie[e][j];
      if (ph[i][1] == 1) {phj = 1;}
      if (ph[i][1] == 3) {phj = 3;}
      if ((ph[i][1] == 2) && (ph[ie[e][o[j]]][1] == 2)
       && (ph[ie[e][o[j+2]]][1] == 2)) {phj = 1;}
    }

    // Loop over sub-control volumes
    for (j = 1; j <= nnpe; ++j) {
      i = ie[e][j];
      n = ie[e][o[j]];
      m = ie[e][o[j+2]];
      ns = n;
      nf = m;

      // Tangential / normal directions: phase interface
      vx1 = vps[i][1]*(x[i][1] -x[n][1]) +vps[i][2]
       *(x[i][2] -x[n][2]);
      vx2 = vps[i][1]*(x[i][1] -x[m][1]) +vps[i][2]
       *(x[i][2] -x[m][2]);
      if (fabs(vx1) > fabs(vx2)) {
        ns = m;
        nf = n;
      }
      t0 = tmin +(tmax -tmin)*tn[nf][1];
      tm = max(tmin +(tmax -tmin)*tn[ns][1], tw);
      if ((vx1 < 0) || (vx2 < 0)) {
        t0 = tmin +(tmax -tmin)*tn[i][1];
        tm = max(tmin +(tmax -tmin)*tn[nf][1], tw);
      }

      // Rime ice
      if (rg == 0) {
        cn[i][5] = 1.0;
      }

      // Inside the solid
      else if (ph[i][1] == 1) {
        ec[i][2] = phs/(phl +1.0E-32);
      }

      // Solid side of phase interface with unfrozen water
      else if ((phj == 1) && (ph[i][1] == 2)) {
        if (cn[m][1] > cn[n][1]) {n = m;}
        dm = lr*sqrt(pow(x[nf][1] -x[i][1], 2.0)
         +pow(x[nf][2] -x[i][2], 2.0));
        at = min(max((hconv*(tp -t0) +0.333*pks*(tp -tm)/dm)
         /(chi*e0*(tp -t0) +0.5*pks*(tp -tm)/dm +beta*u0*g0
         *cw*dt +hconv*dt +1.0E-32), 0.0), 1.0);

        // Add conduction to heat balance (rime; cn[i][5] = 1)
        ds = lr*sqrt(pow(x[i][1] -x[n][1], 2.0)
         +pow(x[i][2] -x[n][2], 2.0));
        di = 0.5*ds*(1.0 +2.0*min(0.5, cn[n][1]));
        ti = at*tf +(1.0 -at)*tn[i][1];
        cn[i][5] = min(max(cn[i][5] +fabs((tmax -tmin)*pks
         *(tf -ti)/(di*u0*g0*pl +1.0E-32)), 0.0), 1.0);
      }

      // Air side of phase interface (incoming droplets)
      else if ((phj == 3) && (ph[i][1] == 2)) {
        ds = lr*sqrt(pow(x[ns][1] -x[i][1],2.0)
         +pow(x[ns][2] -x[i][2],2.0));
        s[1][1] = lc[12+j][1];
        s[1][2] = lc[12+j][2];
        s[2][1] = lc[4+j][1];
        s[2][2] = lc[4+j][2];
        s[3][1] = lc[8+j][1];
        s[3][2] = lc[8+j][2];
        shape(e, nnpm, dx, dy, s, dn, ivo, ar2, dudx, dudy);

        // Film momentum balance for shear-driven runback
        hf = fl[nf][1]*lr*lr*ar2/(0.5*ds +1.0E-32);
        tauf = 0.5*(0.008 +2.0E-05*(u0*hf/vsc))*prl*u0*u0;
        hs = fl[ns][1]*lr*lr*ar2/(0.5*ds +1.0E-32);
        taus = 0.5*(0.008 +2.0E-05*(u0*hs/vsc))*prl*u0*u0;
        mr = 0.5*fabs(tauf*hf*hf -taus*hs*hs)*dtm/vsc;
        flr = mr/(1000.0*lr*lr*ar2 +1.0E-32);
        fl[i][1] = max(min(fl[i][1] +flr, 1.0 -cn[i][1]), 0.0);

        // Laminar flow, smooth surface of arbitrary shape
        //hconv = 0.2926*pkl*pow(u0/(vsc*ds), 0.5);

        // Turbulent flow, smooth surface of arbitrary shape
        // hconv = 0.0287*prl*phl*u0*pow(vsc*prl*phl/pkl, -0.4)
        //  *pow(vsc/(u0*ds), 0.2);

        // Temperature rise due to latent heat released
        ts = t0 +(u0*g0*pl +0.5*hconv*u0*u0/phl)
         /(hconv +u0*g0*phl);
        tn[i][3] = min((ts -tmin)/(tmax -tmin), tf);

        // Criterion for unfrozen liquid layer (glaze ice)
        if (ts >= tp) {
          bt = min(max((u0*g0*cw*dt -0.5*g0*u0*u0*u0
           +chi*e0*dt)/(u0*g0*pl +1.0E-32), 0.0), 1.0);
          cn[i][5] = min(max((cn[i][5] +hconv*(tp -t0))
           /(u0*g0*pl +1.0E-32) +bt, 0.0), 1.0);
        }
        ec[i][1] = 0.0;
        ec[i][2] = 1.0E+08;
        ec[i][3] = 0.0;
      }

      // Droplets in the freestream
      else if (ph[i][1] == 3) {
        ec[i][2] = 1.0;
      }
      ec[i][4] = ec[i][1];
      ec[i][5] = ec[i][2];
      ec[i][6] = ec[i][3];
    }
  }

  // Volume fraction change: surface runoff
  for (i = 1; i <= nnp; ++i) {
    if ((ph[i][1] == 2) && (tn[i][3] >= tf)) {
      cl = (1.0 -cn[i][5])*(cn[i][1] -cn[i][2]);
      fl[i][1] = min(max(cl, 0.0), 1.0);
      cn[i][1] = min(max(cn[i][1] -cl, 0.0), 1.0);
    }
    cn[i][5] = 0.0;
  }
}

// EXE Input / output for F77 executable program
void filedata(int fn)
{
  char c1, c2, c3, c4, c5, c6, c7, c8;
  int i, j, k, n, nnp, nel, nx, ny, n1, n2;
  int nsrf, var1, var2, i1, i2, i3, i4;
  double faca, facb, dval1, dval2, dval3, delt, val;
  double cond, rho, cp, uref, tmn, tmx, lref, worked;

  // Input files
  ifstream file3;
  ifstream file4;
  ifstream file5;
  ifstream file6;
  ifstream file7;
  ifstream file8;
  ifstream file9;
  file3.open("c:/phases/mesh13d.dat");
  file4.open("c:/phases/samp13d.prj");
  file5.open("c:/phases/mesh13d.dat");
  file6.open("c:/phases/bc13d.dat");
  file7.open("c:/phases/ic13d.dat");
  file8.open("c:/phases/output.dat");
  file9.open("c:/phases/cust13d.dat");
  file3 >> nnp >> c1 >> nel>> c2 >> n1 >> c3 >> n2 >> c4 >>
   nx >> c5 >> ny >> c6 >> nsrf;
  file3.close();
  for (i = 1; i <= 30; ++i) {
    file4 >> val;
    worked = Adda(i, 1, 1, 1, val, nnp, nel, nx, ny, nsrf);
    if (i == 6) {delt = val;}
    if (i == 7) {lref = val;}
    if (i == 10) {uref = val;}
    if (i == 17) {tmn = val;}
    if (i == 18) {tmx = val;}
    if (i == 21) {rho = val;}
    if (i == 25) {cond = val;}

    // Set dimensionless time step
    if (i == 27) {
      cp = val;
      if (uref == 0.0) {uref = cond/(rho*cp*lref);}
      worked = Adda(10, 1, 1, 1, uref, nnp, nel, nx, ny, nsrf);
      df = delt/(lref/uref);
      worked = Adda(6, 1, 1, 1, df, nnp, nel, nx, ny, nsrf);
    }
  }

  // Read mesh file (element-node array and x,y node positions)
  file5 >> nnp >> c1 >> nel >> c2 >> n1 >> c3 >> n2 >> c4 >> nx
   >> c5 >> ny >> c6 >> nsrf;
  for (i = 1; i <= nel; ++i) {
    file5 >> ie[i][1] >> c1 >> ie[i][2] >> c2 >> ie[i][3] >> c3
     >> ie[i][4];
    for (j = 1; j <= nnpe; ++j) {
      var1 = ie[i][j];
      worked = Adda(31, i, j, var1, 1, nnp, nel, nx, ny, nsrf);
    }
  }
  for (i = 1; i <= nnp; ++i) {
    file5 >> dval1 >> c1 >> dval2;
    worked = Adda(32, i, 1, 1, dval1, nnp, nel, nx, ny, nsrf);
    worked = Adda(33, i, 1, 1, dval2, nnp, nel, nx, ny, nsrf);
  }
  for (i = 1; i <= nx-1; ++i) {
    for (j = 1; j <= ny-1; ++j) {
      file5 >> i1 >> c1 >> i2 >> c2 >> i3 >> c3 >> i4;
    }
  }

  // Read boundary condition file
  for (i = 1; i <= nsrf; ++i) {
    file6 >> var1 >> c1 >> var2;
    worked = Adda(41, i, 1, var1, 1, nnp, nel, nx, ny, nsrf);
    worked = Adda(41, i, 2, var2, 1, nnp, nel, nx, ny, nsrf);
  }
  for (i = 1; i <= nsrf; ++i) {
    file6 >> var1;
    worked = Adda(42, i, var1, 1, 1, nnp, nel, nx, ny, nsrf);
  }
  for (i = 1; i <= 4; ++i) {
    if (i == 1) {
      faca = 0.0;
      facb = 1.0;
    }
    else if (i == 2) {
      faca = tmn;
      facb = tmx;
    }
    else {
      faca = 0.0;
      facb = uref;
    }
    for (k = 1; k <= 2*nsrf; ++k) {
      file6 >> dval1 >> c1 >> dval2 >> c2 >> val >> c3 >> c4
       >> c5 >> c6 >> c7 >> c8 >> var1;
      dval3 = (val -dval2*faca)/(facb -faca +1.0E-32);
      worked = Adda(43, i, 1, k, dval1, nnp, nel, nx, ny, nsrf);
      worked = Adda(43, i, 2, k, dval2, nnp, nel, nx, ny, nsrf);
      worked = Adda(43, i, 3, k, dval3, nnp, nel, nx, ny, nsrf);
    }
  }

  // Read initial condition file
  for (n = 1; n <= 4; ++n) {
    if (n == 1) {
      faca = 1.0;
      facb = 0.0;
    }
    else if (n == 2) {
      faca = 1.0/(tmx -tmn +1.0E-32);
      facb = -tmin/(tmx -tmn +1.0E-32);
    }
    else if ((n == 3) || (n == 4)) {
      faca = 1.0/(uref +1.0E-32);
      facb = 0.0;
    }
    for (i = 1; i <= nel; ++i) {
      file7 >> var1 >> c1 >> val >> c2 >> val >> c3 >> val
       >> c4 >> val;
      val = faca*val +facb;
      worked = Adda(50 +n, i, 1, 1, val, nnp, nel, nx, ny, nsrf);
    }
  }
  worked = Adda(60, 1, 1, 1, 1, nnp, nel, nx, ny, nsrf);

  // Display output
  //for (n = 1; n <= fn; ++n) {
  //  worked = Adda(0, n, 1, 1, 1, nnp, nel, nx, ny, nsrf);
  //  cout << "step complete: " << n << endl;
  //  for (i = 1; i <= nnp; ++i) {
  //    worked = Adda(200, i, 1, 1, 1, nnp, nel, nx, ny, nsrf);
  //		if (i == 89) {cout << " T = " << worked << endl;}
  //  }
  //  worked = Adda(100, 10, 1, 1, 1, nnp, nel, nx, ny, nsrf);
  //}

	cout << " 1-D Icing Validation Problem " << endl;
  double cs, esn, u0, c0, l0, es, ti, xn, delx, ta, ts, bx;
  faca = 1.0/(tmax -tmin +1.0E-32);
  facb = -tmin/(tmax -tmin +1.0E-32);
  for (n = 1; n <= fn; ++n) {
    worked = Adda(0, n, 1, 1, 1, nnp, nel, nx, ny, nsrf);
    j = 0;
    cs = 0.0*1.0;
    esn = 0.0;
    u0 = 90.0 +0*0.1;
    c0 = 1.0E-06 +0*0.1;
    l0 = 0.01 +0*1.0;
    for (i = 1; i <= nnp; ++i) {
      if (fmod(i, ny) == 1) {
        esn = worked;
        worked = Adda(200, i, 1, 1, 1, nnp, nel, nx, ny, nsrf);
        ti = Adda(201, i, 1, 1, 1, nnp, nel, nx, ny, nsrf);
        j = j++;
        xn = 0.01*(j -1)/(nx -1);
        delx = 0.01/(nx -1);
        if ((worked == 1.0) && (esn < 1.0))
         {cs = 1000.0*(0.01 -(xn -0.5*delx -esn*delx));}
      }
    }
    if (cs == 0.0) {cs = 1000.0*worked*0.5*delx;}
    es = 1000.0*(n*delt -l0/u0)*u0*c0;
    //if (n == 10) {cout << " ... " << endl;}
    cout << n*delt << " " << es << " " << cs << endl;
    //if ((n < 6) || (n > 12)) {cout << n*delt << " "
    // << es << " " << cs << " " << (cs-co)/delt << endl;}
    //co = cs;
  } 
  ta = tmin +(tmax -tmin)*tn[1][1] -273.0;
  ts = tmin +(tmax -tmin)*tn[nnp][1] -273.0;
  if ((ta == -1.0) && (ts == -1.0)) {bx = 1.13;}
  if ((ta == -1.0) && (ts == -10.0)) {bx = 3.40;}
  if ((ta == -10.0) && (ts == -1.0)) {bx = 2.07;}
  if ((ta == -10.0) && (ts == -10.0)) {bx = 4.75;}
  cout << " Ta = " << ta << " Ts = " << ts <<
   " Analytic (Myers): B(100) = " << bx << endl;

  file4.close();
  file5.close();
  file6.close();
  file7.close();
  file8.close();
  file9.close();
}
