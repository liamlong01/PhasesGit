// -------------------------------------------------
// C++ PHASES - Copyright by G. F. Naterer
// Coding format: (i) // comments, (ii) column < 65,
// (ii) line continued -> one column indentation,
// (iv) start block -> two column indentation,
// (v) code optimization: set maximum speed options
// -------------------------------------------------

// This Header file contains all the function declarations.
// Any cpp file that needs any function defined below, needs
// to put -- #include "functInit.h" -- at the top of the code.

double Adda(int w1, int w2, int w3, int w4, double w5, int nnp,
 int nel, int nx, int ny, int nsrf);
void assmbl(int v, int p, int m, int nel, int nelm,
 int band3, double df, double fo, double nn[npe1][npe1],
 double aq[4][npe1][npe1], double c[3*nnpm][6*nym+11],
 double rq[npe1], double r[3*nnpm]);
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
 double c[3*nnpm][6*nym+11], double aq[4][npe1][npe1]);
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