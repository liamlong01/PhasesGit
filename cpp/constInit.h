// -------------------------------------------------
// C++ PHASES - Copyright by G. F. Naterer
// Coding format: (i) // comments, (ii) column < 65,
// (ii) line continued -> one column indentation,
// (iv) start block -> two column indentation,
// (v) code optimization: set maximum speed options
// -------------------------------------------------

// This Header file contains all constants with initialization .
// Any cpp file needing any of the constants needs
// to put -- #include "constInit.h" -- at the top of the code.

const int nnpe = 4; // number of nodal points per element
const int npe1 = nnpe + 1; // C++ arrays start at 0, add one
const int nnpm = 6601; // maximum number of nodal points
const int nelm = 6410; // maximum number of elements
const int nsrfm = 1000; // maximum number of boundary surfaces
const int nym = 100; // maximum number of rows
const double k0 = 1.0E+06; // permeability parameter
const double eps = 0.001; // energy equation of state
const int tph = 0; 
const int fn = 10;
const int stol = 99;

const int o[12] = {0,4,1,2,3,4,1,2,2,1,3,4};
const double lc[17][3] = {0.0,0.0,0.0, 0.0,0.0,0.5,
 0.0,-0.5,0.0, 0.0,0.0,-0.5, 0.0,0.5,0.0, 0.0,0.0,0.0,
 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,1.0,
 0.0,-1.0,0.0, 0.0,0.0,-1.0, 0.0,1.0,0.0, 0.0,0.5,0.5,
 0.0,-0.5,0.5, 0.0,-0.5,-0.5, 0.0,0.5,-0.5};
const double nc[25][3] = {0.0,0.0,0.0, 0.0,0.5,1.0,
 0.0,-0.5,1.0, 0.0,-1.0,0.5, 0.0,-1.0,-0.5, 0.0,-0.5,-1.0,
 0.0,0.5,-1.0, 0.0,1.0,-0.5, 0.0,1.0,0.5, 0.0,1.0,1.0,
 0.0,0.0,1.0, 0.0,-1.0,1.0, 0.0,-1.0,0.0, 0.0,-1.0,-1.0,
 0.0,0.0,-1.0, 0.0,1.0,-1.0, 0.0,1.0,0.0, 0.0,0.0,1.0,
 0.0,-1.0,1.0, 0.0,-1.0,0.0, 0.0,-1.0,-1.0, 0.0,0.0,-1.0,
 0.0,1.0,-1.0, 0.0,1.0,0.0, 0.0,1.0,1.0};