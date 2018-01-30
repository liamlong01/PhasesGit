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
//-- Main control sequence of inter-equation iterations --//
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
  static double aq[4][npe1][npe1], rq[npe1];

  // Store field variables
  ste = phl*dtr/(pl +1.0E-32);
  np3 = nnp*3;
  tdif = 0.0;
  for (i = 1; i <= nnp; ++i) {
											            	 // cout << tn[i][2] << "," << tn[i][1] << "," << tn[i][3] << endl;
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
   //cout<<" velocity loop....."<<vk<<endl;
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
          nullmx(r, z, c, aq);															// calls nullmx
          assmbl(1, 1, 1, nel, nelm, band3, df, fo, nn, aq, c, rq, r);					// calls assmbl
          bndry(1, 1, 1, nsrf, fo, bc, c, r);											// calls bndry
          genul(1, c, nym, nnp, nnpm, band3);											// calls genul	
          forbak(1, c, z, r, nym, nnp, nnpm, band3);									// calls forbak

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
          nullmx(r, z, c, aq);																// calls nullmx
          assmbl(6, 1, 1, nel, nelm, band3, df, fo, nn, aq, c, rq, r);						// calls assmbl						
          assmbl(7, 1, 1, nel, nelm, band3, df, fo, nn, aq, c, rq, r);						// calls assembl
          bndry(1, 1, 1, nsrf, fo, bc, c, r);												// calls bndry
          genul(1, c, nym, nnp, nnpm, band3);												// calls genul
          forbak(1, c, z, r, nym, nnp, nnpm, band3);										// calls forbak
          for (j = 1; j <= nnp; ++j) {
            cn[j][1] = z[j];
          }

          // Surface: three-phase heat, mass / momentum balances
          phase3(nel, nnp, cmf, cmx, cr, nsrf, df, ph, fl, cn, ec);							// calls phase3
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
        nullmx(r, z, c, aq);																	// calls nullmx
        assmbl(2, 1, 1, nel, nelm, band3, df, fo, nn, aq, c, rq, r);							// calls assembl 
        bndry(2, 1, 1, nsrf, fo, bc, c, r);														// calls bndry
        genul(1, c, nym, nnp, nnpm, band3);														// calls genul
        forbak(1, c, z, r, nym, nnp, nnpm, band3);												// calls forbak


      
        // Temperature solution residual
        tdif = 0.0;
        tr = 0;
        for (i = 1; i <= nnp; ++i) {
          tn[i][1] = z[i];
          tdif = tdif +fabs(tn[i][1] -tn[i][3])/nnp;
		  
		  
		  /* double temp;
		  temp = z[i] * (tmax - tmin) + tmin;
		  cout << temp << endl;

		  cout << temp << endl; */
        }

        // Apply phase interface motion rules
        if (ao != 2) {phase(nnp, ph, fl, ec, tr);}													// calls phase
        
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
        nullmx(r, z, c, aq);																		// calls nullmx
        assmbl(3, 1, 3, nel, nelm, band3, df, fo, nn, aq, c, rq, r);								// calls assembl			
        assmbl(4, 2, 3, nel, nelm, band3, df, fo, nn, aq, c, rq, r);								// calls assembl
        assmbl(5, 3, 3, nel, nelm, band3, df, fo, nn, aq, c, rq, r);								// calls assembl

        // Boundary conditions, U-V-P direct simultaneous solution
        bndry(3, 1, 3, nsrf, fo, bc, c, r);															// calls bndry
        bndry(4, 2, 3, nsrf, fo, bc, c, r);															// calls bndry
        bndry(4, 3, 3, nsrf, fo, bc, c, r);															// calls bndry
        genul(3, c, nym, np3, nnpm, band3);															// calls genul
        forbak(3, c, z, r, nym, np3, nnpm, band3);													// calls forbak
        uvcor(nel, vi);																				// calls uvcor
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
    entropy(df, nel, nnp, nsrf, vps);														// calls entropy
    dsmooth(nel, nnp, vps);																	// calls dsmooth
  }
}