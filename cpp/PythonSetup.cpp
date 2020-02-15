


#include "constInit.h"
#include "functInit.h"
#include "variaInit.h"
#include <iostream> 
using namespace std;
#include "PythonSetup.h"


//These are simple wrapper functions that can be called by Python in order to run the phases software

void Pyfiledata(int steps, char* meshDir, char* bcDir, char* icDir, char* outputDir, char* prjDir)
{

	filedata(steps, meshDir, bcDir, icDir, outputDir, prjDir);

}

double PyAdda(int w1, int w2, int w3, int w4, double w5, int nnp,
	int nel, int nx, int ny, int nsrf)
{
	return Adda(w1, w2, w3, w4, w5, nnp, nel, nx, ny, nsrf);
}

double PyCtrl(int nnp, int nel, int nsrf) {

	// Options: se = 1(C), 2(T), 3(UVP), 12(CT)
	// Options: se = 13(CUVP), 23(TUVP), 4(CTUVP)
	int i, j, n, cr, tr, tph, cmf, cmx, vk, np3;
	double cdif, cphs, cphl, tphs, tphl;
	double vdif, tdif, ste, rlx, hbar;

	// Store field variables
	ste = phl * dtr / (pl + 1.0E-32);
	np3 = nnp * 3;
	tdif = 0.0;

	PyCtrlInit(nnp);


	// Start C - T - V iteration sequence
	for (vk = 1; vk <= tvk; ++vk) {
	
		if ((se == 1) || (se == 12) || (se == 13) || (se == 4)) {
			
			PyCtrlC_init(nnp);
			cmx = 2;
			cr = 0;

			for (cmf = 1; cmf <= tk; ++cmf) {
				
				cdif = PyCtrlC(nel,nsrf,nnp,cmf);
				// Quit early for small residual
				if ((cr == 0) && (cdif < tol)) { break; }
			}
		}

		// Start temperature - phase iteration algorithm
		if ((se == 2) || (se == 12) || (se == 23) || (se == 4)) {
			for (tph = 1; tph <= tk; ++tph) {
																	// calls phase
				tdif = PyCtrlT(nel, nsrf, nnp);


				// Quit early if tentative phase matches new phase
				if ((tr == 0) || (tdif < tol)) { break; }
			}

		
		}

		//Phase-Temperature Consistency
		PyPhaseTempCheck(nnp);

		// Start U - V - P solution iterations
		if ((se == 3) || (se == 13) || (se == 23) || (se == 4)) {
			for (n = 1; n <= tk; ++n) {

				vdif = PyCtrlUVP(nel, nsrf, nnp, aq, rq);

				if (vdif <= tol) { break; }
			}
		}


		if ((tdif < tol) && (vdif < tol)) { break; }
	}	


	//entropy extensions
	if (ao == 1) {
		entropy(df, nel, nnp, nsrf, vps);														// calls entropy
		dsmooth(nel, nnp, vps);																	// calls dsmooth
	}

}

void PyPhaseTempCheck(int nnp) {
	double tphs, tphl, cphs, cprl, cphl, hbar, ste;

	ste = phl * dtr / (pl + 1.0E-32);
	if (ao != 2) {
		for (int i = 1; i <= nnp; ++i) {
			tphs = tmlt - eps - cn[i][1] * (tmlt - eps - tsol)
				/ (pcs + 1.0E-32);
			tphl = tmlt + eps - cn[i][1] * (tmlt + eps - tsol)
				/ (pcl + 1.0E-32);
			cphs = pcs * (tmlt - eps - tn[i][1]) / (tmlt - eps - tsol);
			cphl = pcl * (tmlt + eps - tn[i][1]) / (tmlt + eps - tsol);
			hbar = 0.5 * (phs + phl);

			// Solid phase reference states
			if (tn[i][1] <= tphs) {
				ph[i][1] = 1;
				fl[i][1] = 0.0;
				ec[i][1] = 0.0;
				ec[i][2] = phs / (phl + 1.0E-32);
				ec[i][3] = 0.0;

			}

			// Liquid phase reference states
			else if (tn[i][1] >= tphl) {
				ph[i][1] = 3;
				fl[i][1] = 1.0;
				ec[i][1] = hbar * (tphl - tphs) / (phl + 1.0E-32)
					+ (phs * tphs) / (phl + 1.0E-32) + 1.0 / (ste + 1.0E-32);
				ec[i][2] = phl / (phl + 1.0E-32);
				ec[i][3] = tphl;
			}

			// Melt phase reference states
			else {
				ph[i][1] = 2;
				ec[i][1] = (phs * tphs) / (phl + 1.0E-32);
				ec[i][2] = hbar / (phl + 1.0E-32) + 1.0 / (ste * (tphl - tphs)
					+ 1.0E-32);
				ec[i][3] = tphs;
				fl[i][1] = (tn[i][1] - tphs) / (tphl - tphs);
			}
		}
	}






}
void PyCtrlInit(int nnp) {

	for (int i = 1; i <= nnp; ++i) {
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


}


void PyCtrlC_init(int nnp) {
	for (int i = 1; i <= nnp; ++i) {
		cn[i][2] = cn[i][3];
		cn[i][4] = 0.0;
		fl[i][1] = fl[i][2];
	}


}

double PyCtrlC(int nel, int nsrf, int nnp, int cmf) {

	double cdif = 0.0;
	int cmx = 2;
	int cr = 0;

	// C equation: solid / liquid phase change
	if (ao != 2) {
		nullmx(r, z, c, aq);															// calls nullmx
		assmbl(1, 1, 1, nel, nelm, band3, df, fo, nn, aq, c, rq, r);					// calls assmbl
		bndry(1, 1, 1, nsrf, fo, bc, c, r);											// calls bndry
		genul(1, c, nym, nnp, nnpm, band3);											// calls genul	
		forbak(1, c, z, r, nym, nnp, nnpm, band3);									// calls forbak

		// Solution residual
		for (int i = 1; i <= nnp; ++i) {
			cn[i][1] = z[i];
			cdif = cdif + fabs(cn[i][1] - cn[i][3]) / nnp;
			cn[i][3] = cn[i][1];
		}
	}

	// C equation: droplet flow
	else if (ao == 2) {
		for (int i = 1; i <= nnp; ++i) {
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
		for (int j = 1; j <= nnp; ++j) {
			cn[j][1] = z[j];
		}

		// Surface: three-phase heat, mass / momentum balances
		phase3(nel, nnp, cmf, cmx, cr, nsrf, df, ph, fl, cn, ec);							// calls phase3
		if(cmf >= cmx) { cr = 0; }
	}

	return cdif;

}

double PyCtrlT(int nel, int nsrf, int nnp) {

	double tdif;
	int tr;
	for (int i = 1; i <= nnp; ++i) {
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
	for (int i = 1; i <= nnp; ++i) {
		tn[i][1] = z[i];
		tdif = tdif + fabs(tn[i][1] - tn[i][3]) / nnp;


		/* double temp;
		temp = z[i] * (tmax - tmin) + tmin;
		cout << temp << endl;

		cout << temp << endl; */
	}
	// Apply phase interface motion rules

	if (ao != 2) { phase(nnp, ph, fl, ec, tr); }

	return tdif;

}



//TODO how to intialize PYCtrlUVP with array arguments
// aq must persist across the three ctrl loops!!!
// easy (currently used for other vars) solution -> put in variaInit.h
double PyCtrlUVP(int nel, int nsrf, int nnp) {

	//vdif rlx are origninally  defined in ctrl.cpp but only used in UVP solver
	double vdif, rlx;
	int np3 = 3 * nnp;


	// Element assembly: x momentum, y momentum, mass equations

	// 1. nullmx needs r z c aq
		
	/* - r: declared variainit.h as double r[3*nnpm], nnpm = max nodal points (constinit.h)
		//r is modified by assembl
	   - z: input argument to ctrl.cpp, from adda this then defined in variaInit.h (as z), o to use global z as mem space ??
	   -c : declared in variainit.h as double c[3*nnpm][6*nym+11], nym = max number of rows
		- aq: local in cntrl function defined as  static double aq[4][npe1][npe1], an argument here;
	*/
	nullmx(r, z, c, aq);	
		
		
	// assembl calls
		//so far undefined:
		// nel input to Adda
		// nelm in constInit.h
		// band3 declared in variainit as an int, initialized in Adda as value 6*ny+11
		// df, size of timestep, arg of ctrl, if (w1 == 6) {df = w5; } in Adda, declared as int in variaInit.h
		// fo , number of timesteps, variainti.h declared, initialized in Adda on w1 = 0!!! calls
		// nn is in variaint.h
	assmbl(3, 1, 3, nel, nelm, band3, df, fo, nn, aq, c, rq, r);								// calls assembl			
	assmbl(4, 2, 3, nel, nelm, band3, df, fo, nn, aq, c, rq, r);								// calls assembl
	assmbl(5, 3, 3, nel, nelm, band3, df, fo, nn, aq, c, rq, r);								// calls assembl

	// Boundary conditions, U-V-P direct simultaneous solution

	//nsrf -
	//bc - variaInit.h
	//nsrf -> input to Adda function pased to ctrl

	bndry(3, 1, 3, nsrf, fo, bc, c, r);															// calls bndry
	bndry(4, 2, 3, nsrf, fo, bc, c, r);															// calls bndry
	bndry(4, 3, 3, nsrf, fo, bc, c, r);	
	
	// calls bndry
	//np3 - 3*nnp, scopeded locally in ctrl.cpp, nnp is an input to Adda
	genul(3, c, nym, np3, nnpm, band3);															// calls genul
	forbak(3, c, z, r, nym, np3, nnpm, band3);													// calls forbak
	uvcor(nel, vi);																				// calls uvcor
	for (int i = 1; i <= nnp; ++i) {
		vps[i][1] = z[3 * (i - 1) + 1];
		vps[i][2] = z[3 * (i - 1) + 2];
		vps[i][3] = z[3 * (i - 1) + 3];
	}

	// Velocity residual
	vdif = 0.0;
	rlx = 0.5;
	for (int i = 1; i <= nnp; ++i) {
		vdif = vdif + fabs(vps[i][6] - vps[i][1]) / nnp;
		vps[i][6] = vps[i][1] + rlx * (vps[i][6] - vps[i][1]);
		vps[i][7] = vps[i][2] + rlx * (vps[i][7] - vps[i][2]);
	}

	return vdif;
	//TOD: move this check to python/rxternctrl loop
	//if (vdif <= tol) { goto cnt8200; }
}




