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
//-- EXE Input / output for C++ executable program --//
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
 
  /*
  file3.open("C:\\Users\\Peter\\Desktop\\Resource\\mesh2.dat"); // added our own computer's directory
  file4.open("C:\\Users\\Peter\\Desktop\\Resource\\sample2.prj");
  file5.open("C:\\Users\\Peter\\Desktop\\Resource\\mesh2.dat");
  file6.open("C:\\Users\\Peter\\Desktop\\Resource\\bc2.dat");
  file7.open("C:\\Users\\Peter\\Desktop\\Resource\\ic2.dat");
  file8.open("C:\\Users\\Peter\\Desktop\\Resource\\output2.dat");
  file9.open("C:\\Users\\Peter\\Desktop\\Resource\\custom2.dat");
 */
  
  file3.open("C:\\Users\\Peter\\Desktop\\CFD\\Spring 2016\\Resources\\prj2\\mesh2.dat"); // added our own computer's directory
  file4.open("C:\\Users\\Peter\\Desktop\\CFD\\Spring 2016\\Resources\\prj2\\sample2.prj");
  file5.open("C:\\Users\\Peter\\Desktop\\CFD\\Spring 2016\\Resources\\prj2\\mesh2.dat");
  file6.open("C:\\Users\\Peter\\Desktop\\CFD\\Spring 2016\\Resources\\prj2\\bc2.dat");
  file7.open("C:\\Users\\Peter\\Desktop\\CFD\\Spring 2016\\Resources\\prj2\\ic2.dat");
  file8.open("C:\\Users\\Peter\\Desktop\\CFD\\Spring 2016\\Resources\\prj2\\output2.dat");
  file9.open("C:\\Users\\Peter\\Desktop\\CFD\\Spring 2016\\Resources\\prj2\\custom2.dat");
  

  /*
	file3.open("C:/Users/nw8474/Desktop/phases/dat/mesh3.dat"); // added our own computer's directory
  file4.open("C:/Users/nw8474/Desktop/phases/dat/sample3.prj");
  file5.open("C:/Users/nw8474/Desktop/phases/dat/mesh3.dat");
  file6.open("C:/Users/nw8474/Desktop/phases/dat/bc3.dat");
  file7.open("C:/Users/nw8474/Desktop/phases/dat/ic3.dat");
  file8.open("C:/Users/nw8474/Desktop/phases/dat/output31.dat");
  file9.open("C:/Users/nw8474/Desktop/phases/dat/custom3.dat");
  */

  

	// Output files
	ofstream file10;
	//file10.open("c:/phases/output1.dat");
       file10.open("C:\\Users\\Peter\\Desktop\\Resource\\output2.dat");
	   ofstream file11;
	   file11.open("C:\\Users\\Peter\\Desktop\\Resource\\vps2.dat");

  file3 >> nnp >> c1 >> nel >> c2 >> n1 >> c3 >> n2 >> c4 >>
   nx >> c5 >> ny >> c6 >> nsrf;
  file3.close();
  for (i = 1; i <= 30; ++i) {
    file4 >> val;
	//cout << val << endl; 
    worked = Adda(i, 1, 1, 1, val, nnp, nel, nx, ny, nsrf);													// calls Adda
	//cout << worked << endl; 

	if (i == 6) { delt = val;}
    if (i == 7) {lref = val;}
    if (i == 10) {uref = val;}
    if (i == 17) {tmn = val;}
	if (i == 18) { tmx = val; /*cout << tmx << endl;*/ }
	if (i == 21) {rho = val; /*cout << rho << endl;*/ }
	if (i == 25) { cond = val; /*cout << cond << endl; */}
	
	// cout << val << endl;

    // Set dimensionless time step
    if (i == 27) {
      cp = val;
      if (uref == 0.0) {uref = cond/(rho*cp*lref);}
      worked = Adda(10, 1, 1, 1, uref, nnp, nel, nx, ny, nsrf);												// calls Adda
	  df = delt/(lref/uref);
      worked = Adda(6, 1, 1, 1, df, nnp, nel, nx, ny, nsrf);												// calls Adda
	  // cout << df << endl;
	  //cout << uref << endl;
	}

  }

  // Read mesh file (element-node array and x,y node positions)
  file5 >> nnp >> c1 >> nel >> c2 >> n1 >> c3 >> n2 >> c4 >> nx
   >> c5 >> ny >> c6 >> nsrf;
  /*cout << nnp << " " << c1 << " " << nel << " " << c2 
	   << n1 << " " << c3 << " " << n2 << " " << c4 << " " << nx
	   << c5 << " " << ny << " " << c6 << " " << nsrf;*/

  for (i = 1; i <= nel; ++i) {
    file5 >> ie[i][1] >> c1 >> ie[i][2] >> c2 >> ie[i][3] >> c3
     >> ie[i][4];

	//cout << ie[i][1] << c1 << ie[i][2] << c2 << ie[i][3] << c3 << ie[i][4] << endl;


    for (j = 1; j <= nnpe; ++j) {
      var1 = ie[i][j];
      worked = Adda(31, i, j, var1, 1, nnp, nel, nx, ny, nsrf);												// calls Adda
	 // cout << i << " " << j << " " << var1 << endl;	
	  //cout << '$' << worked << endl;
    }
  }
  for (i = 1; i <= nnp; ++i) {
    file5 >> dval1 >> c1 >> dval2; 
		//cout<<"x= "<< dval1 <<" y= "<< dval2<< endl;
    worked = Adda(32, i, 1, 1, dval1, nnp, nel, nx, ny, nsrf);												// calls Adda	
    worked = Adda(33, i, 1, 1, dval2, nnp, nel, nx, ny, nsrf);												// calls Adda
  }
  for (i = 1; i <= nx-1; ++i) {
    for (j = 1; j <= ny-1; ++j) {
      file5 >> i1 >> c1 >> i2 >> c2 >> i3 >> c3 >> i4;

	  //cout << i << " " << j << "  " << i1 << c1 << i2 << c2 << i3 << c3 << i4 << endl;
    }
  }

  // Read boundary condition file
  for (i = 1; i <= nsrf; ++i) {
    file6 >> var1 >> c1 >> var2;

	// cout << var1 << c1 << var2 << endl; 

    worked = Adda(41, i, 1, var1, 1, nnp, nel, nx, ny, nsrf);												// calls Adda
    worked = Adda(41, i, 2, var2, 1, nnp, nel, nx, ny, nsrf);												// calls Adda

	//cout << worked << endl;

  }
  for (i = 1; i <= nsrf; ++i) {
    file6 >> var1;

	//cout << var1 << endl;

    worked = Adda(42, i, var1, 1, 1, nnp, nel, nx, ny, nsrf);												// calls Adda
	
	//cout << '$' << worked << '$' << endl;

  }
  for (i = 1; i <= 4; ++i) {
    if (i == 1) {
      faca = 0.0;
      facb = 1.0;

    }
    else if (i == 2) {
      faca = tmn;
      facb = tmx;

	//  cout << tmn << "&" << tmx << endl;
    }
    else { 
      faca = 0.0;
      facb = uref;
    }
    for (k = 1; k <= 2*nsrf; ++k) {
      file6 >> dval1 >> c1 >> dval2 >> c2 >> val >> c3 >> c4
       >> c5 >> c6 >> c7 >> c8 >> var1;
      dval3 = (val -dval2*faca)/(facb -faca +1.0E-32);
      worked = Adda(43, i, 1, k, dval1, nnp, nel, nx, ny, nsrf);											// calls Adda
      worked = Adda(43, i, 2, k, dval2, nnp, nel, nx, ny, nsrf);											// calls Adda
	  worked = Adda(43, i, 3, k, dval3, nnp, nel, nx, ny, nsrf);											// calls Adda
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
      worked = Adda(50 +n, i, 1, 1, val, nnp, nel, nx, ny, nsrf);												// calls Adda

	  // cout << var1 << c1 << val << c2 << val << c3 << val << c4 << val << c5 << worked << endl;
    }

	//cout << faca << "&" << facb << endl;
  }
  worked = Adda(60, 1, 1, 1, 1, nnp, nel, nx, ny, nsrf);														// calls Adda
   
  // Display output 
  for (n = 1; n <= fn; ++n) {
    worked = Adda(0, n, 1, 1, 1, nnp, nel, nx, ny, nsrf);														// calls Adda
    cout << "step complete: " << n << endl;
    for (i = 1; i <= nnp; ++i) {
		worked = Adda(200, i, 1, 1, 1, nnp, nel, nx, ny, nsrf);													// calls Adda
		double temperature; 
		temperature = worked*(tmax - tmin) + tmin;					// worked = (temperature-tmin)/(tmax-tmin) to give dimensional temperatures 
		double uvelocity;
		uvelocity = vps[i][1] * ur;
		double vvelocity;
		vvelocity = vps[i][2] * ur; 
		double wvelocity;
		wvelocity = vps[i][3] * ur;
	  
	  /* for (int i = 1; i <= nel; i++) {
		  for (int j = 1; j <= 4; j++) {
		  
			  cout << i << "," << j << "," << tn[i][j]<< "," << temperature << endl;
		  }
	  }*/

	  /*
	  for (i = 1; i <= nnp; i++) {
		  file11 << i << "," << vps[i][1] << "," << vps[i][2] << "," << vps[i][3] << "," << vps[i][4] << "," << vps[i][5] << "," << vps[i][6] << "," << vps[i][7] << endl;

	  }
	  */

	  double pressure;
	  pressure = vps[i][3] * (prl*ur*ur);

		if (i == 10) {cout << " T = " << temperature << endl;}
			// Optional print to ouput file 
		file10 << i << "," << temperature << "," << cn[i][6] << "," << rho << "," << uvelocity << "," << vvelocity << "," << wvelocity << "," << '0' << "," << fl[i][1] << "," << ph[i][3] << endl;
	 }

    worked = Adda(100, 10, 1, 1, 1, nnp, nel, nx, ny, nsrf);													// calls Adda
}

  // Display output (samp13d.prj (1-D icing problem))
  //cout << "1-D Icing Problem " << endl;
  //double cs, esn, u0, c0, l0, es, ti, xn, delx, ta, ts, bx;
  //faca = 1.0/(tmax -tmin +1.0E-32);
  //facb = -tmin/(tmax -tmin +1.0E-32);
  //for (n = 1; n <= fn; ++n) {
  //  worked = Adda(0, n, 1, 1, 1, nnp, nel, nx, ny, nsrf);
  //  j = 0;
  //  cs = 0.0*1.0;
  //  esn = 0.0;
  //  u0 = 90.0 +0*0.1;
  //  c0 = 1.0E-06 +0*0.1;
  //  l0 = 0.01 +0*1.0;
  //  delx = 0.01/(nx -1);
  //  for (i = 1; i <= nnp; ++i) {
  //    if (fmod(i, ny) == 1) {
  //      esn = worked;
  //      worked = Adda(200, i, 1, 1, 1, nnp, nel, nx, ny, nsrf);
  //      ti = Adda(201, i, 1, 1, 1, nnp, nel, nx, ny, nsrf);
  //      j = j++;
  //      xn = 0.01*(j -1)/(nx -1);
  //      if ((worked == 1.0) && (esn < 1.0))
  //       {cs = 1000.0*(0.01 -(xn -0.5*delx -esn*delx));}
  //      if ((worked < 1.0) && (i == nnp+1-ny))
  //       {cs = worked*(1000.0*1.5*delx);}
  //    }
  //  }
  //  if (cs == 0.0) {cs = 1000.0*worked*0.5*delx;}
  //  es = 1000.0*(n*delt -l0/u0)*u0*c0;
  //  cout << n*delt << " " << es << " " << cs << endl;
  //} 
  //ta = tmin +(tmax -tmin)*tn[1][1] -273.0;
  //ts = tmin +(tmax -tmin)*tn[nnp][1] -273.0;
  //if ((ta == -1.0) && (ts == -1.0)) {bx = 1.13;}
  //if ((ta == -1.0) && (ts == -10.0)) {bx = 3.40;}
  //if ((ta == -10.0) && (ts == -1.0)) {bx = 2.07;}
  //if ((ta == -10.0) && (ts == -10.0)) {bx = 4.75;}
  //cout << "Ta = " << ta << " Ts = " << ts
  // << " Analytic: B(100) = " << bx << endl;

  file4.close();
  file5.close();
  file6.close();
  file7.close();
  file8.close();
  file9.close();
}