void iterT(){
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
		tdif = tdif + fabs(tn[i][1] - tn[i][3]) / nnp;


		/* double temp;
		temp = z[i] * (tmax - tmin) + tmin;
		cout << temp << endl;

		cout << temp << endl; */
	}

	// Apply phase interface motion rules
	if (ao != 2) { phase(nnp, ph, fl, ec, tr); }

}