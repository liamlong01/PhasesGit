
void iterC() {

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
			cdif = cdif + fabs(cn[i][1] - cn[i][3]) / nnp;
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
		if (cmf >= cmx) { cr = 0; }
	}

}