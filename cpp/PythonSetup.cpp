


#include "constInit.h"
#include "functInit.h"
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


