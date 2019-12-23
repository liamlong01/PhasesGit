

//This header performs necessary setup so that Pyfiledata and PyAdda can be called from other programs such as python

#ifdef PHASES_EXPORTS
	#define PHASES_API  extern "C" __declspec(dllexport) 
#else
	#define PHASES_API __declspec(dllimport) 
#endif
#include <iostream>
using namespace std;


PHASES_API void Pyfiledata(int steps, char* meshDir, char* bcDir, char* icDir, char* outputDir, char* prjDir);

PHASES_API double PyAdda(int w1, int w2, int w3, int w4, double w5, int nnp, int nel, int nx, int ny, int nsrf);