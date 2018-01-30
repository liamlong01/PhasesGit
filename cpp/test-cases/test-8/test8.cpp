
#include "stdafx.h"
#include <iostream>
using namespace std;
void stiffv(double);

// run in 'Release' mode
// Debug - Start without Debugging
void main(int, char**)
{
  double v1; 
  v1 = 1.0;
  stiffv(v1);
  cout << v1 << endl;
  getchar();
}

void stiffv(double v1)
{
  v1 = v1 +1.0;
}