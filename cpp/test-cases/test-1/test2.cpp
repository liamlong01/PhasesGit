
#include <iostream>
using namespace std;

void stiffu(double &v1, double ar2[2][2]);
//void stiffv(double &v1, double ar2[2][2]);

double mainline(double x1)
{
	double v1;
	v1 = 1.0;
	double ar2[2][2] = {1.0, 1.0, 1.0, 1.0};
	stiffu(v1, ar2);
//stiffv(v1, ar2);
//cout << ar2[1][1] << " " << v1 << endl;
  double result = 3.0*x1;
  return result;
}

void stiffu(double &v1, double ar2[2][2])
{
	ar2[1][1] = ar2[1][1] + 1.0;
	v1 = v1 + 1;
}
