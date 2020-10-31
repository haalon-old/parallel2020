#include <math.h>
#include "analytic.h"

// using Real = double

double u_analytical(double l_x, double l_y, double l_z, double x, double y, double z, double t)
{
	double a_t = PI * sqrt(1/(l_x*l_x) + 1/(l_y*l_y) + 4/(l_z*l_z));
	return sin(PI*x/l_x) * sin(PI*y/l_y) * sin(2*PI*z/l_z) * cos(a_t*t);
}

double phi(double l_x, double l_y, double l_z, double x, double y, double z)
{
	// In general:
	// phi == u_analytical(l_x,l_y,l_z,x,y,z,0);
	return sin(PI*x/l_x) * sin(PI*y/l_y) * sin(2*PI*z/l_z);
}