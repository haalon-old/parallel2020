#pragma once
#include <cmath>

#define PI 3.14159265358979323846  /* pi */

#define PERIOD_X false
#define PERIOD_Y false
#define PERIOD_Z true

extern double L_X;
extern double L_Y;
extern double L_Z;
extern double T;

extern int N;
extern int K;

extern double H_X;
extern double H_Y;
extern double H_Z;

extern double TAU;

extern double C_X;
extern double C_Y;
extern double C_Z;

extern int BX;
extern int BY;
extern int BZ;


inline double u_analytical(double l_x, double l_y, double l_z, double x, double y, double z, double t)
{
	double a_t = PI * sqrt(1/(l_x*l_x) + 1/(l_y*l_y) + 4/(l_z*l_z));
	return sin(PI*x/l_x) * sin(PI*y/l_y) * sin(2*PI*z/l_z) * cos(a_t*t);
}

inline double phi(double l_x, double l_y, double l_z, double x, double y, double z)
{
	// In general:
	// phi == u_analytical(l_x,l_y,l_z,x,y,z,0);
	return sin(PI*x/l_x) * sin(PI*y/l_y) * sin(2*PI*z/l_z);
}
