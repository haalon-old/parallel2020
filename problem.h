#pragma once
#include <cmath>

#define PI 3.14159265358979323846  /* pi */

#define PERIOD_X false
#define PERIOD_Y false
#define PERIOD_Z true

#define L_X 1.0
#define L_Y 1.0
#define L_Z 1.0
#define T 0.2

#define N 32
#define K 32

const double H_X = L_X / N;
const double H_Y = L_Y / N;
const double H_Z = L_Z / N;

const double TAU = T / K;

const double C_X = TAU / H_X;
const double C_Y = TAU / H_Y;
const double C_Z = TAU / H_Z;

#define BX 1
#define BY 1
#define BZ 1


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
