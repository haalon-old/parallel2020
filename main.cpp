#include <stdio.h>
#include "analytic.h"

#define L_X 1.0
#define L_Y 1.0
#define L_Z 1.0
#define T 1.0

#define N 10
#define K 2

#define H_X L_X / N
#define H_Y L_Y / N
#define H_Z L_Z / N

#define TAU T / K


int main(int argc, char const *argv[])
{
	for(int t = 0; t < K; t++)
	{
		printf("STEP %3d\n", t);
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				for (int k = 0; k < N; ++k)
				{
					printf("%8.3f", u_analytical(L_X, L_Y, L_Z, H_X*i, H_Y*j, H_Z*k, TAU*t));
				}
				printf("    (%4.2f, %4.2f)\n", H_X*i, H_Y*j);
			}
			printf("\n\n");
		}
	}

	return 0;
}