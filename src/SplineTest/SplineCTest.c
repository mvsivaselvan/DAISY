#include<stdio.h>
#include<stdlib.h>

extern void spcolC(double* knots, int nknots, int k, double* colpts, int ncolpts, int nderiv, double* colmat, int ncolmat);

void main()
{
	double t[] = {0.,0.,0.,0.,0.1429,0.2857,0.4286,0.5714,0.7143,0.8571,1.0,1.0,1.0,1.0,1.0};
	int nt = 14;
	int k = 4;
	double xg[] = {0.0161, 0.0714, 0.1268, 0.1590, 0.2143, 0.2696, 0.3018, 0.3571, 0.4125,
		           0.4447, 0.5000, 0.5553, 0.5875, 0.6429, 0.6982, 0.7304, 0.7857, 0.8410,
		           0.8732, 0.9286, 0.9839};
	int nxg = 21;
	int nderiv = 3;
	int nbasis = nt - k;
	double* colmat = calloc(nxg*nderiv*nbasis,sizeof(double));

	spcolC(t, nt, k, xg, nxg, nderiv, colmat, nxg*nderiv*nbasis);

	printf("%f\n", colmat[2]);

	free(colmat);

	printf("Hello world\n");

}