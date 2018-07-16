#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define N 3
#define eps 5e-2
#define dt 1e-6
//Método Steepest Decent.
double f1( double x[N] )
{	
	return 5* x[0] + pow(x[1], 2) - 4* x[2] - 13;
}

double f2( double x[N] )
{	
	return pow(x[0], 2) +10* x[1] - x[2] -11;
}

double f3( double x[N] )
{	
	return pow(x[1], 3) - 25* x[2] + 22;
}

double G( double (*f[N])(), double x[N] )
{	int i;
	double s=0;
	for( i = 0; i < N; i++)
		s=s+pow(f[i](x), 2);
	return s;
}

double* gradiente( double **J, double (*F[N])(), double x[N], double *normaGrad)
{	int i, j;
	double *grad;

	grad = malloc( N* sizeof(double));
	for( i = 0 ; i < N ; i++ )
	{	for( j = 0 ; j < N ; j++ )
			grad[i] += 2* J[j][i]* F[j](x);
		*normaGrad += pow( grad[i], 2 );
	}
	*normaGrad = sqrt( *normaGrad );
	return grad;
}

double df( double f(), double *x, int i)
{	double diff, x1[N], x2[N];
	int j;

	for( j = 0; j < N; j++ )
	{	x1[j] = x[j];
		x2[j] = x[j];
	}
	x1[i] = x[i] + dt;
	x2[i] = x[i] - dt;
	diff = (f(x1) - f(x2))/ (2* dt);
	return diff;
}

double h( double x[N], double *grad, double a, double (*F[N])() )
{	int i;
	double aux[N] = {0};
	for( i = 0 ; i < N ; i++ )
		aux[i] = x[i] -a* grad[i];
	return G(F, aux);
}

void imprime( double **M, int NL, int NC)
{
	int i, j;
	for( i = 0; i < NL; i++)
	{	for( j = 0; j < NC; j++)
			printf("%lf\t ", M[i][j]);
		puts("");
	}	
}

void main()
{	double *y, x0[N] = {0,0,0}, **L, a[N], gn[N], k[N], alpha;
	double norm, norma;
	double (*F[N])() = {f1,f2,f3};
	double **J, g, *gradg, normaGradg;
	int i, j;

	y = malloc( N* sizeof(double));
	L = malloc( N* sizeof(double));
	for( i = 0; i < N; i++ )
		L[i] = malloc((N+1)* sizeof(double));

	J = malloc( N* sizeof(double));
	for( i = 0; i < N; i++ )
		J[i] = malloc(N* sizeof(double));

	printf("x[0] \t\tx[1] \t\tx[2] \t\t Norma\n%f \t%f \t%f\n",x0[0], x0[1], x0[2]);
	do
	{	norma = norm = 0;
		normaGradg = 0;
		for( i = 0; i < N; i++ )
			gn[i] = 0;
		//criação da matriz
		for( i = 0; i < N; i++ )
		{	for( j = 0; j < N; j++ )
				J[i][j] = df(F[i], x0, j);
		}
		for( i = 0; i < N; i++ )
		{	for( j = 0; j < N; j++ )
				L[i][j] = J[i][j];
			L[i][N] = -F[i](x0);
		}
		g = G(F, x0);
		gradg = gradiente(J, F, x0, &normaGradg );
		for( i = 0; i < N; i++ )
			gradg[i] /= normaGradg;

		a[0] = 0;
		a[2] = 1;
		
		while( h(x0, gradg, a[0], F) < h(x0, gradg, a[2], F) )
		{	a[2] /= 2.;
			if( a[2] < 1e-6)
				break;
		}
		if( h(x0, gradg, a[0], F) < h(x0, gradg, a[2], F) )
		{	a[1] = a[2];
			a[2] = a[0];
			a[0] = a[1];
			a[1] = a[0]/2;
		}
		else			
			a[1] = a[2]/ 2.;
	
		for( i = 0; i < N; i++ )
			gn[i] = h(x0, gradg, a[i], F);
		
		for( i = 0; i < N-1; i++ )  
			k[i] = (gn[i+1] - gn[i])/ (a[i+1] - a[i]);
		k[2] = (k[1] - k[0])/ (a[1] - a[0]);

		alpha = (k[2]* a[1] - k[0])/ (2* k[2]);

		for( i = 0; i < N; i++ )
		{	norma += pow(x0[i], 2);
			x0[i] -= alpha* gradg[i];
			norm += pow(x0[i], 2);
		}
	
		printf("%f \t%f \t%f \t%f\n",x0[0], x0[1], x0[2], sqrt(fabs(norma-norm)));	
	}while( eps < sqrt(fabs(norma-norm)));
}
