
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define N 2
#define eps 1e-5
#define dt 1e-6
//Método Quasi-Newton
double f1( double x[N] )
{	return 2* x[0] - 5* x[1] +8 +1/2;
}

double f2( double x[N] )
{	return 4* pow(x[0],2) -20* x[0] + pow(x[1],2)/4 +8;
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

void imprime( double **M, int NL, int NC)
{
	int i, j;
	for( i = 0; i < NL; i++)
	{	for( j = 0; j < NC; j++)
			printf("%lf\t ", M[i][j]);
		puts("");
	}	
}

void pivoteamento( double **M, int NL, int NC)
{
        double l, pivo, maior, aux;
        int i, j, k, m, n;
        
        for( j = 0; j < NL-1; j++)
        {	pivo = M[j][j];
		maior = pivo;

		for( k = j; k < NL; k++)
		{	if( fabs(maior) < fabs(M[k][j]))
			{	maior = M[k][j];
				m = k;
			}
		}
		if( maior != pivo)
		{	for( n = 0; n < NC; n++)
			{	aux = M[m][n];
				M[m][n] = M[j][n];
				M[j][n] = aux;
			}
		}
		for( i = j+1; i < NL; i++)
		{	l = M[i][j]/M[j][j];       
			for( k = 0; k < NC; k++ )
				M[i][k] -= l* M[j][k]; 
		}
	}       
}

void subreversa( double **M, double *x, int dim)
{
	int i, j;
	double s;

	for( i = dim-1; i >= 0; i--)
	{	s = 0;
		for( j = i+1; j < dim; j++)
			s =s+ M[i][j]* x[j];
		x[i] = ( M[i][dim] - s )/ M[i][i];
	}
}

void main()
{	double *y, x0[N] = {1,1}, **L;
	double norm, norma;
	double (*F[N])() = {f1,f2};
	double **J;
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

		for( i = 0; i < N; i++ )
		{	for( j = 0; j < N; j++ )
				J[i][j] = df(F[i], x0, j);
		}
		for( i = 0; i < N; i++ )
		{	for( j = 0; j < N; j++ )
				L[i][j] = J[i][j];
			L[i][N] = -F[i](x0);
		}

		pivoteamento(L, N, N+1);
		subreversa(L, y, N);
		for( i = 0; i < N; i++ )
		{	norma += pow(x0[i], 2);
			x0[i] = x0[i] + y[i];
			norm += pow(x0[i], 2);
		}
		printf("%f \t%f \t%f \t%f\n",x0[0], x0[1], x0[2], sqrt(fabs(norma-norm)));	
	}while( eps < sqrt(fabs(norma-norm)));
}
