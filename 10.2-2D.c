#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define N 3
#define eps 1e-5
//Método de Newton
double f1( double x[N] )
{	return 8* x[1] + 4* pow(x[2],2) -9;
}

double f2( double x[N] )
{	return 8* x[1]* x[2] + 4;
}

double f3( double x[N] )
{	return 10* x[0] - 2* pow(x[1],2) + 2*x[1] - 2* x[2] - 5;
}

double h1( double x[N] )
{	return 0;
}

double h2( double x[N] )
{	return 8;
}

double h3( double x[N] )
{	return 8* x[2];
}

double h4( double x[N] )
{	return 0;
}

double h5( double x[N] )
{	return 8* x[2];
}

double h6( double x[N] )
{	return 8* x[1];
}

double h7( double x[N] )
{	return 10;
}

double h8( double x[N] )
{	return 4* x[1] +1;
}

double h9( double x[N] )
{	return -2;
}

void imprime( double **M, int NL, int NC)
{
	int i, j;
	for( i = 0; i < NL; i++)
	{	for( j = 0; j < NC; j++)
			printf("%lf\t ", M[i][j]);
		puts("");
	}
	puts("");	
}

void **pivoteamento( double **M, int NL, int NC)
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
	double sum;

	for( i = dim-1; i >= 0; i--)
	{	sum = 0;
		for( j = i+1; j < dim; j++)
			sum += M[i][j]* x[j];
		x[i] = ( M[i][dim] - sum )/ M[i][i];
	}
}

void main()
{	double *y, x0[N] = {0.1,0.1,-0.1}, x[N] = {0}, **L;
	double norm, norma;
	double (*F[N])() = {f1,f2,f3};
	double (*J[N][N])() = {{h1,h2,h3}, {h4,h5,h6}, {h7,h8,h9}};
	int i, j;

	y = malloc( N* sizeof(double));
	L = malloc( N* sizeof(double));
	for( i = 0; i < N; i++ )
		L[i] = malloc((N+1)* sizeof(double));

	printf("x[0] \t\tx[1] \t\tx[2] \t\tNorma \n%f \t%f \t%f \t%f\n",x0[0], x0[1], x0[2], sqrt(fabs(norma-norm)));
	do
	{	norma = norm = 0;
		//criação da matriz
		for( i = 0; i < N; i++ )
		{	for( j = 0; j < N; j++ )
				L[i][j] = J[i][j](x0);
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
