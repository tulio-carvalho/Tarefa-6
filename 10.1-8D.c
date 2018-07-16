#include<stdio.h>
#include<math.h>

#define N 3
#define eps 1e-5
// Método do Ponto Fixo
float f1( float x[N] )
{	return sqrt(x[1] - 2.* x[1] + 2.* x[2]);
}

float f2( float x[N] )
{	return sqrt( (10.* x[2]+ pow(x[0],2))/8.);
}
float f3( float x[N] )
{	return (pow(x[0],2)/7.* x[1]);
}

void main()
{	float xa[N] = {0.1,0.1,0.1};
	float norm, norma;
	float (*equacao[N])() = {f1,f2,f3};
	int i;

	printf("x[0] \t\tx[1] \t\tx[2] \t\tNorma\n%f \t%f \t%f\n",xa[0], xa[1], xa[2]);
	do
	{	norma = norm = 0;
		for( i = 0; i < N; i++ )
		{	norm += pow( xa[i], 2);
			xa[i] = equacao[i](xa);
			norma += pow( xa[i], 2);
		}
		printf("%f \t%f \t%f \t%f\n",xa[0], xa[1], xa[2], sqrt(fabs(norma-norm)));
	}while( eps < sqrt(fabs(norma-norm)));
}
