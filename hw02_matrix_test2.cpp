#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
int main()
{
	
	clock_t t1, t2;				// variables for computing clocks 
	double **A, *x, *b, *c, t, T1, T2;
	int i, j, k, L, N=16;
	int P = 4;
	omp_set_num_threads(P);

	t1 = clock();
	srand(time(NULL));

	for(N=20000;N<=20000;N*=2)
//	for(N=32;N<=32;N*=2)
	{
	
		A = (double **) malloc( N * sizeof(double*) );
		A[0] = (double *) malloc( N*N*sizeof(double));
		x = (double *) malloc( N * sizeof(double) );
		b = (double *) malloc( N * sizeof(double) );
		c = (double *) malloc( N * sizeof(double) );
	
		for(i=1;i<N;++i) A[i] = A[i-1] + N;
		for(i=0;i<N;++i)
		{
			for(j=0;j<N;++j)
			{
				A[i][j] = rand();
			}
			x[i] = rand();
		}

/*		for(i=0;i<N;++i)				//印出矩陣A和向量x 
		{
			for(j=0;j<N;++j)
			{
				printf("%f ",A[i][j]);
			}
			printf("\n");
		}
		printf("-----------------\n");
		for(i=0;i<N;++i)
		{
			printf("%f\n",x[i]);
		}
*/			
		//平行前	
		t1 = clock();
		for(i=0;i<N;++i) 
		{
			c[i] = 0.0;
			for(j=0;j<N;++j)
			{
				c[i] += A[i][j]*x[j];
			}
//			printf("%f\n",c[i]);
		}
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
//		printf("b[N]=%f\n",c[N]);
/*		for(i=0;i<N;++i)
		{
			printf("%f\n",c[i]);
		}
*/		printf("Matrix time vector :%f\n",T1);
		printf("-----------------\n");
		
		
		t1 = clock();
		#pragma omp parallel
		{	
			#pragma omp for reduction(+: t) private(i,j)
			for(i=0;i<N;++i) 
			{
				t = 0.0;
				for(j=0;j<N;++j)
				{
					t += A[i][j]*x[j];
				}
				b[i] = t;
//				printf("%f\n",b[i]);
			}
		}
		t2 = clock();
		T2 = (t2-t1)/(double) CLOCKS_PER_SEC;

//		printf("b[N]=%f\n",b[N]);
/*		for(i=0;i<N;++i)
		{
			printf("%f\n",b[i]);
		}
*/		printf("OMP Matrix time vector :%f\n",T2);
		printf("OMP time: %f\n",T1-T2);
		
		printf("-----------------\n");
/*		//平行後
		t1 = clock();
		#pragma omp parallel for private(i,j,t)
		for(i=0;i<N;++i) 
		{
			k = omp_get_thread_num();
			t = 0.0;
			#pragma omp parallel for 
			for(j=0;j<N;++j)
			{
				L = omp_get_thread_num();
//				printf("%d %d\n",k,L);
				t += A[i][j]*x[j];
			}
			b[i] = t;
			
		}
		t2 = clock();
		T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
//		for(i=0;i<N;++i)
//		{
//			printf("%f\n",fabs(b[i]-c[i]));
//		}
		for(i=0;i<N;++i)
		{
			printf("%f\n",b[i]);
		}
		printf("Matrix time vector :%f\n",T2);
*/		
		
		free(b);
		free(c);
		free(x);
		free(A[0]);
		free(A);
	
	} 

	return 0;
} 
