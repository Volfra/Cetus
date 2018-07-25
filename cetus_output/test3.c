/*
 Kronecker product

*/
/* #include <stdio.h> */
int main()
{
	int A[100][100];
	int B[100][100];
	int C[10000][10000];
	int i, k, j, l;
	int _ret_val_0;
	{
		int c = 0;
		#pragma cetus private(d) 
		#pragma loop name main#0 
		for (; c<100; c ++ )
		{
			{
				int d = 0;
				#pragma loop name main#0#0 
				for (; d<100; d ++ )
				{
					A[c][d]=1;
					B[c][d]=2;
				}
			}
		}
	}
	#pragma cetus private(i, j, k, l) 
	#pragma loop name main#1 
	for (i=0; i<100; i ++ )
	{
		#pragma cetus private(j, k, l) 
		#pragma loop name main#1#0 
		for (k=0; k<100; k ++ )
		{
			#pragma cetus private(j, l) 
			#pragma loop name main#1#0#0 
			#pragma cetus parallel 
			#pragma omp parallel for private(j, l)
			for (j=0; j<100; j ++ )
			{
				#pragma cetus private(l) 
				#pragma loop name main#1#0#0#0 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(l)
				*/
				for (l=0; l<100; l ++ )
				{
					/* printf ("row: %d",i+l+1); */
					/* printf (" col: %d \n",j+k+l); */
					C[(i+l)+1][(j+k)+1]=(A[i][j]*B[k][l]);
				}
			}
		}
	}
	_ret_val_0=0;
	return _ret_val_0;
}
