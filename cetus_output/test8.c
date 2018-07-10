/*
Diagonally dominant matrix

*/
int main()
{
	int A[10005][10005];
	int i, j;
	int boolean = 1;
	int sum = 0;
	int _ret_val_0;
	#pragma cetus private(i, j) 
	#pragma loop name main#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j)
	for (i=0; i<10005; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name main#0#0 
		for (j=0; j<10005; j ++ )
		{
			A[i][j]=(i+j);
		}
	}
	#pragma cetus private(j) 
	#pragma loop name main#1 
	for (i=0; i<10005; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name main#1#0 
		#pragma cetus reduction(+: sum) 
		#pragma cetus parallel 
		#pragma omp parallel for private(j) reduction(+: sum)
		for (j=0; j<10005; j ++ )
		{
			if (i!=j)
			{
				sum+=A[i][j];
			}
		}
		if (A[i][i]<sum)
		{
			boolean=0;
			break;
		}
		sum=0;
	}
	_ret_val_0=0;
	return _ret_val_0;
}
