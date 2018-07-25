/*
Transpose matrix

*/
int main()
{
	int A[1005][1005], T[1005][1005];
	int i, j;
	int _ret_val_0;
	#pragma cetus private(i, j) 
	#pragma loop name main#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j)
	for (i=0; i<1005; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name main#0#0 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(j)
		*/
		for (j=0; j<1005; j ++ )
		{
			A[i][j]=(i+j);
		}
	}
	#pragma cetus private(i, j) 
	#pragma loop name main#1 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j)
	for (i=0; i<1005; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name main#1#0 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(j)
		*/
		for (j=0; j<1005; j ++ )
		{
			T[i][j]=A[j][i];
		}
	}
	_ret_val_0=0;
	return _ret_val_0;
}
