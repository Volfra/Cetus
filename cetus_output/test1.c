/*
 Matrix multiplication

*/
int main()
{
	int first[1000][1000], second[1000][1000], multiply[1000][1000];
	int c, d, k, sum = 0;
	int m, p, q;
	int _ret_val_0;
	m=(p=(q=1000));
	#pragma cetus private(c, d, k) 
	#pragma loop name main#0 
	for (c=0; c<m; c ++ )
	{
		d=0;
		#pragma cetus private(k) 
		#pragma loop name main#0#0 
		#pragma cetus reduction(+: sum) 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(k) reduction(+: sum)
		*/
		for (k=0; k<p; k ++ )
		{
			sum=(sum+(first[c][k]*second[k][d]));
		}
		multiply[c][d]=sum;
		#pragma cetus private(d, k, sum) 
		#pragma loop name main#0#1 
		#pragma cetus parallel 
		#pragma omp parallel for private(d, k, sum)
		for (d=1; d<q; d ++ )
		{
			sum=0;
			#pragma cetus private(k) 
			#pragma loop name main#0#1#0 
			#pragma cetus reduction(+: sum) 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(k) reduction(+: sum)
			*/
			for (k=0; k<p; k ++ )
			{
				sum=(sum+(first[c][k]*second[k][d]));
			}
			multiply[c][d]=sum;
		}
		sum=0;
	}
	_ret_val_0=0;
	return _ret_val_0;
}
