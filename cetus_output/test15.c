/*
Inner product

*/
int dot_product(int * a, int * b, int n)
{
	int i, sum = 0;
	#pragma cetus private(i) 
	#pragma loop name dot_product#0 
	#pragma cetus reduction(+: sum) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*n)))) private(i) reduction(+: sum)
	for (i=0; i<n; i ++ )
	{
		sum+=(a[i]*b[i]);
	}
	return sum;
}

int main()
{
	int a[10000], b[10000];
	int i;
	int _ret_val_0;
	#pragma cetus private(i) 
	#pragma loop name main#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i)
	for (i=0; i<10000; i ++ )
	{
		a[i]=( - 1);
		b[i]=1;
	}
	dot_product(a, b, sizeof a/sizeof a[0]);
	_ret_val_0=0;
	return _ret_val_0;
}
