int main()
{
	int _ret_val_0;
	_ret_val_0=0;
	return _ret_val_0;
}

int foo3()
{
	double a[10000][10000];
	int i, j;
	int _ret_val_0;
	#pragma cetus private(i, j) 
	#pragma loop name foo3#0 
	for (i=0; i<=9998; i+=1)
	{
		#pragma cetus private(j) 
		#pragma loop name foo3#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(j)
		for (j=0; j<=9999; j+=1)
		{
			a[i][j]+=a[i+1][j];
		}
	}
	_ret_val_0=0;
	return _ret_val_0;
}

void foo4(int x, int y)
{
	double a[10000][10000];
	int i, j;
	#pragma cetus private(i, j) 
	#pragma loop name foo4#0 
	for (i=0; i<=9998; i+=1)
	{
		#pragma cetus private(j) 
		#pragma loop name foo4#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(j)
		for (j=0; j<=9999; j+=1)
		{
			a[i][j]+=a[i+1][j];
		}
	}
	return ;
}
