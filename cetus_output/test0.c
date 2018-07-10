/*
Examples

*/
int main()
{
	int a[10000], b[10000], c[10000], d[100000000][100000000];
	int k;
	int i;
	int j;
	int _ret_val_0;
	#pragma cetus private(k) 
	#pragma loop name main#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(k)
	for (k=0; k<10000; k ++ )
	{
		a[k]=k;
		b[k]=(k-10000);
		c[k]=1;
	}
	/* Flow dependence */
	#pragma cetus private(i) 
	#pragma loop name main#1 
	for (i=1; i<10000; i ++ )
	{
		a[i]=b[i];
		c[i]=a[i-1];
	}
	#pragma cetus private(i) 
	#pragma loop name main#2 
	#pragma cetus parallel 
	#pragma omp parallel for private(i)
	for (i=1; i<10000; i ++ )
	{
		a[i]=b[i];
		c[i]=(a[i]+b[i-1]);
	}
	/* Antidependence */
	#pragma cetus private(i) 
	#pragma loop name main#3 
	for (i=1; i<10000; i ++ )
	{
		a[i-1]=b[i];
		c[i]=a[i];
	}
	/* Output dependence */
	#pragma cetus private(i) 
	#pragma loop name main#4 
	for (i=1; i<10000; i ++ )
	{
		a[i]=b[i];
		a[i+1]=c[i];
	}
	#pragma cetus private(i, j) 
	#pragma loop name main#5 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j)
	for (i=0; i<10000; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name main#5#0 
		for (j=0; j<10000; j ++ )
		{
			d[i][j]=(i+j);
		}
	}
	/* loop interchange */
	#pragma cetus private(i, j) 
	#pragma loop name main#6 
	for (i=0; i<10000; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name main#6#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(j)
		for (j=0; j<10000; j ++ )
		{
			d[i+1][j+2]=(d[i][j]+1);
		}
	}
	#pragma cetus private(i, j) 
	#pragma loop name main#7 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j)
	for (i=0; i<10000; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name main#7#0 
		for (j=0; j<10000; j ++ )
		{
			d[i][j+2]=(d[i][j]+1);
		}
	}
	#pragma cetus private(i, j) 
	#pragma loop name main#8 
	for (i=0; i<10000; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name main#8#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(j)
		for (j=0; j<10000; j ++ )
		{
			d[i+1][j-2]=(d[i][j]+1);
		}
	}
	_ret_val_0=0;
	return _ret_val_0;
}
