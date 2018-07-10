/*
Partition problem DP

*/
/* #include <stdio.h> */
/* Returns true if arr[] can be partitioned in two subsets of */
/* equal sum, otherwise false */
int findPartiion(int arr[], int n)
{
	int sum = 0;
	int i, j;
	/* Calculate sun of all elements */
	int part[((sum/2)+1)][(n+1)];
	int _ret_val_0;
	#pragma cetus private(i) 
	#pragma loop name findPartiion#0 
	#pragma cetus reduction(+: sum) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*n)))) private(i) reduction(+: sum)
	for (i=0; i<n; i ++ )
	{
		sum+=arr[i];
	}
	if ((sum%2)!=0)
	{
		_ret_val_0=0;
		return _ret_val_0;
	}
	/* initialize top row as true */
	#pragma cetus private(i) 
	#pragma loop name findPartiion#1 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(4L+(3L*n)))) private(i)
	for (i=0; i<=n; i ++ )
	{
		part[0][i]=1;
	}
	/* initialize leftmost column, except part[0][0], as 0 */
	#pragma cetus private(i) 
	#pragma loop name findPartiion#2 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*(sum/2L))))) private(i)
	for (i=1; i<=(sum/2); i ++ )
	{
		part[i][0]=0;
	}
	/* Fill the partition table in botton up manner  */
	#pragma cetus private(i, j) 
	#pragma loop name findPartiion#3 
	for (i=1; i<=(sum/2); i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name findPartiion#3#0 
		for (j=1; j<=n; j ++ )
		{
			part[i][j]=part[i][j-1];
			if (i>=arr[j-1])
			{
				part[i][j]=(part[i][j]||part[i-arr[j-1]][j-1]);
			}
		}
	}
	/* uncomment this part to print table  */
	/*
	
		for (i = 0; i <= sum2; i++) {
			   for (j = 0; j <= n; j++) {  
					  printf ("%4d", part[i][j]);
				}
				printf("\n");
		} 
		
	*/
	_ret_val_0=part[sum/2][n];
	return _ret_val_0;
}

int main()
{
int arr[] = {3, 1, 1, 2, 2, 1};
	int r;
	int n = sizeof arr/sizeof arr[0];
	int _ret_val_0;
	if (findPartiion(arr, n)==1)
	{
		r=1;
		/* printf("Can be divided into two subsets of equal sum"); */
	}
	else
	{
		r=0;
		/* printf("Can not be divided into two subsets of equal sum"); */
	}
	_ret_val_0=0;
	return _ret_val_0;
}
