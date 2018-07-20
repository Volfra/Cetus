/*
Shear sort algorithm

*/
/* #include <math.h> */
/* #include <stdio.h> */
int n;
int a[1005][1005];
/*

void rowsort(int i){
	    for(int j=0;j<n-1;j++){
		        for(int k=0;k<n-j-1;k++){
			            if(a[i][k]>a[i][k+1]){
				                int temp=a[i][k];
				                a[i][k]=a[i][k+1];
				                a[i][k+1]=temp;
			            }
		        }
	    }
}

void rowrevsort(int i){
	    for(int j=0;j<n-1;j++){
		        for(int k=0;k<n-j-1;k++){
			            if(a[i][k]<a[i][k+1]){
				                int temp=a[i][k];
				                a[i][k]=a[i][k+1];
				                a[i][k+1]=temp;
			            }
		        }
	    }
}

void colsort(int i){
	    for(int j=0;j<n-1;j++){
		        for(int k=0;k<n-j-1;k++){
			            if(a[k][i]>a[k+1][i]){
				                int temp=a[k][i];
				                a[k][i]=a[k+1][i];
				                a[k+1][i]=temp;
			            }
		        }
	    }
}

*/
int main()
{
	int m = (int)ceil(log2(n));
	int _ret_val_0;
	n=1005;
	{
		int j = 0;
		#pragma cetus private(k) 
		#pragma loop name main#0 
		for (; j<n; j ++ )
		{
			{
				int k = 0;
				#pragma loop name main#0#0 
				for (; k<n; k ++ )
				{
					a[j][k]=(n-j);
				}
			}
		}
	}
	{
		int i = 0;
		#pragma cetus private(j, k, t, temp) 
		#pragma loop name main#1 
		for (; i<m; i ++ )
		{
			{
				int j = 0;
				#pragma cetus lastprivate(k, t, temp) 
				#pragma loop name main#1#0 
				for (; j<n; j ++ )
				{
					if ((j%2)==0)
					{
						/* rowsort(j); */
						{
							int t = 0;
							#pragma cetus lastprivate(k, temp) 
							#pragma loop name main#1#0#0 
							for (; t<(n-1); t ++ )
							{
								{
									int k = 0;
									#pragma cetus lastprivate(temp) 
									#pragma loop name main#1#0#0#0 
									for (; k<((n-t)-1); k ++ )
									{
										if (a[j][k]>a[j][k+1])
										{
											int temp = a[j][k];
											a[j][k]=a[j][k+1];
											a[j][k+1]=temp;
										}
									}
								}
							}
						}
					}
					else
					{
						/* rowrevsort(j); */
						{
							int t = 0;
							#pragma cetus lastprivate(k, temp) 
							#pragma loop name main#1#0#1 
							for (; t<(n-1); t ++ )
							{
								{
									int k = 0;
									#pragma cetus lastprivate(temp) 
									#pragma loop name main#1#0#1#0 
									for (; k<((n-t)-1); k ++ )
									{
										if (a[j][k]<a[j][k+1])
										{
											int temp = a[j][k];
											a[j][k]=a[j][k+1];
											a[j][k+1]=temp;
										}
									}
								}
							}
						}
					}
				}
			}
			{
				int j = 0;
				#pragma cetus lastprivate(k, t, temp) 
				#pragma loop name main#1#1 
				for (; j<n; j ++ )
				{
					/* colsort(j); */
					{
						int t = 0;
						#pragma cetus lastprivate(k, temp) 
						#pragma loop name main#1#1#0 
						for (; t<(n-1); t ++ )
						{
							{
								int k = 0;
								#pragma cetus lastprivate(temp) 
								#pragma loop name main#1#1#0#0 
								for (; k<((n-t)-1); k ++ )
								{
									if (a[k][j]>a[k+1][j])
									{
										int temp = a[k][j];
										a[k][j]=a[k+1][j];
										a[k+1][j]=temp;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	{
		int j = 0;
		#pragma cetus private(k, t, temp) 
		#pragma loop name main#2 
		for (; j<n; j ++ )
		{
			if ((j%2)==0)
			{
				/* rowsort(j); */
				{
					int t = 0;
					#pragma cetus lastprivate(k, temp) 
					#pragma loop name main#2#0 
					for (; t<(n-1); t ++ )
					{
						{
							int k = 0;
							#pragma cetus lastprivate(temp) 
							#pragma loop name main#2#0#0 
							for (; k<((n-t)-1); k ++ )
							{
								if (a[j][k]>a[j][k+1])
								{
									int temp = a[j][k];
									a[j][k]=a[j][k+1];
									a[j][k+1]=temp;
								}
							}
						}
					}
				}
			}
			else
			{
				/* rowrevsort(j); */
				{
					int t = 0;
					#pragma cetus lastprivate(k, temp) 
					#pragma loop name main#2#1 
					for (; t<(n-1); t ++ )
					{
						{
							int k = 0;
							#pragma cetus lastprivate(temp) 
							#pragma loop name main#2#1#0 
							for (; k<((n-t)-1); k ++ )
							{
								if (a[j][k]<a[j][k+1])
								{
									int temp = a[j][k];
									a[j][k]=a[j][k+1];
									a[j][k+1]=temp;
								}
							}
						}
					}
				}
			}
		}
	}
	/*
	
		printf("\n Matrix of input data:\n");
		int i,j;
		for(i=0;i<n;i++)  {
		  		for(j=0;j<n;j++) 
		   			printf("%d \t",a[i][j]);
		  		printf("\n");
	 	}
	   
	*/
	_ret_val_0=0;
	return _ret_val_0;
}
