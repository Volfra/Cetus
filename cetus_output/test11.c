/*
LU decomposition

*/
/* #include<stdio.h> */
/* #include<conio.h> */
int main()
{
float A[1000][1000] = {0}, L[1000][1000] = {0}, U[1000][1000];
float B[1000] = {0}, X[1000] = {0}, Y[1000] = {0};
	int i, j, k, n;
	/*
	
	    printf("Enter the order of square matrix: ");
	    scanf("%d",&n);
	    printf("\nEnter matrix element:\n");
	    for(i=0; i<n; i++)
	    {
		        for(j=0; j<n; j++)
		        {
			            printf("Enter A[%d][%d] element: ", i,j);
			            scanf("%f",&A[i][j]);
		        }
	    }
	    
	    printf("\nEnter the constant terms: \n");
	    for(i=0; i<n; i++)
	    {
		        printf("B[%d]",i);
		        scanf("%f",&B[i]);
	    }
	   
	*/
	int _ret_val_0;
	#pragma cetus private(i, j, k) 
	#pragma loop name main#0 
	for (j=0; j<n; j ++ )
	{
		#pragma cetus private(i, k) 
		#pragma loop name main#0#0 
		for (i=0; i<n; i ++ )
		{
			if (i<=j)
			{
				U[i][j]=A[i][j];
				#pragma cetus private(k) 
				#pragma loop name main#0#0#0 
				for (k=0; k<(i-1); k ++ )
				{
					U[i][j]-=(L[i][k]*U[k][j]);
				}
				if (i==j)
				{
					L[i][j]=1;
				}
				else
				{
					L[i][j]=0;
				}
			}
			else
			{
				L[i][j]=A[i][j];
				#pragma cetus private(k) 
				#pragma loop name main#0#0#1 
				for (k=0; k<=(j-1); k ++ )
				{
					L[i][j]-=(L[i][k]*U[k][j]);
				}
				L[i][j]/=U[j][j];
				U[i][j]=0;
			}
		}
	}
	/*
	
	    printf("[L]: \n");
	    for(i=0; i<n; i++)
	    {
		        for(j=0; j<n; j++)
		            printf("%9.3f",L[i][j]);
		        printf("\n");
	    }
	    
	    printf("\n\n[U]: \n");
	    for(i=0; i<n; i++) {
		        for(j=0; j<n; j++)
		            printf("%9.3f",U[i][j]);
		        printf("\n");
	    }
	   
	*/
	#pragma cetus private(i, j) 
	#pragma loop name main#1 
	for (i=0; i<n; i ++ )
	{
		Y[i]=B[i];
		#pragma cetus private(j) 
		#pragma loop name main#1#0 
		for (j=0; j<i; j ++ )
		{
			Y[i]-=(L[i][j]*Y[j]);
		}
	}
	/*
	
	    printf("\n\n[Y]: \n");
	    for(i=0; i<n; i++) {
		        printf("%9.3f",Y[i]);
	    }
	   
	*/
	#pragma cetus private(i, j) 
	#pragma loop name main#2 
	for (i=(n-1); i>=0; i -- )
	{
		X[i]=Y[i];
		#pragma cetus private(j) 
		#pragma loop name main#2#0 
		for (j=(i+1); j<n; j ++ )
		{
			X[i]-=(U[i][j]*X[j]);
		}
		X[i]/=U[i][i];
	}
	/*
	
	    printf("\n\n[X]: \n");
	    for(i=0; i<n; i++) {
		        printf("%9.3f",X[i]);
	    }
	   
	*/
	_ret_val_0=0;
	return _ret_val_0;
}
