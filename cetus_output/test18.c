/*
Determinant of Matrix

*/
/* #include<stdio.h> */
int a[1000][1000], m;
int determinant(int f[1000][1000], int x)
{
	int pr, c[1000], d = 0, b[1000][1000], j, p, q, t;
	int _ret_val_0;
	if (x==2)
	{
		d=0;
		d=((f[1][1]*f[2][2])-(f[1][2]*f[2][1]));
		return d;
	}
	else
	{
		#pragma cetus private(j, p, pr, q, r, s, t) 
		#pragma loop name determinant#0 
		for (j=1; j<=x; j ++ )
		{
			int r = 1, s = 1;
			#pragma cetus private(p, q) 
			#pragma loop name determinant#0#0 
			for (p=1; p<=x; p ++ )
			{
				#pragma cetus private(q) 
				#pragma loop name determinant#0#0#0 
				for (q=1; q<=x; q ++ )
				{
					if ((p!=1)&&(q!=j))
					{
						b[r][s]=f[p][q];
						s ++ ;
						if (s>(x-1))
						{
							r ++ ;
							s=1;
						}
					}
				}
			}
			#pragma cetus private(t) 
			#pragma cetus lastprivate(pr) 
			#pragma loop name determinant#0#1 
			for (((t=1), (pr=1)); t<=(1+j); t ++ )
			{
				pr=(( - 1)*pr);
			}
			c[j]=(pr*determinant(b, x-1));
		}
		#pragma cetus private(j) 
		#pragma cetus lastprivate(d) 
		#pragma loop name determinant#1 
		for (((j=1), (d=0)); j<=x; j ++ )
		{
			d=(d+(f[1][j]*c[j]));
		}
		return d;
	}
	return _ret_val_0;
}

int main()
{
	/*
	
		int m, i, j;
	  	printf("\n\nEnter order of matrix : ");
	  	scanf("%d",&m);
	  	printf("\nEnter the elements of matrix\n");
	  	
	  	for(i=1;i<=m;i++) {
				for(j=1;j<=m;j++) {
			  			printf("a[%d][%d] = ",i,j);
			  			scanf("%d",&a[i][j]);
		  		}
	  	}
	  
		for(i=1;i<=m;i++) {
		          printf("\n");
		          for(j=1;j<=m;j++) {    
			               printf("\t%d \t",a[i][j]);
		          }
	     }
	    
	  	printf("\n Determinant of Matrix A is %d .",determinant(a,m));
	
	*/
	int _ret_val_0;
	_ret_val_0=0;
	return _ret_val_0;
}
