/*
 LCS dynamic programming

*/
int main()
{
	int M[10000][10000];
	char X[10000];
	char Y[10000];
	int _ret_val_0;
	{
		int i = 0;
		#pragma cetus private(j) 
		#pragma loop name main#0 
		for (; i<10000; i ++ )
		{
			int j = 0;
			#pragma loop name main#0#0 
			for (; j<10000; j ++ )
			{
				if ((i==0)||(j==0))
				{
					M[i][j]=0;
				}
				else
				{
					if (X[i]==Y[j])
					{
						M[i][j]=(M[i-1][j-1]+1);
					}
					else
					{
						M[i][j]=((M[i-1][j]>M[i][j-1]) ? M[i-1][j] : M[i][j-1]);
					}
				}
			}
		}
	}
	_ret_val_0=0;
	return _ret_val_0;
}
