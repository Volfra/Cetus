/*
 Sieve of Eratosthenes

*/
/* #include <string.h> */
int main()
{
	int prime[1001];
	/* memset(prime, 1, sizeof(prime)); */
	int _ret_val_0;
	{
		int p = 2;
		#pragma cetus private(i) 
		#pragma loop name main#0 
		for (; (p*p)<=1000; p ++ )
		{
			if (prime[p]==1)
			{
				{
					int i = p*2;
					#pragma loop name main#0#0 
					for (; i<=1000; i=(i+p))
					{
						prime[i]=0;
					}
				}
			}
		}
	}
	_ret_val_0=0;
	return _ret_val_0;
}
