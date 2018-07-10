/*
Convex Hull Algorithm

*/
/* #include<stdio.h> */
/* #include<string.h> */
/* #include<math.h> */
struct node
{
	double x;
	double y;
};

typedef struct node node;
double orientation(node p, node q, node r)
{
	double val = ((q.y-p.y)*(r.x-q.x))-((q.x-p.x)*(r.y-q.y));
	double _ret_val_0;
	if (val==0)
	{
		_ret_val_0=0;
		return _ret_val_0;
	}
	if (val>0)
	{
		_ret_val_0=1;
		return _ret_val_0;
	}
	else
	{
		_ret_val_0=2;
		return _ret_val_0;
	}
	return _ret_val_0;
}

int main()
{
	long long int n;
	node points[10000];
	double x, y;
	/* scanf("%lld",&n); */
	int i;
	int p = 0;
	int l;
	int q = (p+1)%n;
	node result[10001];
	int count = 0;
	int _ret_val_0;
	#pragma cetus private(i) 
	#pragma loop name main#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(4L*n)))) private(i)
	for (i=0; i<n; i ++ )
	{
		/* scanf("%lf",&x); */
		/* scanf("%lf",&y); */
		points[i].x=x;
		points[i].y=y;
	}
	#pragma cetus private(i) 
	#pragma loop name main#1 
	for (i=1; i<n; i ++ )
	{
		if (points[p].x>points[i].x)
		{
			p=i;
		}
		else
		{
			if ((points[p].x==points[i].x)&&(points[p].y>points[i].y))
			{
				p=i;
			}
		}
	}
	l=p;
	if (n<3)
	{
		#pragma cetus private(i) 
		#pragma loop name main#2 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(-2L+(3L*n)))) private(i)
		for (i=0; i<(n-1); i ++ )
		{
			/* printf("%lf %lf\n",points[i].x,points[i].y); */
			result[count+i]=points[i];
		}
		if ((-2+n)>=0)
		{
			count+=(-1+n);
		}
	}
	else
	{
		do
		{
			/* printf("%lf %lf\n",points[p].x,points[p].y); */
			result[count]=points[p];
			count ++ ;
			#pragma cetus private(i) 
			#pragma loop name main#3 
			for (i=0; i<n; i ++ )
			{
				if (orientation(points[p], points[i], points[q])==2)
				{
					q=i;
				}
			}
			p=q;
			q=((p+1)%n);
		}while(p!=l);
		
	}
	_ret_val_0=0;
	return _ret_val_0;
}
