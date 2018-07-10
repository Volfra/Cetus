/*
PI calculation

*/
/* #include <stdlib.h> */
/* #include <stdio.h> */
/* #include <time.h> */
int inside_circle(double x, double y)
{
	int rad = 1;
	int circle_x = 1, circle_y = 1;
	int _ret_val_0;
	if ((((x-circle_x)*(x-circle_x))+((y-circle_y)*(y-circle_y)))<=(rad*rad))
	{
		_ret_val_0=1;
		return _ret_val_0;
	}
	else
	{
		_ret_val_0=0;
		return _ret_val_0;
	}
	return _ret_val_0;
}

int main()
{
	int npoints = 10000;
	int circle_count = 0;
	int j;
	double PI;
	int _ret_val_0;
	srand(time(NULL));
	#pragma cetus private(j, xcoordinate, ycoordinate) 
	#pragma loop name main#0 
	/* #pragma cetus reduction(+: circle_count)  */
	for (j=1; j<=npoints; j ++ )
	{
		double xcoordinate = (((double)rand())/((double)RAND_MAX))*2.0;
		double ycoordinate = (((double)rand())/((double)RAND_MAX))*2.0;
		/*
		
		  	printf("%f ",xcoordinate);
		  	printf("- %f",ycoordinate);
		  	printf(" val: %d\n",inside_circle(xcoordinate, ycoordinate));
		  	
		*/
		if (inside_circle(xcoordinate, ycoordinate)==1)
		{
			circle_count=(circle_count+1);
		}
	}
	PI=((4.0*circle_count)/npoints);
	/* printf ("%f\n",PI); */
	_ret_val_0=0;
	return _ret_val_0;
}
