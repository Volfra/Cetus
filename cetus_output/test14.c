/*
Jacobi method

*/
/* #include <stdio.h> */
/* #include <math.h> */
void driver(void );
void initialize(void );
void jacobi(void );
void error_check(void );
int n, m, mits;
double tol, relax = 1.0, alpha = 0.0543;
double u[200][200], f[200][200], uold[200][200];
double dx, dy;
int main()
{
	/*  float toler; */
	/*
	     printf("Input n,m (< %d) - grid dimension in x,y direction:\n",MSIZE); 
	      scanf ("%d",&n);
	      scanf ("%d",&m);
	      printf("Input tol - error tolerance for iterative solver\n"); 
	      scanf("%f",&toler);
	      tol=(double)toler;
	      printf("Input mits - Maximum iterations for solver\n"); 
	      scanf("%d",&mits);
	
	*/
	int _ret_val_0;
	n=200;
	m=200;
	tol=1.0E-10;
	mits=1000;
	driver();
	_ret_val_0=0;
	return _ret_val_0;
}

/*

* Subroutine driver () 
* This is where the arrays are allocated and initialzed. 
*
* Working varaibles/arrays 
*     dx  - grid spacing in x direction 
*     dy  - grid spacing in y direction 

*/
void driver()
{
	initialize();
	/* Solve Helmholtz equation */
	jacobi();
	/* error_check (n,m,alpha,dx,dy,u,f) */
	error_check();
	return ;
}

/*
     subroutine initialize (n,m,alpha,dx,dy,u,f) 

* Initializes data 
* Assumes exact solution is u(x,y) = (1-x^2)*(1-y^2)
*

*/
void initialize()
{
	int i, j, xx, yy;
	/*  double PI = 3.1415926; */
	dx=(2.0/(n-1));
	/* -->dx@112:2 */
	dy=(2.0/(m-1));
	/* -->dy@113:2 */
	/* Initialize initial condition and RHS */
	#pragma cetus private(i, j, xx, yy) 
	#pragma loop name initialize#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(3L*n))+((6L*m)*n)))) private(i, j, xx, yy)
	for (i=0; i<n; i ++ )
	{
		#pragma cetus private(j, xx, yy) 
		#pragma loop name initialize#0#0 
		for (j=0; j<m; j ++ )
		{
			xx=((int)(( - 1.0)+(dx*(i-1))));
			/* -1 < x < 1 */
			yy=((int)(( - 1.0)+(dy*(j-1))));
			/* -1 < y < 1 */
			u[i][j]=0.0;
			f[i][j]=(((((( - 1.0)*alpha)*(1.0-(xx*xx)))*(1.0-(yy*yy)))-(2.0*(1.0-(xx*xx))))-(2.0*(1.0-(yy*yy))));
		}
	}
	return ;
}

/*
     subroutine jacobi (n,m,dx,dy,alpha,omega,u,f,tol,maxit)

* Subroutine HelmholtzJ
* Solves poisson equation on rectangular grid assuming : 
* (1) Uniform discretization in each direction, and 
* (2) Dirichlect boundary conditions 
* 
* Jacobi method is used in this routine 
*
* Input : n,m   Number of grid points in the X/Y directions 
*         dx,dy Grid spacing in the X/Y directions 
*         alpha Helmholtz eqn. coefficient 
*         omega Relaxation factor 
*         f(n,m) Right hand side function 
*         u(n,m) Dependent variable/Solution
*         tol    Tolerance for iterative solver 
*         maxit  Maximum number of iterations 
*
* Output : u(n,m) - Solution 

*/
void jacobi()
{
	double omega;
	int i, j, k;
	double error, resid, ax, ay, b;
	omega=relax;
	/* Initialize coefficients */
	ax=(1.0/(dx*dx));
	/* X-direction coef */
	ay=(1.0/(dy*dy));
	/* Y-direction coef */
	b=(((( - 2.0)/(dx*dx))-(2.0/(dy*dy)))-alpha);
	/* Central coeff */
	error=(10.0*tol);
	k=1;
	while ((k<=mits)&&(error>tol))
	{
		error=0.0;
		#pragma cetus private(i, j) 
		#pragma loop name jacobi#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(3L*n))+((3L*m)*n)))) private(i, j)
		for (i=0; i<n; i ++ )
		{
			#pragma cetus private(j) 
			#pragma loop name jacobi#0#0 
			for (j=0; j<m; j ++ )
			{
				uold[i][j]=u[i][j];
			}
		}
		#pragma cetus private(i, j, resid) 
		#pragma loop name jacobi#1 
		#pragma cetus reduction(+: error) 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(((15L+(-10L*m))+(-7L*n))+((5L*m)*n)))) private(i, j, resid) reduction(+: error)
		for (i=1; i<(n-1); i ++ )
		{
			#pragma cetus private(j, resid) 
			#pragma loop name jacobi#1#0 
			/* #pragma cetus reduction(+: error)  */
			for (j=1; j<(m-1); j ++ )
			{
				resid=(((((ax*(uold[i-1][j]+uold[i+1][j]))+(ay*(uold[i][j-1]+uold[i][j+1])))+(b*uold[i][j]))-f[i][j])/b);
				u[i][j]=(uold[i][j]-(omega*resid));
				error=(error+(resid*resid));
			}
		}
		/* Error check */
		/* k = k + 1; */
		error=(sqrt(error)/(n*m));
	}
	/*  End iteration loop */
	/* printf ("Total Number of Iterations:%d\n", k); */
	/* printf ("Residual:%E\n", error); */
	return ;
}

void error_check()
{
	int i, j;
	double xx, yy, temp, error;
	dx=(2.0/(n-1));
	dy=(2.0/(m-1));
	error=0.0;
	#pragma cetus private(i, j, temp, xx, yy) 
	#pragma loop name error_check#0 
	#pragma cetus reduction(+: error) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(3L*n))+((6L*m)*n)))) private(i, j, temp, xx, yy) reduction(+: error)
	for (i=0; i<n; i ++ )
	{
		#pragma cetus private(j, temp, xx, yy) 
		#pragma loop name error_check#0#0 
		/* #pragma cetus reduction(+: error)  */
		for (j=0; j<m; j ++ )
		{
			xx=(( - 1.0)+(dx*(i-1)));
			yy=(( - 1.0)+(dy*(j-1)));
			temp=(u[i][j]-((1.0-(xx*xx))*(1.0-(yy*yy))));
			error=(error+(temp*temp));
		}
	}
	error=(sqrt(error)/(n*m));
	/* printf ("Solution Error :%E \n", error); */
	return ;
}
