/*
test Rose
 https:github.com/rose-compiler/rose/tree/master/projects/autoParallelization/tests

*/
typedef double real8;
extern void OtherFunc(int k, real8 * l, real8 * m, real8 * n, real8 * o, real8 * p, real8 q, real8 r, real8 s[3]);
int main()
{
	int _ret_val_0;
	_ret_val_0=0;
	return _ret_val_0;
}

/*

 DelVolBaseLoopAlgorithm

*/
void DelVolBaseLoopAlgorithm()
{
	int j, k, off;
	#pragma cetus private(j, k, off) 
	#pragma loop name DelVolBaseLoopAlgorithm#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(j, k, off)
	for (k=1; k<10000; k ++ )
	{
		#pragma cetus private(j, off) 
		#pragma loop name DelVolBaseLoopAlgorithm#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(j, off)
		for (j=1; j<10000; j ++ )
		{
			off=((j*555)+(k*777));
		}
	}
	return ;
}

/*

 Rose stress

*/
void StressCheckEpsFail(real8 * newSxx, real8 * newSyy, real8 * newSzz, real8 * newTxy, real8 * newTxz, real8 * newTyz, real8 * eps, real8 eps_failure_model, const int * zoneset, int length)
{
	int i;
	int index;
	#pragma cetus private(i, index) 
	#pragma loop name StressCheckEpsFail#0 
	for (i=0; i<=(length-1); i+=1)
	{
		index=zoneset[i];
		if (eps[zoneset[i]]>eps_failure_model)
		{
			newSxx[i]=0.0;
			newSyy[i]=0.0;
			newSzz[i]=0.0;
			newTxy[i]=0.0;
			newTxz[i]=0.0;
			newTyz[i]=0.0;
			eps[zoneset[i]]=(eps_failure_model*1.01);
		}
	}
	return ;
}

void StressStrainWork(real8 * deltz, real8 * delts, const real8 * newSxx, const real8 * newSyy, const real8 * newSzz, const real8 * newTxy, const real8 * newTxz, const real8 * newTyz, const real8 * sxx, const real8 * syy, const real8 * txy, const real8 * txz, const real8 * tyz, const real8 * dxx, const real8 * dyy, const real8 * dzz, const real8 * dxy, const real8 * dxz, const real8 * dyz, real8 deltaTime, const int * zoneset, const real8 * vc, const real8 * vnewc, int length)
{
	int i;
	int index;
	real8 quarterDelta = 0.25*deltaTime;
	real8 szz;
	#pragma cetus private(i, index) 
	#pragma loop name StressStrainWork#0 
	#pragma cetus reduction(+: deltz[zoneset[i]]) 
	for (i=0; i<=(length-1); i+=1)
	{
		index=zoneset[i];
		szz=(( - sxx[zoneset[i]])-syy[zoneset[i]]);
		deltz[zoneset[i]]+=((quarterDelta*(vnewc[i]+vc[i]))*((((((dxx[zoneset[i]]*(sxx[zoneset[i]]+newSxx[i]))+(dyy[zoneset[i]]*(syy[zoneset[i]]+newSyy[i])))+(dzz[zoneset[i]]*(szz+newSzz[i])))+((2.0*dxy[zoneset[i]])*(txy[zoneset[i]]+newTxy[i])))+((2.0*dxz[zoneset[i]])*(txz[zoneset[i]]+newTxz[i])))+((2.0*dyz[zoneset[i]])*(tyz[zoneset[i]]+newTyz[i]))));
		delts[i]+=((quarterDelta*(vnewc[i]+vc[i]))*((((((dxx[zoneset[i]]*sxx[zoneset[i]])+(dyy[zoneset[i]]*syy[zoneset[i]]))+(dzz[zoneset[i]]*szz))+((2.0*dxy[zoneset[i]])*txy[zoneset[i]]))+((2.0*dxz[zoneset[i]])*txz[zoneset[i]]))+((2.0*dyz[zoneset[i]])*tyz[zoneset[i]])));
	}
	return ;
}

void StressStrainHeat(const real8 * deltz, real8 * deltzh, real8 * deltrh, const real8 * shearMod, const real8 * shearRatio, const real8 * shearDer, const real8 * newSxx, const real8 * newSyy, const real8 * newSzz, const real8 * newTxy, const real8 * newTxz, const real8 * newTyz, const real8 * sxx, const real8 * syy, const real8 * txy, const real8 * txz, const real8 * tyz, real8 deltaTime, const int * zoneset, const real8 * vc, const real8 * vnewc, int length)
{
	real8 shearr;
	real8 sheari;
	real8 avgMod;
	int nz;
	int i;
	/* Quiet the compiler - unused argument */
	deltaTime=deltaTime;
	#pragma cetus private(i, nz) 
	#pragma loop name StressStrainHeat#0 
	for (i=0; i<=(length-1); i+=1)
	{
		nz=zoneset[i];
		shearr=(0.5*shearRatio[i]);
		if (shearMod[zoneset[i]]>0.0)
		{
			sheari=(0.5/shearMod[zoneset[i]]);
			deltrh[zoneset[i]]=((0.25*(vnewc[i]+vc[i]))*((((((((newSxx[i]*sheari)-(sxx[zoneset[i]]*shearr))*(sxx[zoneset[i]]+newSxx[i]))+(((newSyy[i]*sheari)-(syy[zoneset[i]]*shearr))*(syy[zoneset[i]]+newSyy[i])))+(((newSzz[i]*sheari)+((syy[zoneset[i]]+sxx[zoneset[i]])*shearr))*((newSzz[i]-sxx[zoneset[i]])-syy[zoneset[i]])))+((2.0*((newTxy[i]*sheari)-(txy[zoneset[i]]*shearr)))*(txy[zoneset[i]]+newTxy[i])))+((2.0*((newTxz[i]*sheari)-(txz[zoneset[i]]*shearr)))*(txz[zoneset[i]]+newTxz[i])))+((2.0*((newTyz[i]*sheari)-(tyz[zoneset[i]]*shearr)))*(tyz[zoneset[i]]+newTyz[i]))));
		}
		else
		{
			deltrh[zoneset[i]]=(((( - 0.25)*(vnewc[i]+vc[i]))*((((((sxx[zoneset[i]]*(sxx[zoneset[i]]+newSxx[i]))+(syy[zoneset[i]]*(syy[zoneset[i]]+newSyy[i])))-((syy[zoneset[i]]+sxx[zoneset[i]])*((newSzz[i]-sxx[zoneset[i]])-syy[zoneset[i]])))+((2.0*txy[zoneset[i]])*(txy[zoneset[i]]+newTxy[i])))+((2.0*txz[zoneset[i]])*(txz[zoneset[i]]+newTxz[i])))+((2.0*tyz[zoneset[i]])*(tyz[zoneset[i]]+newTyz[i]))))*shearr);
		}
		deltzh[zoneset[i]]=(deltz[zoneset[i]]-deltrh[zoneset[i]]);
		avgMod=(0.5*shearMod[zoneset[i]]);
		if (shearRatio[i]>0.0)
		{
			avgMod=(avgMod+(0.5/shearRatio[i]));
		}
		if (avgMod>0.0)
		{
			deltrh[zoneset[i]]=((shearDer[i]*deltrh[zoneset[i]])/avgMod);
		}
		else
		{
			deltrh[zoneset[i]]=0.0;
		}
	}
	return ;
}

/*

 Accumulate force

*/
void AccumulateForce(int * idxBound, int * idxList, int len, double * tmp, double * force)
{
	register int ii;
	#pragma cetus private() 
	#pragma loop name AccumulateForce#0 
	#pragma cetus parallel 
	#pragma omp parallel for
	for (ii=0; ii<=(len-1); ii+=1)
	{
		int count = idxBound[ii+1]-idxBound[ii];
		int * list =  & idxList[idxBound[ii]];
		double sum = 0.0;
		register int jj;
		#pragma cetus private() 
		#pragma loop name AccumulateForce#0#0 
		#pragma cetus reduction(+: sum) 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(4L*count)))) reduction(+: sum)
		for (jj=0; jj<=(count-1); jj+=1)
		{
			int idx = list[jj];
			sum+=tmp[list[jj]];
		}
		force[ii]+=sum;
	}
	return ;
}

/*

 An anti-dependence example

*/
void foo1()
{
	int i;
	int a[1000];
	#pragma cetus private(i) 
	#pragma loop name foo1#0 
	for (i=0; i<=998; i+=1)
	{
		a[i]=(a[i+1]+1);
	}
	return ;
}

/*

 Array Scalar

*/
void foo2()
{
	int i;
	int a[1000];
	#pragma cetus private(i) 
	#pragma loop name foo2#0 
	for (i=0; i<=999; i+=1)
	{
		a[i]=(a[i]+a[0]);
	}
	return ;
}

/*

 C99

*/
int foo3()
{
	double a[10000][10000];
	int i, j;
	int _ret_val_0;
	#pragma cetus private(i, j) 
	#pragma loop name foo3#0 
	for (i=0; i<=9998; i+=1)
	{
		#pragma cetus private(j) 
		#pragma loop name foo3#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(j)
		for (j=0; j<=9999; j+=1)
		{
			a[i][j]+=a[i+1][j];
		}
	}
	_ret_val_0=0;
	return _ret_val_0;
}

void foo4(int i, int j)
{
	double a[10000][10000];
	#pragma cetus private(i, j) 
	#pragma loop name foo4#0 
	for (i=0; i<=9998; i+=1)
	{
		#pragma cetus private(j) 
		#pragma loop name foo4#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(j)
		for (j=0; j<=9999; j+=1)
		{
			a[i][j]+=a[i+1][j];
		}
	}
	return ;
}

/*

 Coefficient subscript

*/
void foo5()
{
	int i;
	int a[1000];
	#pragma cetus private(i) 
	#pragma loop name foo5#0 
	for (i=0; i<=998; i+=1)
	{
		a[(2*i)+1]=(a[i]+1);
	}
	/* a[i+1]=a[i]+1; */
	return ;
}

/*

 Complex condition

*/
void goo(int numAB)
{
	double * c;
	double * bufLoc;
	int k_nom_22;
	#pragma cetus private(k_nom_22) 
	#pragma loop name goo#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+((9L*numAB)*numAB)))) private(k_nom_22)
	for (k_nom_22=0; k_nom_22<=(((numAB*numAB)*3)-1); k_nom_22+=1)
	{
		bufLoc[k_nom_22]=c[k_nom_22];
	}
	return ;
}

/*

 Deep distance

*/
void foo6()
{
	int i;
	int j;
int b[1000][1000] = {234};
	#pragma cetus private(i, j) 
	#pragma loop name foo6#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j)
	for (i=1; i<=99; i+=1)
	{
		#pragma cetus private(j) 
		#pragma loop name foo6#0#0 
		for (j=1; j<=99; j+=1)
		{
			b[i][j]=(b[i][j-1]+1);
		}
	}
	return ;
}

/*

 Doall

*/
void foo7()
{
	int i;
	int a[10000];
	#pragma cetus private(i) 
	#pragma loop name foo7#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i)
	for (i=0; i<=9999; i+=1)
	{
		a[i]=(a[i]+1);
	}
	return ;
}

/*

 first private

*/
void foo8(double * o1, double * c, int * * idx, int len)
{
	int i;
	#pragma cetus private() 
	#pragma loop name foo8#0 
	#pragma cetus reduction(+: o1[lidx[ii]]) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(30L*len)))) reduction(+: o1[lidx[ii]])
	for (i=0; i<=(len-1); i+=1)
	{
		int ii;
		const int * lidx = idx[i];
		double volnew_o8 = 0.5*c[i];
		#pragma cetus private() 
		#pragma loop name foo8#0#0 
		#pragma cetus reduction(+: o1[lidx[ii]]) 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for reduction(+: o1[lidx[ii]])
		*/
		for (ii=0; ii<=5; ii+=1)
		{
			int llidx = lidx[ii];
			o1[lidx[ii]]+=volnew_o8;
		}
	}
	return ;
}

/*

 first private 3

*/
void error_check()
{
	int i;
	int j;
	int n;
	int m;
	double xx;
	double yy;
	double temp;
	double error;
	double dx = 2.0/(n-1);
	double dy = 2.0/(m-1);
	double u[200][200];
	error=0.0;
	#pragma cetus private(i, j, temp, xx, yy) 
	#pragma loop name error_check#0 
	#pragma cetus reduction(+: error) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(3L*n))+((6L*m)*n)))) private(i, j, temp, xx, yy) reduction(+: error)
	for (i=0; i<=(n-1); i+=1)
	{
		#pragma cetus private(j, temp, xx, yy) 
		#pragma loop name error_check#0#0 
		#pragma cetus reduction(+: error) 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(6L*m)))) private(j, temp, xx, yy) reduction(+: error)
		for (j=0; j<=(m-1); j+=1)
		{
			xx=(( - 1.0)+(dx*(i-1)));
			yy=(( - 1.0)+(dy*(j-1)));
			temp=(u[i][j]-((1.0-(xx*xx))*(1.0-(yy*yy))));
			error=(error+(temp*temp));
		}
	}
	error=(sqrt(error)/(n*m));
	return ;
}

/*

 foo

*/
void foo9(double o1[], double c[], int len)
{
	int i;
	#pragma cetus private() 
	#pragma loop name foo9#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(4L*len))))
	for (i=0; i<=(len-1); i+=1)
	{
		double volnew_o8 = 0.5*c[i];
		o1[i]=volnew_o8;
	}
	return ;
}

/*

 for i j

*/
int ax[10000];
int bx[10000];
void foo10()
{
	int i;
	int j = 0;
	int i_ub = 10000;
	#pragma cetus private(i) 
	#pragma loop name foo10#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i)
	for (i=0; i<=(i_ub-1); i ++ )
	{
		bx[i+j]=ax[i];
	}
	j+=10000;
	return ;
}

/*

 function call

*/
void foo11(int istart, int iend, real8 * a, real8 * b, real8 * c, int k, real8 * l, real8 * m, real8 * n, real8 * o, real8 * p)
{
	int i;
	#pragma cetus private(k) 
	#pragma loop name foo11#0 
	for (i=istart; i<=(iend-1); i+=1)
	{
		real8 s[3];
		real8 afi = a[i];
		real8 bfi = b[i];
		int k;
		OtherFunc(k, l, m, n, o, p, afi, bfi, s);
		#pragma cetus lastprivate(k) 
		#pragma loop name foo11#0#0 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for lastprivate(k)
		*/
		for (k=0; k<=2; k+=1)
		{
			c[(3*i)+k]=s[k];
		}
	}
	return ;
}

/*

 Global arrays

*/
double u[1000][1000];
double f[1000][1000];
int n;
int m;
void initialize()
{
	int i;
	int j;
	int xx;
	double dx = 2.0/(n-1);
	n=1000;
	m=1000;
	#pragma cetus private(i, j, xx) 
	#pragma loop name initialize#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j, xx)
	for (i=0; i<=(n-1); i+=1)
	{
		#pragma cetus private(j, xx) 
		#pragma loop name initialize#0#0 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(j, xx)
		*/
		for (j=0; j<=(m-1); j+=1)
		{
			xx=((int)(( - 1.0)+(dx*(i-1))));
			u[i][j]=0.0;
			f[i][j]=(( - 1.0)*(1.0-(xx*xx)));
		}
	}
	return ;
}

/*

 if and for

*/
void foo12(int j)
{
	int i;
	int a[10000];
	if (j!=( - 1))
	{
		#pragma cetus private(i) 
		#pragma loop name foo12#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(i)
		for (i=0; i<=9999; i+=1)
		{
			a[i]=(a[i]+1);
		}
	}
	return ;
}

/*

 Indirect index

*/
double eps[1000];
int zoneset[1000];
void StressCheckEpsFail1(double eps_failure_model)
{
	int i;
	int index;
	#pragma cetus private(i, index) 
	#pragma loop name StressCheckEpsFail1#0 
	for (i=0; i<=999; i+=1)
	{
		index=zoneset[i];
		eps[zoneset[i]]=(eps_failure_model*1.01);
		eps[zoneset[i]]=1.01;
	}
	return ;
}

void StressCheckEpsFaili2(double eps_failure_model)
{
	int i;
	int index1;
	#pragma cetus private(i, index2) 
	#pragma loop name StressCheckEpsFaili2#0 
	for (i=0; i<=999; i+=1)
	{
		int index2 = index1;
		index1=zoneset[i];
		eps[zoneset[i]]=(eps_failure_model*1.01);
		eps[zoneset[i]]=1.01;
	}
	return ;
}

void foo13()
{
	int n = 1000;
	int m = 1000;
	double b[n][m];
	int i;
	int j;
	int index;
	int zoneset[m];
	#pragma cetus private(i, index, j) 
	#pragma loop name foo13#0 
	for (i=0; i<=(n-1); i+=1)
	{
		#pragma cetus private(index, j) 
		#pragma loop name foo13#0#0 
		for (j=0; j<=(m-1); j+=1)
		{
			index=zoneset[j];
			b[i][zoneset[j]]=b[i-1][index-1];
		}
	}
	return ;
}

/*

 indirect index transfered

*/
void foo14(int * indexSet, int N, int ax)
{
	double * xa3[N];
	int idx;
	#pragma cetus private(idx) 
	#pragma loop name foo14#0 
	#pragma cetus reduction(+: xa3[indexSet[idx]]) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(4L*N)))) private(idx) reduction(+: xa3[indexSet[idx]])
	for (idx=0; idx<=(N-1); idx+=1)
	{
		xa3[indexSet[idx]]+=ax;
		xa3[indexSet[idx]]+=ax;
	}
	return ;
}

void foo15(int * indexSet, int N, int ax)
{
	double * xa3[N];
	int idx;
	#pragma cetus private() 
	#pragma loop name foo15#0 
	#pragma cetus reduction(+: xa3[indexSet[idx]]) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(5L*N)))) reduction(+: xa3[indexSet[idx]])
	for (idx=0; idx<=(N-1); idx+=1)
	{
		const int i = indexSet[idx];
		xa3[indexSet[idx]]+=ax;
		xa3[indexSet[idx]]+=ax;
	}
	return ;
}

/*

 inner only loop

*/
void foo16()
{
	int n = 10000;
	int m = 10000;
	double b[n][m];
	int i;
	int j;
	#pragma cetus private(i, j) 
	#pragma loop name foo16#0 
	for (i=0; i<=(n-1); i+=1)
	{
		#pragma cetus private(j) 
		#pragma loop name foo16#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(j)
		for (j=0; j<=(m-1); j+=1)
		{
			b[i][j]=b[i-1][j-1];
		}
	}
	return ;
}

/*

 last private

*/
void foo17()
{
	int i;
	int x;
	#pragma cetus private(i, x) 
	#pragma loop name foo17#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, x)
	for (i=0; i<=99999; i+=1)
	{
		x=i;
	}
	return ;
}

void foo18()
{
	int a[1000];
	int i;
	int x = 10;
	#pragma cetus private(i) 
	#pragma loop name foo18#0 
	for (i=0; i<=999; i+=1)
	{
		a[i]=x;
		x=i;
	}
	return ;
}

/*

 liveness test

*/
void foo19(real8 * y, real8 * d__, real8 * d11, real8 * d12, real8 * d13, real8 * d22, real8 * d23, real8 * d33, real8 * m, int * nell, real8 * p, int t, int flagB, int flagA, int ub)
{
	int l;
	int nel;
	int t1 = t-1;
	if (flagB==0)
	{
		#pragma cetus private() 
		#pragma loop name foo19#0 
		#pragma cetus reduction(+: y[l36+0], y[l36+10], y[l36+11], y[l36+12], y[l36+13], y[l36+14], y[l36+15], y[l36+16], y[l36+17], y[l36+18], y[l36+19], y[l36+1], y[l36+20], y[l36+21], y[l36+22], y[l36+23], y[l36+24], y[l36+25], y[l36+26], y[l36+27], y[l36+28], y[l36+29], y[l36+2], y[l36+30], y[l36+31], y[l36+32], y[l36+33], y[l36+34], y[l36+35], y[l36+3], y[l36+4], y[l36+5], y[l36+6], y[l36+7], y[l36+8], y[l36+9]) 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(65L*ub)))) reduction(+: y[l36+0], y[l36+10], y[l36+11], y[l36+12], y[l36+13], y[l36+14], y[l36+15], y[l36+16], y[l36+17], y[l36+18], y[l36+19], y[l36+1], y[l36+20], y[l36+21], y[l36+22], y[l36+23], y[l36+24], y[l36+25], y[l36+26], y[l36+27], y[l36+28], y[l36+29], y[l36+2], y[l36+30], y[l36+31], y[l36+32], y[l36+33], y[l36+34], y[l36+35], y[l36+3], y[l36+4], y[l36+5], y[l36+6], y[l36+7], y[l36+8], y[l36+9])
		for (l=0; l<=(ub-1); l+=1)
		{
			int l8 = l*8;
			int l36 = l*36;
			real8 h12 = m[((l8+0)*4)+1];
			real8 h13 = m[((l8+0)*4)+2];
			real8 h14 = m[((l8+0)*4)+3];
			real8 h22 = m[((l8+1)*4)+1];
			real8 h23 = m[((l8+1)*4)+2];
			real8 h24 = m[((l8+1)*4)+3];
			real8 h32 = m[((l8+2)*4)+1];
			real8 h33 = m[((l8+2)*4)+2];
			real8 h34 = m[((l8+2)*4)+3];
			real8 h42 = m[((l8+3)*4)+1];
			real8 h43 = m[((l8+3)*4)+2];
			real8 h44 = m[((l8+3)*4)+3];
			real8 h52 = m[((l8+4)*4)+1];
			real8 h53 = m[((l8+4)*4)+2];
			real8 h54 = m[((l8+4)*4)+3];
			real8 h62 = m[((l8+5)*4)+1];
			real8 h63 = m[((l8+5)*4)+2];
			real8 h64 = m[((l8+5)*4)+3];
			real8 h72 = m[((l8+6)*4)+1];
			real8 h73 = m[((l8+6)*4)+2];
			real8 h74 = m[((l8+6)*4)+3];
			real8 h82 = m[((l8+7)*4)+1];
			real8 h83 = m[((l8+7)*4)+2];
			real8 h84 = m[((l8+7)*4)+3];
			real8 ddd = d__[l];
			y[l36+0]+=(ddd*(((h12*h12)+(h13*h13))+(h14*h14)));
			y[l36+1]+=(ddd*(((h12*h22)+(h13*h23))+(h14*h24)));
			y[l36+2]+=(ddd*(((h22*h22)+(h23*h23))+(h24*h24)));
			y[l36+3]+=(ddd*(((h12*h32)+(h13*h33))+(h14*h34)));
			y[l36+4]+=(ddd*(((h22*h32)+(h23*h33))+(h24*h34)));
			y[l36+5]+=(ddd*(((h32*h32)+(h33*h33))+(h34*h34)));
			y[l36+6]+=(ddd*(((h12*h42)+(h13*h43))+(h14*h44)));
			y[l36+7]+=(ddd*(((h22*h42)+(h23*h43))+(h24*h44)));
			y[l36+8]+=(ddd*(((h32*h42)+(h33*h43))+(h34*h44)));
			y[l36+9]+=(ddd*(((h42*h42)+(h43*h43))+(h44*h44)));
			y[l36+10]+=(ddd*(((h12*h52)+(h13*h53))+(h14*h54)));
			y[l36+11]+=(ddd*(((h22*h52)+(h23*h53))+(h24*h54)));
			y[l36+12]+=(ddd*(((h32*h52)+(h33*h53))+(h34*h54)));
			y[l36+13]+=(ddd*(((h42*h52)+(h43*h53))+(h44*h54)));
			y[l36+14]+=(ddd*(((h52*h52)+(h53*h53))+(h54*h54)));
			y[l36+15]+=(ddd*(((h12*h62)+(h13*h63))+(h14*h64)));
			y[l36+16]+=(ddd*(((h22*h62)+(h23*h63))+(h24*h64)));
			y[l36+17]+=(ddd*(((h32*h62)+(h33*h63))+(h34*h64)));
			y[l36+18]+=(ddd*(((h42*h62)+(h43*h63))+(h44*h64)));
			y[l36+19]+=(ddd*(((h52*h62)+(h53*h63))+(h54*h64)));
			y[l36+20]+=(ddd*(((h62*h62)+(h63*h63))+(h64*h64)));
			y[l36+21]+=(ddd*(((h12*h72)+(h13*h73))+(h14*h74)));
			y[l36+22]+=(ddd*(((h22*h72)+(h23*h73))+(h24*h74)));
			y[l36+23]+=(ddd*(((h32*h72)+(h33*h73))+(h34*h74)));
			y[l36+24]+=(ddd*(((h42*h72)+(h43*h73))+(h44*h74)));
			y[l36+25]+=(ddd*(((h52*h72)+(h53*h73))+(h54*h74)));
			y[l36+26]+=(ddd*(((h62*h72)+(h63*h73))+(h64*h74)));
			y[l36+27]+=(ddd*(((h72*h72)+(h73*h73))+(h74*h74)));
			y[l36+28]+=(ddd*(((h12*h82)+(h13*h83))+(h14*h84)));
			y[l36+29]+=(ddd*(((h22*h82)+(h23*h83))+(h24*h84)));
			y[l36+30]+=(ddd*(((h32*h82)+(h33*h83))+(h34*h84)));
			y[l36+31]+=(ddd*(((h42*h82)+(h43*h83))+(h44*h84)));
			y[l36+32]+=(ddd*(((h52*h82)+(h53*h83))+(h54*h84)));
			y[l36+33]+=(ddd*(((h62*h82)+(h63*h83))+(h64*h84)));
			y[l36+34]+=(ddd*(((h72*h82)+(h73*h83))+(h74*h84)));
			y[l36+35]+=(ddd*(((h82*h82)+(h83*h83))+(h84*h84)));
		}
		if (flagA>0)
		{
			#pragma cetus private(nel) 
			#pragma loop name foo19#1 
			#pragma cetus reduction(+: p[nell[l]]) 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+(8L*ub)))) private(nel) reduction(+: p[nell[l]])
			for (l=0; l<=(ub-1); l+=1)
			{
				int l8 = l*8;
				real8 h1 = m[((t1+l8)*4)+1];
				real8 h2 = m[((t1+l8)*4)+2];
				real8 h3 = m[((t1+l8)*4)+3];
				nel=nell[l];
				p[nell[l]]+=((d__[l]*64.0)*(((h1*h1)+(h2*h2))+(h3*h3)));
			}
		}
	}
	else
	{
		#pragma cetus private() 
		#pragma loop name foo19#2 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(70L*ub))))
		for (l=0; l<=(ub-1); l+=1)
		{
			int l8 = l*8;
			int l36 = l*36;
			real8 d_11 = d11[l];
			real8 d_12 = d12[l];
			real8 d_13 = d13[l];
			real8 d_22 = d22[l];
			real8 d_23 = d23[l];
			real8 d_33 = d33[l];
			real8 h12 = m[((l8+0)*4)+1];
			real8 h13 = m[((l8+0)*4)+2];
			real8 h14 = m[((l8+0)*4)+3];
			real8 h22 = m[((l8+1)*4)+1];
			real8 h23 = m[((l8+1)*4)+2];
			real8 h24 = m[((l8+1)*4)+3];
			real8 h32 = m[((l8+2)*4)+1];
			real8 h33 = m[((l8+2)*4)+2];
			real8 h34 = m[((l8+2)*4)+3];
			real8 h42 = m[((l8+3)*4)+1];
			real8 h43 = m[((l8+3)*4)+2];
			real8 h44 = m[((l8+3)*4)+3];
			real8 h52 = m[((l8+4)*4)+1];
			real8 h53 = m[((l8+4)*4)+2];
			real8 h54 = m[((l8+4)*4)+3];
			real8 h62 = m[((l8+5)*4)+1];
			real8 h63 = m[((l8+5)*4)+2];
			real8 h64 = m[((l8+5)*4)+3];
			real8 h72 = m[((l8+6)*4)+1];
			real8 h73 = m[((l8+6)*4)+2];
			real8 h74 = m[((l8+6)*4)+3];
			real8 h82 = m[((l8+7)*4)+1];
			real8 h83 = m[((l8+7)*4)+2];
			real8 h84 = m[((l8+7)*4)+3];
			y[l36+0]=(((y[l36+0]+(h12*(((d_11*h12)+(d_12*h13))+(d_13*h14))))+(h13*(((d_12*h12)+(d_22*h13))+(d_23*h14))))+(h14*(((d_13*h12)+(d_23*h13))+(d_33*h14))));
			y[l36+1]=(((y[l36+1]+(h22*(((d_11*h12)+(d_12*h13))+(d_13*h14))))+(h23*(((d_12*h12)+(d_22*h13))+(d_23*h14))))+(h24*(((d_13*h12)+(d_23*h13))+(d_33*h14))));
			y[l36+2]=(((y[l36+2]+(h22*(((d_11*h22)+(d_12*h23))+(d_13*h24))))+(h23*(((d_12*h22)+(d_22*h23))+(d_23*h24))))+(h24*(((d_13*h22)+(d_23*h23))+(d_33*h24))));
			y[l36+3]=(((y[l36+3]+(h32*(((d_11*h12)+(d_12*h13))+(d_13*h14))))+(h33*(((d_12*h12)+(d_22*h13))+(d_23*h14))))+(h34*(((d_13*h12)+(d_23*h13))+(d_33*h14))));
			y[l36+4]=(((y[l36+4]+(h32*(((d_11*h22)+(d_12*h23))+(d_13*h24))))+(h33*(((d_12*h22)+(d_22*h23))+(d_23*h24))))+(h34*(((d_13*h22)+(d_23*h23))+(d_33*h24))));
			y[l36+5]=(((y[l36+5]+(h32*(((d_11*h32)+(d_12*h33))+(d_13*h34))))+(h33*(((d_12*h32)+(d_22*h33))+(d_23*h34))))+(h34*(((d_13*h32)+(d_23*h33))+(d_33*h34))));
			y[l36+6]=(((y[l36+6]+(h42*(((d_11*h12)+(d_12*h13))+(d_13*h14))))+(h43*(((d_12*h12)+(d_22*h13))+(d_23*h14))))+(h44*(((d_13*h12)+(d_23*h13))+(d_33*h14))));
			y[l36+7]=(((y[l36+7]+(h42*(((d_11*h22)+(d_12*h23))+(d_13*h24))))+(h43*(((d_12*h22)+(d_22*h23))+(d_23*h24))))+(h44*(((d_13*h22)+(d_23*h23))+(d_33*h24))));
			y[l36+8]=(((y[l36+8]+(h42*(((d_11*h32)+(d_12*h33))+(d_13*h34))))+(h43*(((d_12*h32)+(d_22*h33))+(d_23*h34))))+(h44*(((d_13*h32)+(d_23*h33))+(d_33*h34))));
			y[l36+9]=(((y[l36+9]+(h42*(((d_11*h42)+(d_12*h43))+(d_13*h44))))+(h43*(((d_12*h42)+(d_22*h43))+(d_23*h44))))+(h44*(((d_13*h42)+(d_23*h43))+(d_33*h44))));
			y[l36+10]=(((y[l36+10]+(h52*(((d_11*h12)+(d_12*h13))+(d_13*h14))))+(h53*(((d_12*h12)+(d_22*h13))+(d_23*h14))))+(h54*(((d_13*h12)+(d_23*h13))+(d_33*h14))));
			y[l36+11]=(((y[l36+11]+(h52*(((d_11*h22)+(d_12*h23))+(d_13*h24))))+(h53*(((d_12*h22)+(d_22*h23))+(d_23*h24))))+(h54*(((d_13*h22)+(d_23*h23))+(d_33*h24))));
			y[l36+12]=(((y[l36+12]+(h52*(((d_11*h32)+(d_12*h33))+(d_13*h34))))+(h53*(((d_12*h32)+(d_22*h33))+(d_23*h34))))+(h54*(((d_13*h32)+(d_23*h33))+(d_33*h34))));
			y[l36+13]=(((y[l36+13]+(h52*(((d_11*h42)+(d_12*h43))+(d_13*h44))))+(h53*(((d_12*h42)+(d_22*h43))+(d_23*h44))))+(h54*(((d_13*h42)+(d_23*h43))+(d_33*h44))));
			y[l36+14]=(((y[l36+14]+(h52*(((d_11*h52)+(d_12*h53))+(d_13*h54))))+(h53*(((d_12*h52)+(d_22*h53))+(d_23*h54))))+(h54*(((d_13*h52)+(d_23*h53))+(d_33*h54))));
			y[l36+15]=(((y[l36+15]+(h62*(((d_11*h12)+(d_12*h13))+(d_13*h14))))+(h63*(((d_12*h12)+(d_22*h13))+(d_23*h14))))+(h64*(((d_13*h12)+(d_23*h13))+(d_33*h14))));
			y[l36+16]=(((y[l36+16]+(h62*(((d_11*h22)+(d_12*h23))+(d_13*h24))))+(h63*(((d_12*h22)+(d_22*h23))+(d_23*h24))))+(h64*(((d_13*h22)+(d_23*h23))+(d_33*h24))));
			y[l36+17]=(((y[l36+17]+(h62*(((d_11*h32)+(d_12*h33))+(d_13*h34))))+(h63*(((d_12*h32)+(d_22*h33))+(d_23*h34))))+(h64*(((d_13*h32)+(d_23*h33))+(d_33*h34))));
			y[l36+18]=(((y[l36+18]+(h62*(((d_11*h42)+(d_12*h43))+(d_13*h44))))+(h63*(((d_12*h42)+(d_22*h43))+(d_23*h44))))+(h64*(((d_13*h42)+(d_23*h43))+(d_33*h44))));
			y[l36+19]=(((y[l36+19]+(h62*(((d_11*h52)+(d_12*h53))+(d_13*h54))))+(h63*(((d_12*h52)+(d_22*h53))+(d_23*h54))))+(h64*(((d_13*h52)+(d_23*h53))+(d_33*h54))));
			y[l36+20]=(((y[l36+20]+(h62*(((d_11*h62)+(d_12*h63))+(d_13*h64))))+(h63*(((d_12*h62)+(d_22*h63))+(d_23*h64))))+(h64*(((d_13*h62)+(d_23*h63))+(d_33*h64))));
			y[l36+21]=(((y[l36+21]+(h72*(((d_11*h12)+(d_12*h13))+(d_13*h14))))+(h73*(((d_12*h12)+(d_22*h13))+(d_23*h14))))+(h74*(((d_13*h12)+(d_23*h13))+(d_33*h14))));
			y[l36+22]=(((y[l36+22]+(h72*(((d_11*h22)+(d_12*h23))+(d_13*h24))))+(h73*(((d_12*h22)+(d_22*h23))+(d_23*h24))))+(h74*(((d_13*h22)+(d_23*h23))+(d_33*h24))));
			y[l36+23]=(((y[l36+23]+(h72*(((d_11*h32)+(d_12*h33))+(d_13*h34))))+(h73*(((d_12*h32)+(d_22*h33))+(d_23*h34))))+(h74*(((d_13*h32)+(d_23*h33))+(d_33*h34))));
			y[l36+24]=(((y[l36+24]+(h72*(((d_11*h42)+(d_12*h43))+(d_13*h44))))+(h73*(((d_12*h42)+(d_22*h43))+(d_23*h44))))+(h74*(((d_13*h42)+(d_23*h43))+(d_33*h44))));
			y[l36+25]=(((y[l36+25]+(h72*(((d_11*h52)+(d_12*h53))+(d_13*h54))))+(h73*(((d_12*h52)+(d_22*h53))+(d_23*h54))))+(h74*(((d_13*h52)+(d_23*h53))+(d_33*h54))));
			y[l36+26]=(((y[l36+26]+(h72*(((d_11*h62)+(d_12*h63))+(d_13*h64))))+(h73*(((d_12*h62)+(d_22*h63))+(d_23*h64))))+(h74*(((d_13*h62)+(d_23*h63))+(d_33*h64))));
			y[l36+27]=(((y[l36+27]+(h72*(((d_11*h72)+(d_12*h73))+(d_13*h74))))+(h73*(((d_12*h72)+(d_22*h73))+(d_23*h74))))+(h74*(((d_13*h72)+(d_23*h73))+(d_33*h74))));
			y[l36+28]=(((y[l36+28]+(h82*(((d_11*h12)+(d_12*h13))+(d_13*h14))))+(h83*(((d_12*h12)+(d_22*h13))+(d_23*h14))))+(h84*(((d_13*h12)+(d_23*h13))+(d_33*h14))));
			y[l36+29]=(((y[l36+29]+(h82*(((d_11*h22)+(d_12*h23))+(d_13*h24))))+(h83*(((d_12*h22)+(d_22*h23))+(d_23*h24))))+(h84*(((d_13*h22)+(d_23*h23))+(d_33*h24))));
			y[l36+30]=(((y[l36+30]+(h82*(((d_11*h32)+(d_12*h33))+(d_13*h34))))+(h83*(((d_12*h32)+(d_22*h33))+(d_23*h34))))+(h84*(((d_13*h32)+(d_23*h33))+(d_33*h34))));
			y[l36+31]=(((y[l36+31]+(h82*(((d_11*h42)+(d_12*h43))+(d_13*h44))))+(h83*(((d_12*h42)+(d_22*h43))+(d_23*h44))))+(h84*(((d_13*h42)+(d_23*h43))+(d_33*h44))));
			y[l36+32]=(((y[l36+32]+(h82*(((d_11*h52)+(d_12*h53))+(d_13*h54))))+(h83*(((d_12*h52)+(d_22*h53))+(d_23*h54))))+(h84*(((d_13*h52)+(d_23*h53))+(d_33*h54))));
			y[l36+33]=(((y[l36+33]+(h82*(((d_11*h62)+(d_12*h63))+(d_13*h64))))+(h83*(((d_12*h62)+(d_22*h63))+(d_23*h64))))+(h84*(((d_13*h62)+(d_23*h63))+(d_33*h64))));
			y[l36+34]=(((y[l36+34]+(h82*(((d_11*h72)+(d_12*h73))+(d_13*h74))))+(h83*(((d_12*h72)+(d_22*h73))+(d_23*h74))))+(h84*(((d_13*h72)+(d_23*h73))+(d_33*h74))));
			y[l36+35]=(((y[l36+35]+(h82*(((d_11*h82)+(d_12*h83))+(d_13*h84))))+(h83*(((d_12*h82)+(d_22*h83))+(d_23*h84))))+(h84*(((d_13*h82)+(d_23*h83))+(d_33*h84))));
		}
		if (flagA>0)
		{
			#pragma cetus private(nel) 
			#pragma loop name foo19#3 
			#pragma cetus reduction(+: p[nell[l]]) 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+(8L*ub)))) private(nel) reduction(+: p[nell[l]])
			for (l=0; l<=(ub-1); l+=1)
			{
				int l8 = l*8;
				real8 h1 = m[((t1+l8)*4)+1];
				real8 h2 = m[((t1+l8)*4)+2];
				real8 h3 = m[((t1+l8)*4)+3];
				nel=nell[l];
				p[nell[l]]+=((((h1*(((d11[l]*h1)+((d12[l]*2.0)*h2))+((d13[l]*2.0)*h3)))+(h2*((d22[l]*h2)+((d23[l]*2.0)*h3))))+((h3*d33[l])*h3))*64.0);
			}
		}
	}
	return ;
}

/*

 Matrix vector multiplication

*/
int mmm()
{
	int i;
	int j;
	int k;
	double a[1000][1000];
	double v[1000];
	double v_out[1000];
	int _ret_val_0;
	#pragma cetus private(i, j) 
	#pragma loop name mmm#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j)
	for (i=0; i<=999; i+=1)
	{
		float sum = 0.0;
		#pragma cetus private(j) 
		#pragma loop name mmm#0#0 
		#pragma cetus reduction(+: sum) 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(j) reduction(+: sum)
		*/
		for (j=0; j<=999; j+=1)
		{
			sum+=(a[i][j]*v[j]);
		}
		v_out[i]=sum;
	}
	_ret_val_0=0;
	return _ret_val_0;
}

/*

 minus minus

*/
void foo20(int numNodes, int numNodes2, int * x, int * nodelist)
{
	int j;
	#pragma cetus private(j) 
	#pragma loop name foo20#0 
	for (j=(numNodes-1); j>=0; j+=( - 1))
	{
		if (x[j]<=0.0)
		{
			numNodes2 -- ;
			nodelist[j]=nodelist[numNodes2];
			nodelist[numNodes2]=j;
		}
	}
	return ;
}

/*

 Outer only loop

*/
void foo21()
{
	int n = 1000;
	int m = 1000;
	double b[n][m];
	int i;
	int j;
	#pragma cetus private(i, j) 
	#pragma loop name foo21#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j)
	for (i=0; i<=(n-1); i+=1)
	{
		#pragma cetus private(j) 
		#pragma loop name foo21#0#0 
		for (j=1; j<=(m-1); j+=1)
		{
			b[i][j]=b[i][j-1];
		}
	}
	return ;
}

/*

 Output dependence

*/
void foo22()
{
	int i;
	int x;
	#pragma cetus private(i, x) 
	#pragma loop name foo22#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, x)
	for (i=0; i<=99999; i+=1)
	{
		x=i;
	}
	return ;
}

/*

 Output dependence2

*/
void foo23()
{
	int i;
	int x;
	int y;
	#pragma cetus private(i, x, y) 
	#pragma loop name foo23#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, x, y)
	for (i=0; i<=99999; i+=1)
	{
		x=i;
		y=i;
		y=(i+1);
	}
	return ;
}

/*

 Output dependence3

*/
void foo24()
{
	int i;
	int j;
	int x;
	int y;
	#pragma cetus private(i, j, x, y) 
	#pragma loop name foo24#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j, x, y)
	for (i=0; i<=9999; i+=1)
	{
		#pragma cetus private(j, x, y) 
		#pragma loop name foo24#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(j, x, y)
		for (j=0; j<=9999; j+=1)
		{
			x=i;
			y=x;
			y=(i+j);
			y=(y+1);
		}
	}
	return ;
}

/*

 Plus assign

*/
void foo25()
{
	int i;
	int j;
	double a[10000][10000];
	#pragma cetus private(i, j) 
	#pragma loop name foo25#0 
	for (i=0; i<=9998; i+=1)
	{
		#pragma cetus private(j) 
		#pragma loop name foo25#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(j)
		for (j=0; j<=9999; j+=1)
		{
			a[i][j]+=a[i+1][j];
		}
	}
	return ;
}

/*

 Plus plus op

*/
int HighPassFilter(int * input, int inLen, int * output, int threshold)
{
	int outLen = 0;
	int i;
	#pragma cetus private(i) 
	#pragma loop name HighPassFilter#0 
	for (i=0; i<=(inLen-1); i+=1)
	{
		if (input[i]>threshold)
		{
			output[outLen ++ ]=input[i];
		}
	}
	return outLen;
}

/*

 Pointer dereference

*/
int * nd_array[10000];
int * gr_array[10000];
void foo26(int rlenmix)
{
	int i;
	#pragma cetus private() 
	#pragma loop name foo26#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(6L*rlenmix))))
	for (i=0; i<=(rlenmix-1); i+=1)
	{
		int nn =  * (nd_array[1]+i);
		int gg =  * (gr_array[1]+i);
		nn=( * (nd_array[2]+i));
		gg=( * (gr_array[2]+i));
	}
	return ;
}

/*

 Pointers

*/
void foo27(double * x, int jp, int begin, int end, double rh1)
{
	double * x1;
	double * x2;
	int i;
	x1=x;
	x2=(x1+jp);
	#pragma cetus private(i) 
	#pragma loop name foo27#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(-4L*begin))+(4L*end)))) private(i)
	for (i=begin; i<=(end-1); i+=1)
	{
		x1[i]+=rh1;
		x2[i]-=rh1;
	}
	return ;
}

/*

 Private

*/
int g;
void foo28()
{
	int i;
	int x;
	int a[10000];
	int b[10000];
	#pragma cetus private(x) 
	#pragma loop name foo28#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(x)
	for (i=0; i<=9999; i+=1)
	{
		int y = i+1;
		/*   g = y; */
		x=(a[i]+g);
		/* b[i]=x+1+y; */
	}
	return ;
}

/*

 Reduction

*/
int a[10000];
int sum;
void foo29()
{
	int i;
	int sum2;
	int xx;
	int yy;
	int zz;
	sum=0;
	#pragma cetus private(i) 
	#pragma loop name foo29#0 
	#pragma cetus reduction(*: zz) reduction(+: sum) 
	#pragma cetus parallel 
	#pragma omp parallel for private(i) reduction(*: zz+: sum)
	for (i=0; i<=9999; i+=1)
	{
		a[i]=i;
		sum=(a[i]+sum);
		zz*=a[i];
	}
	yy+=-10000;
	xx+=10000;
	sum2=(((sum+xx)+yy)+zz);
	a[1]=1;
	return ;
}

/*

 Reduction 2

*/
float uu[1000][1000];
float foo30()
{
	int i;
	int j;
	float temp;
	float error;
	#pragma cetus private(i, j, temp) 
	#pragma loop name foo30#0 
	#pragma cetus reduction(+: error) 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j, temp) reduction(+: error)
	for (i=0; i<=999; i+=1)
	{
		#pragma cetus private(j, temp) 
		#pragma loop name foo30#0#0 
		#pragma cetus reduction(+: error) 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(j, temp) reduction(+: error)
		*/
		for (j=0; j<=999; j+=1)
		{
			temp=uu[i][j];
			error=(error+(temp*temp));
		}
	}
	return error;
}

/*

 Reduction fake

*/
extern int bar();
int a[10000];
int sum;
void foo31()
{
	int i;
	int sum2;
	int xx;
	int yy;
	int zz;
	sum=0;
	#pragma cetus private(i) 
	#pragma loop name foo31#0 
	for (i=0; i<=9999; i+=1)
	{
		a[i]=i;
		sum=((a[i]+sum)+bar());
		/*    sum = a[i]+ sum ;     */
	}
	sum2=sum;
	a[1]=1;
	return ;
}

/*

 Reduction max

*/
double aa[10000];
void foo32()
{
	double max_val =  - 1.0E99;
	double min_val = 1.0E99;
	int i;
	#pragma cetus private(i) 
	#pragma loop name foo32#0 
	for (i=0; i<=9999; i+=1)
	{
		if (aa[i]>max_val)
		{
			max_val=aa[i];
		}
		if (aa[i]<min_val)
		{
			min_val=aa[i];
		}
	}
	return ;
}

/*

 Regression

*/
void foo33(real8 * a, real8 * b, real8 * c, real8 * d, int len)
{
	int icol;
	int jrow;
	int l;
	#pragma cetus private(icol, l8) 
	#pragma loop name foo33#0 
	#pragma cetus reduction(+: a[icol+((jrow+l8)*8)]) 
	for (l=0; l<=(len-1); l+=1)
	{
		int l8 = l*8;
		real8 e = d[(l*3)+0];
		real8 f = d[(l*3)+1];
		real8 g = d[(l*3)+2];
		real8 h = b[l];
		real8 tmp[10000];
		#pragma cetus private(icol) 
		#pragma loop name foo33#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(icol)
		for (icol=0; icol<=9999; icol+=1)
		{
			tmp[icol]=(((e*c[((icol+l8)*4)+1])+(f*c[((icol+l8)*4)+2]))+(g*c[((icol+l8)*4)+3]));
		}
		#pragma cetus private(icol) 
		#pragma loop name foo33#0#1 
		#pragma cetus reduction(+: a[icol+((jrow+l8)*8)]) 
		#pragma cetus parallel 
		#pragma omp parallel for private(icol) reduction(+: a[icol+((jrow+l8)*8)])
		for (jrow=0; jrow<=9999; jrow+=1)
		{
			real8 hj1 = h*c[(jrow+l8)*4];
			#pragma cetus private(icol) 
			#pragma loop name foo33#0#1#0 
			#pragma cetus reduction(+: a[icol+((jrow+l8)*8)]) 
			#pragma cetus parallel 
			#pragma omp parallel for private(icol) reduction(+: a[icol+((jrow+l8)*8)])
			for (icol=0; icol<=9999; icol+=1)
			{
				a[icol+((jrow+l8)*8)]+=(hj1*tmp[icol]);
			}
		}
	}
	return ;
}

/*

 Scalar-anti

*/
void foo34()
{
	int i;
	int tmp;
	tmp=10;
	#pragma cetus private(i) 
	#pragma loop name foo34#0 
	for (i=0; i<=9999; i+=1)
	{
		a[i]=tmp;
		tmp=(a[i]+i);
	}
	return ;
}

void foo35()
{
	int i;
	int tmp;
	tmp=10;
	#pragma cetus private(i) 
	#pragma loop name foo35#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i)
	for (i=0; i<=9999; i+=1)
	{
		a[i]=tmp;
	}
	i=tmp;
	return ;
}

/*

 Scalar-output

*/
void foo36()
{
	int i;
	int tmp;
	#pragma cetus private(i, tmp) 
	#pragma loop name foo36#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, tmp)
	for (i=0; i<=9999; i+=1)
	{
		tmp=(a[i]+i);
	}
	return ;
}

void foo37()
{
	int i;
	int tmp;
	#pragma cetus private(i) 
	#pragma cetus lastprivate(tmp) 
	#pragma loop name foo37#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i) lastprivate(tmp)
	for (i=0; i<=9999; i+=1)
	{
		tmp=(a[i]+i);
	}
	i=tmp;
	return ;
}

/*

 Scalar privatization

*/
int a[10000];
int b[10000];
void foo38()
{
	int i;
	#pragma cetus private(i) 
	#pragma loop name foo38#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i)
	for (i=0; i<=9999; i+=1)
	{
		int tmp;
		tmp=(a[i]+i);
		b[i]=tmp;
	}
	return ;
}

/*

 Scalar true

*/
void foo39()
{
	int i;
	int tmp;
	#pragma cetus private(i, tmp) 
	#pragma loop name foo39#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, tmp)
	for (i=0; i<=9999; i+=1)
	{
		tmp=(a[i]+i);
		a[i]=tmp;
	}
	return ;
}

void foo40()
{
	int i;
	int tmp;
	#pragma cetus private(i) 
	#pragma cetus lastprivate(tmp) 
	#pragma loop name foo40#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i) lastprivate(tmp)
	for (i=0; i<=9999; i+=1)
	{
		tmp=(a[i]+i);
		a[i]=tmp;
	}
	i=tmp;
	return ;
}

/*

 Shared

*/
void foo41()
{
	int i;
	int x;
	int a[10000];
	#pragma cetus private(i) 
	#pragma loop name foo41#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i)
	for (i=0; i<=9999; i+=1)
	{
		a[i]=(a[i]+1);
	}
	return ;
}

/*

 Slow input

*/
void StressZero(real8 * newSxx, real8 * newSyy, real8 * newSzz, real8 * newTxy, real8 * newTxz, real8 * newTyz, const real8 * fun2j, const real8 * shearMod, real8 eosvmax, real8 stresscut, const int * zoneset, const real8 * vc, int length)
{
	int i;
	int index;
	/*
	This value 1.e-20 is used to prevent underflow. It is NOT a
	     cuttoff. DO NOT TOUCH THIS VALE.
	*/
	real8 stress2 = stresscut*1.0E-20;
	real8 nstres2 =  - stress2;
	#pragma cetus private(i, index) 
	#pragma loop name StressZero#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(22L*length)))) private(i, index)
	for (i=0; i<=(length-1); i+=1)
	{
		index=zoneset[i];
		if (((shearMod[index]==0.0)||(fun2j[i]<stresscut))||(vc[i]>=eosvmax))
		{
			newSxx[i]=0.0;
			newSyy[i]=0.0;
			newSzz[i]=0.0;
			newTxy[i]=0.0;
			newTxz[i]=0.0;
			newTyz[i]=0.0;
		}
		if ((newSxx[i]<stress2)&&(newSxx[i]>nstres2))
		{
			newSxx[i]=0.0;
		}
		if ((newSyy[i]<stress2)&&(newSyy[i]>nstres2))
		{
			newSyy[i]=0.0;
		}
		if ((newSzz[i]<stress2)&&(newSzz[i]>nstres2))
		{
			newSzz[i]=0.0;
		}
		if ((newTxy[i]<stress2)&&(newTxy[i]>nstres2))
		{
			newTxy[i]=0.0;
		}
		if ((newTxz[i]<stress2)&&(newTxz[i]>nstres2))
		{
			newTxz[i]=0.0;
		}
		if ((newTyz[i]<stress2)&&(newTyz[i]>nstres2))
		{
			newTyz[i]=0.0;
		}
	}
	return ;
}

/*

 Struct

*/
struct VectorXY
{
	double x;
	double y;
};

struct VectorXY v1[10000];
struct VectorXY v2[10000];
void applyVelocity()
{
	int in;
	#pragma cetus private(in) 
	#pragma loop name applyVelocity#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(in)
	for (in=0; in<=9999; in+=1)
	{
		v1[in].y=v2[3].y;
	}
	return ;
}

/*

 true dependence

*/
void foo42()
{
	int i;
	int a[1000];
	#pragma cetus private(i) 
	#pragma loop name foo42#0 
	for (i=0; i<=998; i+=1)
	{
		a[i+1]=(a[i]+1);
	}
	return ;
}

/*

 true dependence both levels

*/
void foo43()
{
	int n = 1000;
	int m = 1000;
	double a[n][m];
	int i;
	int j;
	#pragma cetus private(i, j) 
	#pragma loop name foo43#0 
	for (i=1; i<=(n-1); i+=1)
	{
		#pragma cetus private(j) 
		#pragma loop name foo43#0#0 
		for (j=1; j<=(m-1); j+=1)
		{
			a[i][j]=(a[i][j-1]+a[i-1][j]);
		}
	}
	return ;
}

/*

 true dependence scalar

*/
void foo44()
{
	int i;
	int a[1000];
	int temp;
	int t2;
	temp=1;
	t2=temp;
	#pragma cetus private(i, temp) 
	#pragma loop name foo44#0 
	for (i=0; i<=998; i+=1)
	{
		temp=a[i];
		a[i+1]+=(temp+1);
	}
	return ;
}

/*

 true I2

*/
int ii;
int jj;
int cc[1000][1000];
void foo45()
{
	#pragma cetus lastprivate(ii, jj) 
	#pragma loop name foo45#0 
	#pragma cetus parallel 
	#pragma omp parallel for lastprivate(ii, jj)
	for (ii=1; ii<=999; ii+=1)
	{
		#pragma cetus lastprivate(jj) 
		#pragma loop name foo45#0#0 
		for (jj=1; jj<=999; jj+=1)
		{
			cc[ii][jj]=(cc[ii][jj-1]+1);
		}
	}
	return ;
}

/*

 Uniform Indirect Indexed ArrayRefs

*/
double dnym1;
double dnzm1;
int grid_points[3];
double ux[24][(((24/2)*2)+1)][(((24/2)*2)+1)][5];
void initialize1()
{
	int i;
	int j;
	int k;
	int m;
	int ix;
	int iy;
	int iz;
	double xi;
	double eta;
	double zeta;
	double Pface[2][3][5];
	double Pxi;
	double Peta;
	double Pzeta;
	double temp[5];
	/* --------------------------------------------------------------------- */
	/* west face */
	/* --------------------------------------------------------------------- */
	i=0;
	xi=0.0;
	#pragma cetus private(eta, j, k, m, zeta) 
	#pragma loop name initialize1#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(4L*grid_points[2L]))+((19L*grid_points[1L])*grid_points[2L])))) private(eta, j, k, m, zeta)
	for (k=0; k<=(grid_points[2]-1); k+=1)
	{
		zeta=(((double)k)*dnzm1);
		#pragma cetus private(eta, j, m) 
		#pragma loop name initialize1#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(19L*grid_points[1L])))) private(eta, j, m)
		for (j=0; j<=(grid_points[1]-1); j+=1)
		{
			eta=(((double)j)*dnym1);
			/* exact_solution(xi,eta,zeta,temp); */
			#pragma cetus private(m) 
			#pragma loop name initialize1#0#0#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(m)
			*/
			for (m=0; m<=4; m+=1)
			{
				ux[k][j][i][m]=temp[m];
			}
		}
	}
	return ;
}
