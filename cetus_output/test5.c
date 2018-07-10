/*
bitonic sort

*/
/* #include <stdio.h> */
void merge_down(int * arr, int n)
{
	int step = n/2, i, j, k, temp;
	while (step>0)
	{
		#pragma loop name merge_down#0 
		for (i=0; i<n; i+=(step*2))
		{
			#pragma loop name merge_down#0#0 
			for (((j=i), (k=0)); k<step; ((j ++ ), (k ++ )))
			{
				if (arr[j]<arr[j+step])
				{
					/* swap */
					temp=arr[j];
					arr[j]=arr[j+step];
					arr[j+step]=temp;
				}
			}
		}
		step/=2;
	}
	return ;
}

void mergeup(int * arr, int n)
{
	int step = n/2, i, j, k, temp;
	while (step>0)
	{
		#pragma loop name mergeup#0 
		for (i=0; i<n; i+=(step*2))
		{
			#pragma loop name mergeup#0#0 
			for (((j=i), (k=0)); k<step; ((j ++ ), (k ++ )))
			{
				if (arr[j]>arr[j+step])
				{
					/* swap */
					temp=arr[j];
					arr[j]=arr[j+step];
					arr[j+step]=temp;
				}
			}
		}
		step/=2;
	}
	return ;
}

/*

        void merge_up(intarr, int n) {      	
	          int step=n/2,i,j,k,temp;
	          while (step > 0) {
		            for (i=0; i < n; i+=step*2) {
			              for (j=i,k=0; k < step; j++,k++) {
				            	if (arr[j] > arr[j+step]) {
					              	swap
					              	temp = arr[j];
					              	arr[j]=arr[j+step];
					              	arr[j+step]=temp;
				            	}
			              }
		            }
		            step /= 2;
	          }
        }
       
*/
/*

        void printArray(intarr, int n) {
	          int i;
	
	          printf("[%d",arr[0]);
	          for (i=1; i < n;i++) {
		            printf(",%d",arr[i]);
	          }
	          printf("]\n");
        }

*/
int main()
{
	int n, * arr, i, s;
int arr1[16] = {3, 5, 8, 9, 10, 12, 14, 20, 95, 90, 60, 40, 35, 23, 18, 0};
	int _ret_val_0;
	arr=arr1;
	n=16;
	/* print array before  */
	/* printArray(arr,n); */
	/* do merges */
	#pragma loop name main#0 
	for (s=2; s<=n; s*=2)
	{
		#pragma loop name main#0#0 
		for (i=0; i<n; i+=(s*2))
		{
			merge_up(arr+i, s);
			merge_down((arr+i)+s, s);
		}
	}
	/* printArray(arr,n); */
	_ret_val_0=0;
	return _ret_val_0;
}
