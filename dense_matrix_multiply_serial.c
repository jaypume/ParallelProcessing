#include <stdio.h>
#include "time.h"
#define LEN 1000
int num1[LEN][LEN];
int num2[LEN][LEN];
int result[LEN][LEN];
time_t c_start,c_end;
main ( int argc, char *argv[])
{
	int i,j,k;
	for(i=0;i <LEN;i++)
	{
		for(j=0;j <LEN;j++)
		{
			num1[i][j] = rand() % 100;
		}
	}
	for(i=0;i<LEN;i++)
	{
		for(j=0;j<LEN;j++)
		{
			num2[i][j] = rand() % 100;
		}
	}
	for(i=0;i<LEN;i++)
	{
		for(j=0;j<LEN;j++)
		{
			result[i][j] =0;
		}
	}

	c_start = clock();
	for(i=0;i<LEN;i++)
	{
		for(j=0;j<LEN;j++)
		{
			for(k=0;k<LEN;k++)
			{
				result[i][j]+=num1[i][k]*num2[k][j];
			}
		} 
	}
	c_end =clock();
	printf("串行方式执行时间: %fms\n",difftime(c_end,c_start));

}
