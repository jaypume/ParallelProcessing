#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

#define LEN      1000
#define TAG_A     1
#define TAG_B     2
#define TAG_C     3
int      get_rank(int row, int col, int sp);
void     produce_matrix_AB();
void     scatter_AB_to_rankx();
void     recive_AB_from_rank0();
void     matrix_AB_alignment();
void     multiply_AB_by_cannon();
void     recive_C_from_rankx();
void     send_C_to_rank0();


float     **A,**B,**C;
float     *a,*b,*c,*buf_a,*buf_b;
int       d,d2,p,sp;//p is the number of processor
int       myRank,myRow,myCol;
MPI_Status     status;

int main(int argc, char *argv[])
{
	int      i,nameLen;
	char     nameProcessor[MPI_MAX_PROCESSOR_NAME];
	double      runtime;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	MPI_Get_processor_name(nameProcessor,&nameLen);
	sp = sqrt(p);
	if((sp*sp!=p)&&(myRank==0)){
		printf("处理器数目必须为平方数！\n");
		MPI_Finalize();
		exit(1);
	}
	myRow      = myRank/sp;
	myCol      = myRank%sp;
	d     = LEN/sp;
	d2    = d * d;
	a     = (float *)malloc(d2 * sizeof(float));
	b     = (float *)malloc(d2 * sizeof(float));
	c     = (float *)malloc(d2 * sizeof(float));
	buf_a     = (float *)malloc(d2 * sizeof(float));
	buf_b     = (float *)malloc(d2 * sizeof(float));
	for(i=0;i<d2;i++)c[i] = 0.0;
	if(myRank == 0){
		produce_matrix_AB();
	}
	runtime     = - MPI_Wtime();
	if(myRank == 0){
		scatter_AB_to_rankx();
		printf("scatter AB over!\n");
	}
	else{
		recive_AB_from_rank0();
	}
	printf("Rank %d has received AB!\n",myRank);
	matrix_AB_alignment();
	multiply_AB_by_cannon();
	printf("cannon multiply AB on rank %d is over\n",myRank);
	if(myRank == 0){
		recive_C_from_rankx();
		printf("Rank0 have recived C from ALL other ranks!\n");
	}
	else{
		send_C_to_rank0();
		printf("have sent c to rank0 from %d\n",myRank);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	runtime     += MPI_Wtime();
	printf("Process %d on %s 's runtime is %lf ms\n",myRank,nameProcessor,runtime*1000);
	MPI_Finalize();
	return 0;          
}
int          get_rank(int row, int col, int sp){
	return ((row+sp)%sp)*sp + (col+sp)%sp;
}

void     produce_matrix_AB(){
	int i,j;
	A     = (float **)malloc(LEN * sizeof(float*));
	B     = (float **)malloc(LEN * sizeof(float*));
	C     = (float **)malloc(LEN * sizeof(float*));
	for(i=0;i<LEN;i++){
		A[i]     = (float *)malloc(LEN * sizeof(float));
		B[i]     = (float *)malloc(LEN * sizeof(float));
		C[i]     = (float *)malloc(LEN * sizeof(float));
	}
	srand((unsigned int)time(NULL));
	for(i=0;i<LEN;i++)
	for(j=0;j<LEN;j++){
		A[i][j]     = rand();
		B[i][j]     = rand();
		C[i][j]     = 0.0;
	}
}

void     scatter_AB_to_rankx(){
	int     i,j,k,dest;
	int     minRow,maxRow,minCol,maxCol;
	for(dest=0;dest<p;dest++){
		minCol     = (dest % sp) * d;
		maxCol     = (dest % sp) * d + d - 1;
		minRow     = (int)(dest / sp) * d;
		maxRow     = (int)(dest / sp) * d + d - 1;
		printf("rank%d分配到的矩阵坐标范围,minCol:%d,maxCol:%d,minRow:%d,maxRow:%d\n",dest,minCol,maxCol,minRow,maxRow);
		k     = 0;
		for(i=minRow;i<maxRow;i++)
		for(j=minCol;j<maxCol;j++){
			buf_a[k]     = A[i][j];
			buf_b[k]     = B[i][j];
			k++;
		}
		if(dest==0){
			memcpy(a, buf_a, d2 * sizeof(float));
			memcpy(b, buf_b, d2 * sizeof(float));
		}
		else{
			MPI_Send(buf_a, d2, MPI_FLOAT, dest, TAG_A, MPI_COMM_WORLD);
			MPI_Send(buf_b, d2, MPI_FLOAT, dest, TAG_B, MPI_COMM_WORLD); 
		}
	}
}

void     recive_AB_from_rank0(){
	MPI_Recv(a, d2, MPI_FLOAT, 0, TAG_A, MPI_COMM_WORLD, &status);
	MPI_Recv(b, d2, MPI_FLOAT, 0, TAG_B, MPI_COMM_WORLD, &status);     
}

void     matrix_AB_alignment(){
	MPI_Sendrecv(a, d2, MPI_FLOAT, get_rank(myRow, myCol-myRow, sp), 1,
	buf_a, d2, MPI_FLOAT, get_rank(myRow, myCol+myRow, sp), 1, MPI_COMM_WORLD, &status);
	memcpy(a, buf_a, d2*sizeof(float));
	MPI_Sendrecv(b, d2, MPI_FLOAT, get_rank(myRow-myCol, myCol, sp), 1,
	buf_b, d2, MPI_FLOAT, get_rank(myRow+myCol, myCol, sp), 1, MPI_COMM_WORLD, &status);
	memcpy(b, buf_b, d2*sizeof(float));     
}

void     multiply_AB_by_cannon(){
	int i,j,k,l,bufSize;
	float *bufTemp;
	for(l=0;l<sp;l++){
		for(i=0;i<d;i++)
		for(j=0;j<d;j++)
		for(k=0;k<d;k++){
			c[i*d+j] += a[i*d+k]*b[k*d+j];
		}
		printf("calculation of c on process %d is over\n",myRank);	
		MPI_Pack_size(d2, MPI_FLOAT, MPI_COMM_WORLD, &bufSize);
		bufTemp	= (float *)malloc(bufSize + 2*MPI_BSEND_OVERHEAD);
		if(!bufTemp){
			printf("申请缓存失败");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		MPI_Buffer_attach(bufTemp, bufSize+MPI_BSEND_OVERHEAD);
		MPI_Bsend(a, d2, MPI_FLOAT, get_rank(myRow, myCol-1, sp), TAG_A, MPI_COMM_WORLD);
		MPI_Recv(a, d2, MPI_FLOAT, get_rank(myRow, myCol+1, sp), TAG_A, MPI_COMM_WORLD, &status);
		MPI_Buffer_detach(bufTemp, &bufSize);
		MPI_Buffer_attach(bufTemp, bufSize+MPI_BSEND_OVERHEAD);
		MPI_Bsend(b, d2, MPI_FLOAT, get_rank(myRow-1, myCol, sp), TAG_B, MPI_COMM_WORLD);
		MPI_Recv(b, d2, MPI_FLOAT, get_rank(myRow+1, myCol, sp), TAG_B, MPI_COMM_WORLD, &status);
		MPI_Buffer_detach(bufTemp, &bufSize);
	}

}

void     recive_C_from_rankx(){
	int     i,j,x,y,dest;
	int     minRow,maxRow,minCol,maxCol;
	for(i=0;i<d;i++)
	for(j=0;j<d;j++)
	C[i][j]=c[i*d+j];
	printf("rank%d's matrix c is received\n",myRank);
	for(dest=1;dest<p;dest++){
		MPI_Recv(c, d2, MPI_FLOAT, dest, TAG_C, MPI_COMM_WORLD, &status);
		minCol     = (dest % sp) * d;
		maxCol     = (dest % sp) * d + d - 1;
		minRow     = (int)(dest / sp) * d;
		maxRow     = (int)(dest / sp) * d + d - 1;
		x     = 0;
		for(i=minRow;i<maxRow;i++){
			y     = 0;
			for(j=minCol;j<maxCol;j++){
				C[i][j]=c[x*d+y];
				y++;
			}
			x++;
		}
		printf("rank%d's matrix c is received\n",dest);
	}     
}

void     send_C_to_rank0(){
	MPI_Send(c, d2, MPI_FLOAT, 0, TAG_C, MPI_COMM_WORLD);
}







