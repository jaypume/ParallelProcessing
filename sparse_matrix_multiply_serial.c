#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define LEN 1000
#define NNZ 50000
typedef float dataType;
typedef struct{
	int row;
	int col;
	dataType	val;
}SNode;
typedef struct{
	int nRow,nCol,nNum;
	SNode data[NNZ];
}SM;
typedef struct{
	int nRow,nCol,nNum;
	SNode data[LEN*LEN];
}DM;
void sparse_matrix_produce(SM *p);
void sparse_matrix_multiply(SM *X,SM *Y,DM *Z);
void sparse_matrix_print(SM *p);



int main(int argc, char *argv){
	double time1,time2; 
	srand((unsigned int)time(NULL));
	SM *X=(SM *)malloc(sizeof(SM));
	SM *Y=(SM *)malloc(sizeof(SM));
	DM *Z=(DM *)malloc(sizeof(DM));
	sparse_matrix_produce(X);
	sparse_matrix_produce(Y);
	time1	=clock();
	sparse_matrix_multiply(X,Y,Z);			//Z=XY
	time2	=clock();
	//	sparse_matrix_print(X);
	//	sparse_matrix_print(Y);
	//	sparse_matrix_print(Z);
	printf("计算稀疏矩阵耗时：%6.3fs\n",(time2-time1)/CLOCKS_PER_SEC);
	
	
}

void sparse_matrix_produce(SM *p){
	int i,j,k;
	int	*A;
	if((A=(int *)malloc(sizeof(int)*LEN*LEN))==NULL){
		printf("error");
	}
	for(i=0;i<LEN;i++)
	for(j=0;j<LEN;j++){
		A[i*LEN+j]=0;
	}
	for(k=0;k<NNZ;k++){
		i	= rand()%LEN;
		j	= rand()%LEN;
		A[i*LEN+j]	= 1;
	}
	k=0;
	for(i=0;i<LEN;i++)
	for(j=0;j<LEN;j++){
		if(A[i*LEN+j]==1){
			(*p).data[k].row	=i;
			(*p).data[k].col	=j;
			(*p).data[k].val	=rand()%100;
			k++;
		}
	}
	(*p).nRow	= LEN;
	(*p).nCol	= LEN;
	(*p).nNum	= k;
}

void sparse_matrix_multiply(SM *X,SM *Y,DM *Z){
	int i,j,k,q;
	int index_of_node_on_X;
	int	count_of_node_on_Z;
	int	index_of_last_nonzero;
	int	number_of_nonzero_Y_row[(*Y).nRow];
	int	first_nonzero_on_Y_row[(*Y).nRow];
	int	adder[LEN];
	(*Z).nRow=(*X).nRow;
	(*Z).nCol=(*Y).nCol;
	for(i=0;i<(*Y).nRow;i++)number_of_nonzero_Y_row[i]=0;
	for(k=0;k<(*Y).nNum;k++){
		i=(*Y).data[k].row;
		number_of_nonzero_Y_row[i]++;
	}
	first_nonzero_on_Y_row[0]=0;
	for(i=1;i<(*Y).nRow;i++){
		first_nonzero_on_Y_row[i]=first_nonzero_on_Y_row[i-1]+number_of_nonzero_Y_row[i-1];
	}
	count_of_node_on_Z=0;
	index_of_node_on_X=0;
	for(i=0;i<(*X).nRow;i++){//X的每一行进行一次循环
		for(j=0;j<(*Y).nCol;j++){
			adder[j]=0;
		}
		while((*X).data[index_of_node_on_X].row==i){
			k=(*X).data[index_of_node_on_X].col;
			if(k<(*Y).nRow)
			index_of_last_nonzero=first_nonzero_on_Y_row[k+1];
			else 
			index_of_last_nonzero=(*Y).nNum;
			for(q=first_nonzero_on_Y_row[k];q<index_of_last_nonzero;q++){
				j=(*Y).data[q].col;
				adder[j]+=(*X).data[index_of_node_on_X].val*(*Y).data[q].val;
			}
			index_of_node_on_X++;
		}
		for(j=0;j<(*Y).nCol;j++){
			if(adder[j]){
				count_of_node_on_Z++;
				(*Z).data[count_of_node_on_Z].row	= i;
				(*Z).data[count_of_node_on_Z].col	= j;
				(*Z).data[count_of_node_on_Z].val	= adder[j];
			}
		}
	}
	(*Z).nNum=count_of_node_on_Z;
}

void sparse_matrix_print(SM *p){
	int i,j,k=0;
	printf("总元素的个数为%d,展示如下\n",(*p).nNum);
	for(i=0;i<LEN;i++){
		for(j=0;j<LEN;j++){
			if((i==(*p).data[k].row)&&(j==(*p).data[k].col)){
				printf("%0.0f\t",(*p).data[k].val);
				k++;
			}
			else {
				printf("0\t");
			}
		}
		printf("\n");
	}
}










