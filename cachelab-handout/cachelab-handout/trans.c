/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"
#define BLOCKSIZE 8

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    //temp var
    int a1,a2,a3,a4,a5,a6,a7,a8;
    if(M==32)
    {
        for(int k = 0;k < M ;k+=8)
        {
            for(int l =0 ;l < N; l++)
            {
                a1=A[l][k];
                a2=A[l][k+1];
                a3=A[l][k+2];
                a4=A[l][k+3];
                a5=A[l][k+4];
                a6=A[l][k+5];
                a7=A[l][k+6];
                a8=A[l][k+7];
                B[k][l]=a1;
                B[k+1][l]=a2;
                B[k+2][l]=a3;
                B[k+3][l]=a4;
                B[k+4][l]=a5;
                B[k+5][l]=a6;
                B[k+6][l]=a7;
                B[k+7][l]=a8;  
            }
        }
    }
    else if(M==64)
    {
        for(int i=0;i<64;i+=8)
        {
            for(int j=0;j<64;j+=8)
            {
            //将A矩阵前4行拷贝到B矩阵，同时完成B矩阵左上角4*4区域的转制操作
            //total miss: 8 times
                for(int k=j; k<j+4; k++)
                {
               		//first: cache line 1
               		//then: cache line 6, 7, 8
                    a1=A[k][i];                   
                    a2=A[k][i+1];
                    a3=A[k][i+2];
                    a4=A[k][i+3];
                    a5=A[k][i+4];
                    a6=A[k][i+5];
                    a7=A[k][i+6];
                    a8=A[k][i+7];
               		//line 2    
                    B[i][k]=a1;
                    B[i][k+4]=a5;
                    //line 3
                    B[i+1][k]=a2;   
                    B[i+1][k+4]=a6;  
                    //line 4               
                    B[i+2][k]=a3;  
                    B[i+2][k+4]=a7;    
                    //line 5            
                    B[i+3][k]=a4;                                      
                    B[i+3][k+4]=a8;                         
                }
                //完成B矩阵左下角4*4以及右上角4*4区域的转制
                //total miss: 8 times
                for(int k=i; k<i+4; k++)
                {
                	//hit cache line 2,3,4,5
         			a1=B[k][j+4];
                    a2=B[k][j+5];
                    a3=B[k][j+6];
                    a4=B[k][j+7];
                    
                    /*evivts line 1: A[0][0] -- A[0][7]
                    *	     line 6: A[1][0] -- A[1][7]
                    *        line 7: A[2][0] -- A[2][7]
                    *        line 8: A[3][0] -- A[3][7]
                    *to A[4][0] -- A[4][7]
                    *   A[5][0] -- A[5][7]
                    *   A[6][0] -- A[6][7]
                    *   A[7][0] -- A[7][7]
                    * miss: 4 times
                    */
                    a5=A[j+4][k];
                    a6=A[j+5][k];
                    a7=A[j+6][k];
                    a8=A[j+7][k];
                    
                    //hit cache line 2,3,4,5
                    //完成B矩阵右上角4*4区域的转制
                    B[k][j+4]=a5;
                    B[k][j+5]=a6;
                    B[k][j+6]=a7;
                    B[k][j+7]=a8;
                    
                    /*evivts line 1: A[4][0] -- A[4][7]
                    *	     line 6: A[5][0] -- A[5][7]
                    *        line 7: A[6][0] -- A[6][7]
                    *        line 8: A[7][0] -- A[7][7]
                    *to B[4][0] -- A[4][7]
                    *   B[5][0] -- A[5][7]
                    *   B[6][0] -- A[6][7]
                    *   B[7][0] -- A[7][7]
                    * miss: 4 times
                    */
                    //完成B矩阵左下角4*4区域的转制
                    B[k+4][j]=a1;
                    B[k+4][j+1]=a2;
                    B[k+4][j+2]=a3;
                    B[k+4][j+3]=a4;
                   
					
                }
                //B矩阵右下角4*4区域转制
                //total miss: 4 times
                for(int k=i+4; k<i+8; k++)
                {
                	/*evivts line 2,3,4,5
                	*from B[0][0] -- B[0][7]
                	*	  B[1][0] -- B[1][7]
                	*
                	*TO A[4][0] -- A[4][7]
                	*	A[5][0] -- A[5][7]
                	*
                	*miss: 4 times
                	*/
                    a1=A[j+4][k];
                    a2=A[j+5][k];
                    a3=A[j+6][k];
                    a4=A[j+7][k];

					//hit cache line: 1,6,7,8
                    B[k][j+4]=a1;
                    B[k][j+5]=a2;
                    B[k][j+6]=a3;
                    B[k][j+7]=a4;
                }
            }
        }
    }
    else
    {
        for(int k = 0;k < M ;k+=17)
        {
            for(int l =0 ;l < N; l+=17)
            {
                for(int i=l ; i<l+17 && (i<N);i++)
                {
                    for(int j=k;j<k+17 &&(j<M);j++)
                    {
                       B[j][i]=A[i][j];
                    }
                }
            }
        }
    }
}

/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;
    for (i = 0; i < N; i++) 
    {
        for (j = 0; j < M; j++) 
        {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }
}

char trans_desc2[] = "column-wise scan transpose";
void trans2(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;
    for (i = 0; i < M; i++) 
    {
        for (j = 0; j < N; j++) 
        {
            tmp = A[j][i];
            B[i][j] = tmp;
        }
    }
}

char trans_desc3[] = "using a zig-zag access pattern";
void trans3(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;
    for (i = 0; i < N; i++) 
    {
        for (j = 0; j < M; j++) 
        {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }
}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc);
    registerTransFunction(trans2, trans_desc2); 
    registerTransFunction(trans3, trans_desc3); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

