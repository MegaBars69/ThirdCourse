#include <iostream>
#include "MatMul.hpp"
#include <cmath>
#include <ctime>

int main(void){
    double *A = new double[4];
    double *B = new double[4];
    double *C = new double[4];
    int n = 8;
    double *D = new double[n*n];
    double *E = new double[n*n];

    double *A1 = new double[n*n];
    double *A2 = new double[n*n];



    A[0] = 1; A[1] = 0; A[2] = 0; A[3] = 1;
    B[0] = 2; B[1] = 4; B[2] = 6; B[3] = 8;
    for (int i = 0; i < n*n; i++)
    {
        D[i] = i;
    }
    for (int i = 0; i < n; i++)
    {
       for (int j = 0; j < n; j++)
       {
            E[i*n + j] = (i != j ? 0:1);
       }
       
    }
    
    /*Print(A, 2); Print(B, 2);*/ Print(D, n); Print(E, n);
    
    
    /*Test(A,B,C, 2, MatMul1);
    Print(C, 2);

    Test(A,B,C, 2, MatMul2);
    Print(C, 2);

    Test(A,B,C, 2, MatMul3);
    Print(C, 2);*/
/*  Test(D,D,E, n, MatMul1);
    Print(E, n);
*/
    Test(D,E,A1, n, MatMul3);
    Print(A1, n);
    Test(D,E,A2, n, MatMul4);
    Print(A2, n);

    delete A; delete B; delete C, delete D, delete E;
    return 0;
}