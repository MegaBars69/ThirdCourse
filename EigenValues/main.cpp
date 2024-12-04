#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <cmath>
#include "initialize_matrix.hpp"
#include "algorithm.hpp"
#define EPSILON pow(10, -15)

using namespace std;

int main(int argc, char* argv[])
{
    if(argc < 5)
    {
        printf("Usage ./a.out n m eps s filename\n");
        return 1;
    }

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    double eps = atof(argv[3]);
    int s = atoi(argv[4]);

    if (n <= 0 || m < 0 || eps < 0 || s < 0 || s>4)
    {
        printf("Usage ./a.out n m eps s filename\n");
        return 1;
    }
    
    double* A = new double[n*n];
    int res_of_read = 0;
    if (s == 0)
    {
        res_of_read = ReadMatrixFromFile(argv[5], A, n);
        if (res_of_read > 0)
        {
            delete[] A;
            return res_of_read;
        }
        else
        {
            if(CheckMatrix(A, n) > 0)
            {
                printf("Matrix is'n simmetric\n");
                delete[] A;
                return 4;
            }
        }
        
    }
    else
    {
        FormulaMatrixInitialization(A, n, s);
    }

    PrintMatrix(A, n, m);

    double* U = new double[n];
    memset(U, 0, n* sizeof(double));

    double norm = Norm(A, U, n);
    double trace = Trace(A, n);
    cout<<trace<<endl;
    double mera = norm*EPSILON;   

    TriDiagonalize(A, U, n, mera);

    PrintMatrix(A, n, m); 
    cout<<Trace(A, n)<<endl;

    delete[] A;
    delete[] U;
    return 0;
}