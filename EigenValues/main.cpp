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
    double res1 = 0, res2 = 0, t1 = 0, t2 = 0;
    double sum = 0, len = 0;
    int n = 1, its = 0;
    if(argc < 5)
    {
        printf("Usage ./a.out n m eps s filename\n");
        printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n", argv[0], res1, res2, its, its, t1, t2);
        return 1;
    }

    n = atoi(argv[1]);
    int m = atoi(argv[2]);
    double eps = atof(argv[3]);
    int s = atoi(argv[4]);

    if (n <= 0 || m < 0 || eps < 0 || s < 0 || s>4 || (s != 0 && argv[5] != nullptr))
    {
        printf("Usage ./a.out n m eps s filename\n");
        printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n", argv[0], res1, res2, its, its , t1, t2);

        return 1;
    }
    
    double* A = new double[n*n];
    int res_of_read = 0;
    if (s == 0)
    {
        res_of_read = ReadMatrixFromFile(argv[5], A, n);
        if (res_of_read > 0)
        {
            printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n", argv[0], res1, res2, its, its / n, t1, t2);
            
            delete[] A;
            return res_of_read;
        }
        else
        {
            if(CheckMatrix(A, n) > 0)
            {
                printf("Matrix is'n simmetric\n");
                printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n", argv[0], res1, res2, its, its / n, t1, t2);
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
    double Length = LengthOfMatrix(A, n);
    cout<<"trA = "<<trace<<endl;
    cout<<"||A|| = "<<Length<<endl;
    double mera = norm*EPSILON;   
    
    clock_t start1 = clock();

    TriDiagonalize(A, U, n, mera);

    clock_t end1 = clock();
    
    t1 = static_cast<double>(end1 - start1) / CLOCKS_PER_SEC;

    PrintMatrix(A, n, m); 
    trace = Trace(A, n);
    Length = LengthOfMatrix(A, n);
    cout<<"trA = "<<trace<<endl;
    cout<<"||A|| = "<<Length<<endl;

    clock_t start2 = clock();
    
    its = FindEigenValues(A, n, U, eps);

    clock_t end2 = clock();

    t2 = static_cast<double>(end2 - start2) / CLOCKS_PER_SEC;

    for (int i = 0; i < n; i++)
    {
        sum += U[i];
        len += U[i]*U[i];
    }
    for (int i = 0; i < min(n, m); i++)
    {
        cout<<U[i]<<" ";
    }
    
    cout<<endl<<"Sum = "<<sum<<endl;
    cout<<"Length = "<<sqrt(len)<<endl;
    res1 = fabs(sum - trace)/norm;
    res2 = fabs(Length - sqrt(len))/norm;
    cout<<endl;


    printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n", argv[0], res1, res2, its, its / n, t1, t2);

    delete[] A;
    delete[] U;
    return 0;
}