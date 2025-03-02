#include "mpi.h"
#include "initialize_matrix.hpp"
#include "algorithm.hpp"
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fenv.h>
#define EPSILON pow(10,-16)

using namespace std;

int main(int argc, char* argv[])
{
    int n = 0, m = 0, r = 0, s = 0, p = 0, proc_num = 0;
    double r1 = -1, r2 = -1, t1 = 0, t2 = 0;
    int er_l = 0, er_g = 0;
    int inv_res = 0;
    Args a;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &proc_num);

    
    if (argc != 5 && argc != 6) 
    {
        er_l = 1;
    } 
    else if (!(sscanf(argv[1], "%d", &n) == 1 && sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%d", &r) == 1 && sscanf(argv[4], "%d", &s) == 1)) 
    {
        er_l = 1;
    } 
    else if (n <= 0 || m <= 0 || r < 0 || s < 0 || s > 4 || n < m || p < 0 || (s==0 && argc==5)) 
    {
        er_l = 1;
    }

    MPI_Allreduce(&er_l, &er_g, 1, MPI_INT, MPI_MAX, comm);

    r = r > n ? n : r;

    if(er_g != 0) 
    {
        if (proc_num == 0)
        {
            printf("Errors with arguments\n");
            printf("Usage: mpirun -np p ./a.out n m r s filename\n");
        }
        MPI_Finalize();
        return 0;
    }
    int k = n/m, l = n%m;
    int max_rows = get_max_rows(n,m,p);

    a.n = n;
    a.m = m;
    a.active_procceses = a.p = (n/m + (l == 0 ? 0 : 1) >= p ? p : n/m + (l == 0 ? 0 : 1));
    a.k = proc_num;
    a.l = l;
    a.K = k;
    a.r = r;
    a.nomer_v_okne = proc_num;
    a.rows = get_rows(n, m, p, proc_num);
    a.max_rows = max_rows;
    a.last_line_isnt_fool = (l != 0 && (k%p) == proc_num);

    double* A = new double/*[n*m*(max_rows+p)]*/ [(max_rows + 1) * m * n + p * m * m];
    double* B = new double/*[n*m*(max_rows+p)]*/ [(max_rows + 1) * m * n + p * m * m];
    double* buf = new double[(2*n + m)*m];
    double* U = new double[(m+1)*(m+1)];
    bool* ZerosMatrix = new bool[(k+1)*(k+1)];
    double* ZeroMatrix = new double[m*m];
    double* ProductResult = new double[m*m];
    double* results = new double[n];
    double* tmp_row_matrix = new double[m * n + m * m];
    double* tmp_row_inverse = new double[m * n + m * m];

    if(A == nullptr || B == nullptr || buf == nullptr || U == nullptr || ZerosMatrix == nullptr || ProductResult == nullptr)
    {
        er_l = 1;
    }
    
    MPI_Allreduce(&er_l, &er_g, 1, MPI_INT, MPI_MAX, comm);

    if (er_g != 0)
    {
        if(proc_num == 0)
        {
            printf("Error with allocationg memmory\n");
        }
        
        if(proc_num == 0)
        {
            printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m, p);
        }
        
        delete[] A;
        delete[] B;
        delete[] buf;
        delete[] U;
        delete[] ZeroMatrix;
        delete[] ZerosMatrix;
        delete[] ProductResult;
        delete[] results;
        
        MPI_Finalize();
        return 0;
    }
    memset(A, 0, ((max_rows + 1) * m * n + p * m * m) * sizeof(double));
    memset(B, 0, ((max_rows + 1) * m * n + p * m * m) * sizeof(double));
    memset(U, 0, (m+1)*(m+1)*sizeof(double));
    memset(ProductResult, 0, m*m*sizeof(double));
    memset(results, 0, n*sizeof(double));
    memset(buf, 0, (2*n + m)*m*sizeof(double));
    memset(ZeroMatrix, 0, m*m*sizeof(double));

    memset(tmp_row_matrix, 0, (n * m + m * m) * sizeof(double));
    memset(tmp_row_inverse, 0, (n * m  + m * m) * sizeof(double));
    //memset(A, 0, n*m*(max_rows + p)*sizeof(double));

    for (int i = 0; i < k+1; i++)
    {
        for (int j = 0; j < k + 1; j++)
        {
            ZerosMatrix[i*(k+1) + j] = (i == j);
        }
        
    }

    if (s == 0)
    {
        er_l = ReadMatrixFromFile(A, n, m, p, proc_num, argv[5], buf, comm);
    }
    else
    {
        FormulaMatrixInitialization(A, n, m, p, proc_num, s);
    }
    
    MPI_Allreduce(&er_l, &er_g, 1, MPI_INT, MPI_MAX, comm);
    if (er_g != 0)
    {
        if(proc_num == 0)
        {
            if (er_g == 1)
            {
                printf("Error: opening file %s\n", argv[5]);
            }
            else if (er_g == 2)
            {
                printf("Error: Not enough elements from file %s\n",argv[5]);
            }
            else if (er_g ==3)
            {
                printf("Error: with reading element from file %s\n",argv[5]);
            }
            else if (er_g == 4)
            {
                printf("Error: Too many elements in file %s\n",argv[5]);
            }
            else if (er_g == 5)
            {
                printf("Error: Invalide format of file %s\n",argv[5]);
            }
            printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m, p);
        }
        
        delete[] A;
        delete[] B;
        delete[] buf;
        delete[] U;
        delete[] ZeroMatrix;
        delete[] ZerosMatrix;
        delete[] ProductResult;
        delete[] results;
        
        delete[] tmp_row_matrix;
        delete[] tmp_row_inverse;
        MPI_Finalize();

        return 0;
    }
    
    BuildE(B, n, m, p, proc_num);
    MPI_Barrier(comm);
    
    if(proc_num == 0)
    {
        cout<<"A"<<endl;
    }

    PrintMatrix(A, n, m, p, proc_num, r, buf, comm);
    MPI_Barrier(comm);

    a.A = A;
    a.B = B;
    a.buf = buf;
    a.ZerosMatrix = ZerosMatrix;
    a.ZeroMatrix = ZeroMatrix;
    a.ProductResult = ProductResult;
    a.U = U;
    a.results = results;
    a.norm = Norm(&a);
    a.norm *= EPSILON;

    /*if (proc_num == 0)
    {
        a.PrintInfo();
    }*/
    MPI_Barrier(comm);
    t1 = get_full_time();
    inv_res = InverseMatrixParallel(&a);
    t1 = get_full_time() - t1;
    if (inv_res == 0 && n <= 11000) 
    {
        if (s == 0)
        {
            er_l = ReadMatrixFromFile(A, n, m, p, proc_num, argv[5], buf, comm);
        }
        else
        {
            FormulaMatrixInitialization(A, n, m, p, proc_num, s);
        }

        MPI_Barrier(comm);
        t2 = get_full_time();
        r1 = calculate_discrepancy(A, B, n, m, tmp_row_matrix, tmp_row_inverse, proc_num, p, comm, tmp_row_matrix);
        r2 = calculate_discrepancy(B, A, n, m, tmp_row_matrix, tmp_row_inverse, proc_num, p, comm, tmp_row_matrix);
        t2 = get_full_time() - t2;
    }
    if(proc_num == 0)
    {
        cout<<"Inversed A"<<endl;
    }

    PrintMatrix(B, n, m, p, proc_num, r,buf, comm);
    //PrintMatrix(A, n, m, p, proc_num, r,buf, comm);

    if(proc_num == 0) 
    {
        printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m, p);
    }

    delete[] A;
    delete[] B;
    delete[] buf;
    delete[] U;
    delete[] ZeroMatrix;
    delete[] ZerosMatrix;
    delete[] ProductResult;
    delete[] results;

    delete[] tmp_row_matrix;
    delete[] tmp_row_inverse;

    MPI_Finalize();
    
    return 0;
}