#include "mpi.h"
#include "initialize_matrix.hpp"
#include "algorithm.hpp"
#include <stdio.h>
#include <iostream>
#include <fenv.h>
#include <new>

using namespace std;

int main(int argc, char* argv[])
{
    int n = 0, m = 0, r = 0, s = 0, p = 0, proc_num = 0;
    double r1 = -1, r2 = -1, t1 = 0, t2 = 0;
    MPI_Comm comm = MPI_COMM_WORLD;
    int er_l = 0, er_g = 0;
    MPI_Init(&argc, &argv);

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    er_g = er_l;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &proc_num);


    if (argc != 5 && argc != 6) {
        er_l = 1;
        if (proc_num == 0) {
            printf("Wrong amount of arguments\n");
        }
    } else if (!(sscanf(argv[1], "%d", &n) == 1 && sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%d", &r) == 1 && sscanf(argv[4], "%d", &s) == 1)) 
    {
        er_l = 1;
        if (proc_num == 0) {
            printf("Wrong type of arguments\n");
        }
    } else if (n <= 0 || m <= 0 || r < 0 || s < 0 || s > 4 || (s==0 && argc==5) || n < m || p <= 0) {
        er_l = 1;
        if (proc_num == 0) {
            printf("Wrong arguments\n");
        }
    }
    MPI_Allreduce(&er_l, &er_g, 1, MPI_INT, MPI_MAX, comm);
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

    if(er_g == 0) {
        if (proc_num == 0) {
            printf("Inverse matrix :\n");
        }
        //print_matrix(inverse, n, m, proc_num, p, r, tmp_row_matrix, comm);
    }

    int max_rows = n/p + 1;
    double* A = new(std::nothrow) double[n*max_rows];
    double* B = new(std::nothrow) double[n*max_rows];
    double* buf = new(std::nothrow) double[n*m];
    
    memset(A, 0, n*max_rows * sizeof(double));
    memset(B, 0, n*max_rows * sizeof(double));
    memset(buf, 0, n*m * sizeof(double));

    if (s == 0)
    {
        er_l = ReadMatrixFromFile(A, n, m, p, proc_num, argv[5], buf, comm);

    }
    else
    {
        FormulaMatrixInitialization(A, n, m, p, proc_num, s);
    }
    MPI_Barrier(comm);

    BuildE(B, n, m, p, proc_num);
    MPI_Barrier(comm);

    PrintMatrix(A, n, m, p, proc_num, r,buf, comm);
    MPI_Barrier(comm);

    PrintMatrix(B, n, m, p, proc_num, r,buf, comm);

    if(proc_num == 0) 
    {
        printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m, p);
    }

    delete[] A;
    delete[] B;
    delete[] buf;
    
    MPI_Barrier(comm);
    MPI_Finalize();

    return 0;
}