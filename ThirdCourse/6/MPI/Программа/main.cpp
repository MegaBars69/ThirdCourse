#include "mpi.h"
#include "initialize_matrix.hpp"
#include "algorithm.hpp"
#include <stdio.h>
#include <iostream>
#include <fenv.h>

using namespace std;

int main(int argc, char* argv[])
{
    int n = 0, m = 0, r = 0, s = 0, p = 0, proc_num = 0;
    double r1 = -1, r2 = -1, t1 = 0, t2 = 0;
    int er_l = 0, er_g = 0;
    
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
    
    int max_rows = get_max_rows(n,m,p);
    double* A = new double[n*m*(max_rows+1)];
    double* B = new double[n*m*(max_rows+1)];
    double* buf = new double[(n + m)*m];

    if(A == nullptr || B == nullptr || buf == nullptr)
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
        
        MPI_Finalize();
        return 0;
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
        
        MPI_Finalize();

        return 0;
    }
    
    BuildE(B, n, m, p, proc_num);
    MPI_Barrier(comm);
    
    if(proc_num == 0)
    {
        cout<<"A"<<endl;
    }

    PrintMatrix(A, n, m, p, proc_num, r,buf, comm);
    MPI_Barrier(comm);

    if(proc_num == 0)
    {
        cout<<"Inversed A"<<endl;
    }

    PrintMatrix(B, n, m, p, proc_num, r,buf, comm);
    
    if(proc_num == 0) 
    {
        printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m, p);
    }

    delete[] A;
    delete[] B;
    delete[] buf;
    
    MPI_Finalize();
    
    return 0;
}