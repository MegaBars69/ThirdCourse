#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cstring>
#include "initialize_matrix.hpp"
#include "algorithm.hpp"
#define EPSILON pow(10,-16)

using namespace std;

int main(int argc, char* argv[]) 
{
    Args* a;
    double* A, *B, **U, **ProductResult, **ZeroMatrix;
    bool *ZerosMatrix;
    int k;
    int n = 0, m = 0, r = 0, s = 0, p = 0, M = 0, ostatok = 0, Kk = 0;

    double t1 = 0.0, t2 =0, r1 = -1,r2 = -1;

    if (argc >=6)
    {
        n = atoi(argv[1]);
        m = atoi(argv[2]);
        r = atoi(argv[4]);
        p = atoi(argv[3]);
        s = atoi(argv[5]);

        if (p <= 0 || n <= 0 || m <=0 || r<0)
        {
            printf("Usage ./a.out n m p r s filename\n");
           
            return 4;
        }
        M = n/p;
        p = (n/m < p && n/m > 0 ? n/m : p);
        Kk = n/m;
        ostatok = (p > 0 ? n%p : 1);
        if (p <= 0 || n <= 0 || m <=0 || r<0)
        {
            printf("Usage ./a.out n m p r s filename\n");
           
            return 4;
        }
    }
    else
    {
        printf("Usage ./a.out n m p r s filename\n");
        return 1;
    }
    //std::string filename = argv[6];

    
    

    if(argc < 6)
    {
        printf("Usage ./a.out n m p r s filename\n");
        return 1;
    }

    if (p > 0)
    {
        a = new Args[p];
        if (n > 0)
        {
            A = new double[n*n];
            B = new double[n*n];
            ZerosMatrix = new bool[(Kk + 1)*(Kk + 1)];
            U = new double*[p];

            ProductResult = new double*[p];
            ZeroMatrix = new double*[p];
            for (int i = 0; i < p; i++)
            {
                U[i] = new double[(m+1)*(m+1)];
                ProductResult[i] = new double[m*m];
                ZeroMatrix[i] = new double[m*m];
            }
        }     
    }
    else
    {
        printf(" p >=0");
        return 4;
    }    
    
    memset(ZerosMatrix, 0, (Kk + 1)*(Kk + 1)*sizeof(bool));
    for (int i = 0; i < Kk +1 ; i++)
    {      
        ZerosMatrix[i*(Kk + 1) + i] = 1;        
    }
    
    //Обработка массива. Распределение памяти.
    for (k = 0; k < p; k++)
    {
        a[k].s = s;
        a[k].r = r;
        a[k].k = k;
        a[k].p = p;
        a[k].nomer_v_okne = k;
 
        a[k].n = n;
        a[k].m = m;
        a[k].M = M + (ostatok > 0 ? 1 : 0);
        ostatok--;

        if (s == 0 && argc == 7)
        {
            a[k].name = argv[6];
        }

        a[k].A = A;
        a[k].B = B;
        a[k].ZerosMatrix = ZerosMatrix;
        a[k].U = U[k];
        a[k].ProductResult = ProductResult[k];
        a[k].ZeroMatrix = ZeroMatrix[k];
    }

    //Запуск потоков.
    for (k = 1; k < p; k++)
    {
        if (pthread_create(&a[k].tid, nullptr, thread_func, a+k))
        {
            std::cerr << "Error creating thread " << k <<std::endl;
            for (int i = 0; i < p; i++)
            {
                delete[] U[i];
                delete[] ProductResult[i];
                delete[] ZeroMatrix[i];
            }
            delete[] A;
            delete[] B;
            delete[] ZerosMatrix;
            delete[] U;
            delete[] a;
            delete[] ProductResult;
            delete[] ZeroMatrix;
            return 1;
        }
    }
    a[0].tid = pthread_self();

    thread_func(a+0);

    //Прибиваем потоки молотком.
    for (k = 1; k < p; k++)
    {
        pthread_join(a[k].tid, nullptr);
    }
    if (a[0].res > 0)
    {
        printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m, atoi(argv[3]));
        for (int i = 0; i < p; i++)
        {
            delete[] U[i];
            delete[] ProductResult[i];
            delete[] ZeroMatrix[i];
        }
        delete[] A;
        delete[] B;
        delete[] ZerosMatrix;
        delete[] a;
        delete[] U;
        delete[] ProductResult;
        delete[] ZeroMatrix;

        return 2;
    }
    
    //Инициализируем матрицу для подсчета невязки
    
    if (s == 0) {
        string filename = argv[6];
        ReadMatrixFromFile(filename, A, n, m);
    } else {
        FormulaMatrixInitialization(A, n, m, s, 1, 0);
    }
    
    
    PrintMatrix(A,n,m,r,0,0,true,true,false);
     
    t1 = a[0].astr_time;
    double* Sum = new double[m*m];
    double* Column = new double[m];

    for (int i = 0; i < m; i++)
    {
        Column[i] = 0;
    }
    for (int i = 0; i < m*m; i++)
    {
        Sum[i] = 0;
    }
    
    

    //Подсчет невязки
    clock_t start2 = clock();
    if (n <= 11000 && fabs(a[0].res) <= EPSILON)
    {
        r1 = Discrepancy(A, B, Column, ProductResult[0], Sum, n, m);
        r2 = Discrepancy(B, A, Column, ProductResult[0], Sum, n, m);
        PrintMatrix(B, n, m, r, 0, 0, true, true, false);
    }
    else if(a[0].res > 0)
    {
        r1 = -1;
        r2 = -1;
    }
    else  
    {
        PrintMatrix(B, n, m, r, 0, 0, true, true, false);
        r1 = 0;
        r2 = 0; 
    }
    clock_t end2 = clock();

    t2 = static_cast<double>(end2 - start2) / CLOCKS_PER_SEC;
    
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m, atoi(argv[3]));
    
    delete[] A;
    delete[] B;
    delete[] ZerosMatrix;
    delete[] a;
    for (int i = 0; i < p; i++)
    {
        delete[] U[i];
        delete[] ProductResult[i];
        delete[] ZeroMatrix[i];
    }
    delete[] U;
    delete[] ProductResult;
    delete[] ZeroMatrix;
    delete[] Sum;
    delete[] Column;

    return 0;
}