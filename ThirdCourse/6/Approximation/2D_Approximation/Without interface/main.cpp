#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cstring>
#include "initialize_matrix.hpp"
#include "algorithm.hpp"
#include "thread_function.hpp"
#define EPSILON 1e-16

using namespace std;

int main(int argc, char* argv[]) 
{
    Args* aA;
    double* A, *B, *x,*u,*v,*r;
    int *I;
    int len_msr;
    int k, func_id;
    double a, b, c, d, eps;
    int nx, ny, max_it, p, its; 
    double t1 = 0.0, t2 =0, r1 = -1,r2 = -1, r3 = -1, r4 = -1;

    if (argc >= 11)
    {
        a = atoi(argv[1]);
        b = atoi(argv[2]);
        c = atoi(argv[4]);
        d = atoi(argv[3]);
        nx = atoi(argv[5]);
        ny = atoi(argv[6]);
        func_id = atoi(argv[7]);
        eps = atoi(argv[8]);
        max_it = atoi(argv[9]);
        p = atoi(argv[10]);

        if (nx <= 0 || ny<= 0 || func_id<0 || func_id > 8 || eps < 0 || p <= 0 || argc < 11)
        {
            printf("Usage ./a.out a b c d nx ny func_id epsilon max_iterations p");
           
            return 4;
        }
    }
    else
    {
        printf("Usage ./a.out a b c d nx ny func_id epsilon max_iterations p");
        return 1;
    }
    
    int N = (nx+1)*(ny+1);

    if (p > 0)
    {
        aA = new Args[p];
        N = (nx + 1)*(ny + 1);
        len_msr =  N + 1 + get_len_msr (nx, ny);
        A = new double[len_msr];
        I = new int[len_msr];
        B = new double[N];
        x = new double[N];
        r = new double[N];
        u = new double[N];
        v = new double[N];
    }
    else
    {
        printf(" p >=0");
        return 4;
    }    
    
    //Обработка массива. Распределение памяти.
    for (k = 0; k < p; k++)
    {
        aA[k].A = A;
        aA[k].B = B;
        aA[k].r = r;
        aA[k].u = u;
        aA[k].v = v;
        aA[k].x = x;
        aA[k].a = a;
        aA[k].b = b;
        aA[k].c = c;
        aA[k].d = d;
        aA[k].func_id = func_id;
        aA[k].maxit = max_it;
        aA[k].p = p;
        aA[k].k = k;
        aA[k].eps = eps;
        aA[k].nx = nx;
        aA[k].ny = ny;
        aA[k].N = (nx + 1)*(ny + 1);

    }

    //Запуск потоков.
    for (k = 1; k < p; k++)
    {
        if (pthread_create(&aA[k].tid, nullptr, thread_func, aA+k))
        {
            std::cerr << "Error creating thread " << k <<std::endl;
            delete[] A;
            delete[] B;
            delete[] x;
            delete[] u;
            delete[] v;
            delete[] r;
            delete[] I;
        
            return 1;
        }
    }
    its = aA->its;
    aA[0].tid = pthread_self();

    thread_func(aA+0);

    //Прибиваем потоки молотком.
    for (k = 1; k < p; k++)
    {
        pthread_join(aA[k].tid, nullptr);
    }

    if (aA[0].res > 0)
    {
        printf (
            "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\
            It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
            argv[0], 5, r1, r2, r3, r4, t1, t2, its, eps, k, nx, ny, p);
        
        delete[] A;
        delete[] B;
        delete[] x;
        delete[] u;
        delete[] v;
        delete[] r;
        delete[] I;
        return 2;
    }

    printf (
        "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\
        It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
        argv[0], 5, r1, r2, r3, r4, t1, t2, its, eps, k, nx, ny, p);
    
    delete[] A;
    delete[] B;
    delete[] x;
    delete[] u;
    delete[] v;
    delete[] r;
    delete[] I;
    return 0;
}