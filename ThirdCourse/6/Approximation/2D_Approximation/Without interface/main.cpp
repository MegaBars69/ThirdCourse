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
#include <fenv.h>
#define EPSILON 1e-16

using namespace std;

int main(int argc, char* argv[]) 
{ 
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    Args* aA;
    double* A, *B, *x,*u,*v,*r;
    int *I;
    int len_msr;
    int k, func_id;
    double a, b, c, d, eps;
    int nx, ny, max_it, p, its; 
    double t1 = 0.0, t2 =0, r1 = -1,r2 = -1, r3 = -1, r4 = -1;

    if (argc != 11 || sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 || sscanf(argv[3], "%lf", &c) != 1 || sscanf(argv[4], "%lf", &d) != 1 || sscanf(argv[5], "%d", &nx) != 1 || sscanf(argv[6], "%d", &ny) != 1|| sscanf(argv[7], "%d", &func_id) != 1 || sscanf(argv[8], "%lf", &eps) != 1 || sscanf(argv[9], "%d", &max_it) != 1 || sscanf(argv[10], "%d", &p) != 1) 
    {
        printf("Usage1 ./a.out a b c d nx ny func_id epsilon max_iterations p\n");
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
        init_reduce_sum (p);
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
        aA[k].I = I;
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
        aA[k].len_msr = len_msr;
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
            delete[] aA;
            free_reduce_sum();
        
            return 1;
        }
    }
    
    aA[0].tid = pthread_self();

    thread_func(aA+0);

    //Прибиваем потоки молотком.
    for (k = 1; k < p; k++)
    {
        pthread_join(aA[k].tid, nullptr);
    }

    its = aA->its;
    r1 = aA->r1;
    r2 = aA->r2;
    r3 = aA->r3;
    r4 = aA->r4;
    t1 = aA->t1;
    t2 = aA->t2;

    printf (
        "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\
        It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
        argv[0], 5, r1, r2, r3, r4, t1, t2, its, eps, func_id, nx, ny, p);
    
    delete[] A;
    delete[] B;
    delete[] x;
    delete[] u;
    delete[] v;
    delete[] r;
    delete[] I;
    delete[] aA;
    free_reduce_sum();
    return 0;
}