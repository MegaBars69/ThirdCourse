#include <iostream>
#include <cmath>
#include "algorithm.hpp"
#include "initialize_matrix.hpp"
#include "string.h"
#include <cstring>
#include "thread_function.hpp"

double func0(double, double)
{
    return 1.0;
}

double func1(double x, double) 
{
    return x;
}

double func2(double, double y) 
{
    return y;
}

double func3(double x, double y) 
{
    return x + y;
}

double func4(double x, double y) 
{
    return sqrt(x * x + y * y);
}

double func5(double x, double y) 
{
    return x * x + y * y;
}

double func6(double x, double y) 
{
    return exp(x * x - y * y);
}

double func7(double x, double y) 
{
    return 1.0 / (25.0 * (x * x + y * y) + 1);
}

typedef double (*FuncPtr)(double, double);

FuncPtr retrieve_function(int index) {
    switch (index) {
        case 0: return func0;
        case 1: return func1;
        case 2: return func2;
        case 3: return func3;
        case 4: return func4;
        case 5: return func5;
        case 6: return func6;
        case 7: return func7;
        default: return nullptr;
    }
}
void PrintMatrix(double*A, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout<<A[i]<<" ";
    }
    std::cout<<"\n";
}

void PrintMatrixInt(int*A, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout<<A[i]<<" ";
    }
    std::cout<<"\n";
}

void* thread_func(void *arg)
{  
    Args *  aa = (Args *) arg; 
    bool* working = aa->working;
    bool* threads_quiting = aa->quiting_app;
    pthread_mutex_t* p_mutex = aa->p_mutex;
	pthread_cond_t* p_cond = aa->p_cond;
    while (1)
    {
        pthread_mutex_lock(p_mutex);
        // Ожидание, пока есть работа или требуется выход
        while (!(*working) && !(*threads_quiting)) {
            pthread_cond_wait(p_cond, p_mutex);
        }
        // Проверка на необходимость завершения
        if (*threads_quiting) {
            pthread_mutex_unlock(p_mutex);
            return nullptr;
        }
        pthread_mutex_unlock(p_mutex);
    
        int nx = aa->nx;
        int ny = aa->ny;
        double a = aa->a;
        double b = aa->b;
        double c = aa->c;
        double d = aa->d;
        int func_id = aa->func_id;
        int p = aa->p;
        double eps = aa->eps;
        int maxit = aa->maxit;
        int k = aa->k;
        int its;
        double* A = aa->A;
        double* B = aa->B;
        double*x =aa->x;
        double*r =aa->r;
        double*u =aa->u;
        double*v =aa->v;
        int *I = aa->I;
        cpu_set_t cpu;
        CPU_ZERO(&cpu);
        pthread_t tid = aa->tid;
        CPU_SET(p - 1 - (k % (p)), &cpu);
        pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
        FuncPtr f = retrieve_function(func_id);

        pthread_mutex_t* p_mutex = aa->p_mutex;
        //pthread_cond_t* p_cond = aa->p_cond;

        double hx = (b - a)/nx, hy = (d - c)/ny;
        aa->hx = hx, aa->hy = hy;
        int N = aa->N;
        double t1 = -1, t2 = -1, r1 = -1, r2 = -1, r3 = -1, r4 = -1, err = 0;

        int point = aa->point;
        double norm = aa->norm;

        if(k==0)
        {
            pthread_mutex_lock(p_mutex);
            *working = true;
            pthread_mutex_unlock(p_mutex);
        }

        memset(x, 0, N*sizeof(double));
        t1 = get_full_time();

        if(k == 0)
        {
            fill_I(nx, ny, I);
        }
        reduce_sum(p);

        err = fill_IA(nx, ny, hx, hy, I, A, p, k);
        reduce_sum(p, &err, 1);
        if(fabs(err) > 0)
        {
            return nullptr;
        }
        fill_B (N, nx, ny, hx, hy, B, a, c, p, k, f, point, norm);
        reduce_sum(p);

        its = minimal_error_msr_matrix_full(N, A, I, B, x, r, u, v, eps, maxit, 100, p, k); 
        t1 = get_full_time() - t1;

        reduce_sum(p);
        t2 = get_full_time();
        ResidualCalculation(&r1, &r2, &r3, &r4, x, a, c, hx, hy, nx, ny, p, k, f);
        t2 = get_full_time() - t2;
        reduce_sum(p);

        aa->t1 = t1;
        aa->t2 = t2;
        aa->r1 = r1;
        aa->r2 = r2;
        aa->r3 = r3;
        aa->r4 = r4;
        aa->its = its;
        if(k == 0)
        {
            printf ("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\
                It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
                "./a.out", 5, r1, r2, r3, r4, t1, t2, its, eps, func_id, nx, ny, p);

        }   

    
    
        *working = false;
    } 
    return nullptr;
}