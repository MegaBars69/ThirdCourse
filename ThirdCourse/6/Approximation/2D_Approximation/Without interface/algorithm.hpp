#ifndef HEADER1  
#define HEADER1
#include <iostream>
#include <sys/sysinfo.h>
#include "initialize_matrix.hpp"

void reduce_sum(int p,double * a = nullptr, int n = 0)
{
    static pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static double* r = nullptr;
    int i;
    if (p <= 1)
    {
        return;
    }

    pthread_mutex_lock(&my_mutex);

    if (r==nullptr)
    {
        r = a;
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            r[i] += a[i];
        }
    }

    t_in++;

    if (t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else
    {
        while (t_in < p)
        {
            pthread_cond_wait(&c_in, &my_mutex);
        }
    }

    if(r != a)
    {
        for (i = 0; i < n; i++)
        {
            a[i] = r[i];
        }
    }
    t_out++;
    if (t_out >= p)
    {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    }
    else
    {
        while (t_out < p)
        {
            pthread_cond_wait(&c_out, &my_mutex);
        }
        
    }
    pthread_mutex_unlock(&my_mutex);
}

void reduce_sum_int(int p,int * a = nullptr, int n = 0)
{
    static pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static int* r = nullptr;
    int i;
    if (p <= 1)
    {
        return;
    }

    pthread_mutex_lock(&my_mutex);

    if (r==nullptr)
    {
        r = a;
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            r[i] += a[i];
        }
    }

    t_in++;

    if (t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else
    {
        while (t_in < p)
        {
            pthread_cond_wait(&c_in, &my_mutex);
        }
    }

    if(r != a)
    {
        for (i = 0; i < n; i++)
        {
            a[i] = r[i];
        }
    }
    t_out++;
    if (t_out >= p)
    {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    }
    else
    {
        while (t_out < p)
        {
            pthread_cond_wait(&c_out, &my_mutex);
        }
        
    }
    pthread_mutex_unlock(&my_mutex);
}


void reduce_sum(int p,double * a = nullptr, int n = 0)
{
    static pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static double* r = nullptr;
    int i;
    if (p <= 1)
    {
        return;
    }

    pthread_mutex_lock(&my_mutex);

    if (r==nullptr)
    {
        r = a;
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            r[i] += a[i];
        }
    }

    t_in++;

    if (t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else
    {
        while (t_in < p)
        {
            pthread_cond_wait(&c_in, &my_mutex);
        }
    }

    if(r != a)
    {
        for (i = 0; i < n; i++)
        {
            a[i] = r[i];
        }
    }
    t_out++;
    if (t_out >= p)
    {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    }
    else
    {
        while (t_out < p)
        {
            pthread_cond_wait(&c_out, &my_mutex);
        }
        
    }
    pthread_mutex_unlock(&my_mutex);
}

void reduce_sum_int(int p,int * a = nullptr, int n = 0)
{
    static pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static int* r = nullptr;
    int i;
    if (p <= 1)
    {
        return;
    }

    pthread_mutex_lock(&my_mutex);

    if (r==nullptr)
    {
        r = a;
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            r[i] += a[i];
        }
    }

    t_in++;

    if (t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else
    {
        while (t_in < p)
        {
            pthread_cond_wait(&c_in, &my_mutex);
        }
    }

    if(r != a)
    {
        for (i = 0; i < n; i++)
        {
            a[i] = r[i];
        }
    }
    t_out++;
    if (t_out >= p)
    {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    }
    else
    {
        while (t_out < p)
        {
            pthread_cond_wait(&c_out, &my_mutex);
        }
        
    }
    pthread_mutex_unlock(&my_mutex);
}

class Args{
    public:
        int nx = 0, ny = 0, N = 0, len_msr = 0, func_id;
        double a = 0, b = 0, c = 0, d = 0, hx = 0, hy = 0;
        int p = 0, k = 0;

        double* A = nullptr;
        int* I = nullptr;
        double* B = nullptr;
        double *x = nullptr;
        double *r = nullptr;
        double *u = nullptr;
        double *v = nullptr;

        int its = 0, maxit = 0;
        double eps = 0;
        double t1 = 0, t2 = 0,r1 = 0, r2 = 0, r3 = 0, r4 = 0;

        pthread_t tid = 0;
        double cpu_time = 0;
        double cpu_time_of_all_threads = 0;
        double astr_time = 0.0;
        io_status status = io_status::not_working;
        double res = 0;
        /*
        void free_memory ()
        {
            if (A) delete [] A;
            if (I) delete [] I;
            if (B) delete [] B;
            if (x) delete [] x;
            if (r) delete [] r;
            if (u) delete [] u;
            if (v) delete [] v;
        }
        int allocate_memory()
        {
            free_memory();
            N = (nx + 1)*(ny + 1);
            len_msr =  N + 1 + get_len_msr(nx, ny);
            A = new double[len_msr];
            I = new int[len_msr];
            B = new double[N];
            x = new double[N];
            r = new double[N];
            u = new double[N];
            v = new double[N];
        }
        */
};

enum class io_status
{
    working,
    not_working,
    succes
};

class Args{
    public:
        int nx = 0, ny = 0, N = 0, len_msr = 0, func_id;
        double a = 0, b = 0, c = 0, d = 0, hx = 0, hy = 0;
        int p = 0, k = 0;

        double* A = nullptr;
        int* I = nullptr;
        double* B = nullptr;
        double *x = nullptr;
        double *r = nullptr;
        double *u = nullptr;
        double *v = nullptr;

        int its = 0, maxit = 0;
        double eps = 0;
        double t1 = 0, t2 = 0,r1 = 0, r2 = 0, r3 = 0, r4 = 0;

        pthread_t tid = 0;
        double cpu_time = 0;
        double cpu_time_of_all_threads = 0;
        double astr_time = 0.0;
        io_status status = io_status::not_working;
        double res = 0;

        /*
        void free_memory ()
        {
          if (A) delete [] A;
          if (I) delete [] I;
          if (B) delete [] B;
          if (x) delete [] x;
          if (r) delete [] r;
          if (u) delete [] u;
          if (v) delete [] v;
        }
        int allocate_memory()
        {
            free_memory();
            N = (nx + 1)*(ny + 1);
            len_msr =  N + 1 + get_len_msr(nx, ny);
            A = new double[len_msr];
            I = new int[len_msr];
            B = new double[N];
            x = new double[N];
            r = new double[N];
            u = new double[N];
            v = new double[N];
        }
        */
};

double get_cpu_time();
double get_full_time();
int minimal_resid_msr_matrix_full (int n, double *A, int *I, double *b,	double *x /*Начальное, а в конце будет ответ*/,	double *r,	double *u,	double *v,	double eps,	int max_it,	int max_step, int p, int k);
double Pf(double* res, double x, double y, double a, double c, double hx, double hy, int nx, int ny);
double calc_r1(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double));
double calc_r2(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double));
double calc_r3(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double));
double calc_r4(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double));

#endif 