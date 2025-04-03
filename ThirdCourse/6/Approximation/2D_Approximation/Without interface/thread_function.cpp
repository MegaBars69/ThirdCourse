#include <iostream>
#include <cmath>
#include "algorithm.hpp"
#include "initialize_matrix.hpp"

#include "thread_function.hpp"

double (*get_funk(int num))(double, double) 
{
    static auto f0 = [](double /*x*/, double /*y*/) {
        return 1.;
    };
    static auto f1 = [](double x, double /*y*/) {
        return x;
    };
    static auto f2 = [](double /*x*/, double y) {
        return y;
    };
    static auto f3 = [](double x, double y) {
        return x + y;
    };
    static auto f4 = [](double x, double y) {
        return std::sqrt(x * x + y * y);
    };
    static auto f5 = [](double x, double y) {
        return x * x + y * y;
    };
    static auto f6 = [](double x, double y) {
        return std::exp(x * x - y * y);
    };
    static auto f7 = [](double x, double y) {
        return 1./(25. * (x * x + y * y) + 1);
    };

    switch (num) {
        case 0: 
            return f0;
        case 1:
            return f1;
        case 2:
            return f2;
        case 3:
            return f3;
        case 4:
            return f4;
        case 5:
            return f5;
        case 6:
            return f6;
        case 7:
            return f7;
        default :
            return nullptr;
    }
    return nullptr;
}

void* thread_func(void *arg)
{  
    Args *  aa = (Args *) arg;

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

    double (*f)(double, double) = get_funk(func_id);

    double hx = (b - a)/nx, hy = (d - c)/ny;
    aa->hx = hx, aa->hy = hy;
    int N = aa->N;
    double t1 = -1, t2 = -1, r1 = -1, r2 = -1, r3 = -1, r4 = -1, err = 0;

    if(k == 0)
    {
        fill_I(nx, ny, I);
    }

    err = fill_IA(nx, ny, hx, hy, I, A, p, k);
    reduce_sum(p, &err, 1);
    if(fabs(err) > 0)
    {
        return nullptr;
    }

    fill_B(a, c, nx, ny, hx, hy, B, p, k, f);
    
    reduce_sum(p);

    t1 = get_full_time();
    its = minimal_resid_msr_matrix_full(N, A, I, B, x, r, u, v, eps, maxit, maxit, p, k); 
    reduce_sum(p);
    t1 = get_full_time() - t1;

    reduce_sum(p);
    t2 = get_full_time();
    r1 = calc_r1(x, a, c, hx, hy, nx, ny, p, k, f);
    r2 = calc_r2(x, a, c, hx, hy, nx, ny, p, k, f);
    r3 = calc_r3(x, a, c, hx, hy, nx, ny, p, k, f);
    r4 = calc_r4(x, a, c, hx, hy, nx, ny, p, k, f);
    reduce_sum(p);
    t2 = get_full_time() - t2;

    aa->t1 = t1;
    aa->t2 = t2;
    aa->r1 = r1;
    aa->r2 = r2;
    aa->r3 = r3;
    aa->r4 = r4;
    aa->its = its;

    return nullptr;
}