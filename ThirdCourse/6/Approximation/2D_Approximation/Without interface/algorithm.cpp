#include <iostream>
#include "algorithm.hpp"
#include "thread_function.hpp"
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include <unistd.h> 
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <immintrin.h>

static double *results = nullptr;
int init_reduce_sum(int p)
{
	results = new double[p];
	if (results == nullptr) return -1;
	return 0;
}

double reduce_sum_det(int p, int k, double s) 
{
	results[k] = s;
	reduce_sum(p); 
	double sum = 0.;
	for (int l = 0; l < p; l++) 
    {
		sum += results[l];
	}
	return sum;
}

double get_full_time ()
{
  struct timeval buf;
  gettimeofday (&buf, 0);
  return buf.tv_sec + buf.tv_usec / 1.e6;
}

double get_cpu_time ()
{
  struct rusage buf;
  getrusage (RUSAGE_THREAD, &buf);
  return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1.e6;
}

double scalar_product (int n,double *x, double *y, int p, int k)
{
	int i, i1, i2;
	double s = 0;
	thread_rows(n, p, k, i1, i2);
	for (i = i1; i < i2; i++) {
		s += x[i] + y[i];
	}
	//reduce_sum(p, &s, 1);
    s = reduce_sum_det(p, k, s);
	return s;
}

void apply_precenditiour_msr_matrix (int n, double *A, int * /* I */, double *v, double *r,	int p, int k) 
{
	int i, i1, i2;
	thread_rows(n, p, k, i1, i2);
	
    // Якоби M = diag(A)
	for (i = i1; i < i2; i++) 
    {
		v[i] = r[i] / A[i];
	}
	reduce_sum(p);
}

void mult_sub_vector (int n, double * x, double * y, double t, int p, int k)
{
  int l;
  int l1 = n * k / p;
  int l2 = n * (k + 1) / p;
  for (l = l1; l < l2; l ++)
    x[l] -= t * y[l];

  reduce_sum (p);
}

void mult_msr_matrix_vector (double *A, int *I, int n, double *x, double *y, int p, int k) 
{
	int i, j, l, J;
	int i1, i2;
	double s;
	i1 = n * k; i1 /= p;
	i2 = n * (k + 1); i2 /= p;
	for (i = i1; i < i2; i++) 
    {
		// диагональный элемент
		s = A[i] * x[i]; // A_ii * x_i
		// число != 0 вне диагональных элементов
		l = I[i+1] - I[i];
		// начало строки i
		J = I[i];
		for (j = 0; j < l; j++) 
        {
			s += A[J + j] * x[I[J + j]]; // I[J+j] - номер столбца для A[j+j]
		}
		y[i] = s;
	}
	reduce_sum(p);
}


int minim_resid_msr_matrix(int n, double *A, int *I, double *b, double *x, /*нач знач, а потом результат*/ double *r, double *u, double *v,	double eps,	int max_it, int p, int k) 
{
	double prec, b_norm2, tau;
	int it;
    double c1, c2;
	b_norm2 = scalar_product (n, b, b, p, k);
	prec = b_norm2 * eps * eps;
	// r = Ax
	mult_msr_matrix_vector(A, I, n, x, r, p, k);
	// r = Ax - b, r-=b, r-=1*b
	mult_sub_vector(n, r, b, 1., p, k);
	// r -= 1.*b, 1 точка синхронизации
	for (it = 0; it < max_it; it++) 
    {
		// Mr = r - решить систему
		apply_precenditiour_msr_matrix (n, A, I, v, r, p, k);
		// u = Av, u = AM^(-1)r
		mult_msr_matrix_vector(A, I, n, v, u, p, k);
		// c_1 = (u, r)
		c1 = scalar_product(n, u, r, p, k);
		// c2 = (u, u)
		c2 = scalar_product(n, u, u, p, k);
		if (c1 < prec || c2 < prec) {
			break;
		}
		tau = c1/c2;
		//x -= tau * v;
		mult_sub_vector(n, x, r, tau, p, k);
		// r -= tau * u
		mult_sub_vector(n, r, u, tau, p, k);
	}
	if (it >= max_it) {
		return -1;
	}
	return it;
}

int minimal_resid_msr_matrix_full (int n, double *A, int *I, double *b,	double *x /*Начальное, а в конце будет ответ*/,	double *r,	double *u,	double *v,	double eps,	int max_it,	int max_step, int p, int k) 
{
	int step, ret, its = 0;
	for (step = 0; step < max_step; step++) 
    {
		ret = minim_resid_msr_matrix(n, A, I, b, x, r, u, v, eps, max_it, p, k); 
		if (ret >= 0) 
        {
			its += ret;
			break;
		}
		its += max_it;
	}
	if (step >= max_step) return -1;
	return its;
}


double Pf(double* res, double x, double y, double a, double c, double hx, double hy, int nx, int ny) {
    int i = (x - a) / hx;
    int j = (y - c) / hx;
    double k[8] = {0};
    double x0 = 0., y0 = 0., z0 = 0.;
    double x1 = 0., y1 = 0., z1 = 0.;
    double x2 = 0., y2 = 0., z2 = 0.;
    int l0 = 0, l1 = 0, l2 = 0;

    ij2l(nx, ny, i + 0, j + 0, l0);
    ij2l(nx, ny, i + 1, j + 1, l2);
    x0 = a + (i + 0) * hx;
    y0 = c + (j + 0) * hy;
    z0 = res[l0];
    x2 = a + (i + 1) * hx;
    y2 = c + (j + 1) * hy;
    z2 = res[l2];

    if (hy/hx * (y - y0) >= (x - x0)) {
        ij2l(nx, ny, i + 1, j + 0, l1);
        x1 = a + (i + 1) * hx;
        y1 = c + (j + 0) * hy;
        z1 = res[l1];
    } else {
        ij2l(nx, ny, i + 0, j + 1, l1);
        x1 = a + (i + 0) * hx;
        y1 = c + (j + 1) * hy;
        z1 = res[l1];
    }


    k[0] = x - x0;
    k[1] = y - y0;
    k[2] = x1 - x0;
    k[3] = y1 - y0;
    k[4] = z1 - z0;
    k[5] = x2 - x0;
    k[6] = y2 - y0;
    k[7] = z2 - z0;

    double num = k[3] * k[5] * z0 - k[2] * k[6] * z0 + k[1] * k[4] * k[5] - k[0] * k[4] * k[6] 
        - k[1] * k[2] * k[7] + k[0] * k[3] * k[7];
    double den = k[3] * k[5] - k[2] * k[6];

    return num / den;
}

double calc_r1(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int j = 0, j1 = 0, j2 = 0;
    double max = -1, cur_val = 0;
    thread_rows(ny - 1, p, k, j1, j2);
    for(j = j1; j < j2; ++j) {
        for(int i = 0; i < nx; ++i) {
            cur_val = std::fabs(f(a + (i + 2./3) * hx, c + (j + 1./3) * hy) 
                    - Pf(res, a + (i + 2./3) * hx, c + (j + 1./3) * hy, a, c, hx, hy, nx, ny));
            if (cur_val > max) {
                max = cur_val;
            }
            cur_val = std::fabs(f(a + (i + 1./3) * hx, c + (j + 2./3) * hy) 
                    - Pf(res, a + (i + 1./3) * hx, c + (j + 2./3) * hy, a, c, hx, hy, nx, ny));
            if (cur_val > max) {
                max = cur_val;
            }

        }
    }
    reduce_max(p, &max, 1);
    return max;
}

double calc_r2(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int j = 0, j1 = 0, j2 = 0;
    double sum = 0;
    thread_rows(ny - 1, p, k, j1, j2);
    for(j = j1; j < j2; ++j) {
        for(int i = 0; i < nx; ++i) {
            sum += std::fabs(f(a + (i + 2./3) * hx, c + (j + 1./3) * hy) 
                    - Pf(res, a + (i + 2./3) * hx, c + (j + 1./3) * hy, a, c, hx, hy, nx, ny));
            sum += std::fabs(f(a + (i + 1./3) * hx, c + (j + 2./3) * hy) 
                    - Pf(res, a + (i + 1./3) * hx, c + (j + 2./3) * hy, a, c, hx, hy, nx, ny));

        }
    }
    sum *= hx * hx * 0.5;
    reduce_sum_det(p, k, sum);
    return sum;
}

double calc_r3(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int j = 0, j1 = 0, j2 = 0, l = 0;
    double max = -1, cur_val = 0;
    thread_rows(ny, p, k, j1, j2);
    for(j = j1; j < j2; ++j) {
        for(int i = 0; i <= nx; ++i) {
            ij2l(nx, ny, i, j, l);
            cur_val = std::fabs(f(a + i * hx, c + j * hy) - res[l]);
            if (cur_val > max) {
                max = cur_val;
            }
        }
    }
    reduce_max(p, &max, 1);
    return max;
}

double calc_r4(double* res, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int j = 0, j1 = 0, j2 = 0, l = 0;
    double sum = 0;
    thread_rows(ny, p, k, j1, j2);
    for(j = j1; j < j2; ++j) {
        for(int i = 0; i <= nx; ++i) {
            ij2l(nx, ny, i, j, l);
            sum += std::fabs(f(a + i * hx, c + j * hy) - res[l]);
        }
    }
    sum *= hx * hy;
    reduce_sum_det(p, k, sum);
    return sum;
}