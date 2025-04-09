#include <iostream>
#include "algorithm.hpp"
#include "thread_function.hpp"
//#include <cmath>
#include <time.h>
#include <sys/time.h>
#include <unistd.h> 
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <immintrin.h>

double fabs(double a)
{
    return (a > 0 ? a : -a);
}

void reduce_sum(int p,double * a, int n)
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

void reduce_sum_int(int p,int * a , int n)
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
        double sum = 0.;
	reduce_sum(p); 
	for (int l = 0; l < p; l++) 
    {
		sum += results[l];
	}
    reduce_sum(p); 
	return sum;
}

void free_reduce_sum ()
{
  if (results)
    delete [] results;
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

#include <emmintrin.h> // Для использования SSE2

double scalar_product(int n, double *x, double *y, int p, int k) {
    int i, i1, i2;
    double s = 0;
    thread_rows(n, p, k, i1, i2);

    __m128d sum_vec = _mm_setzero_pd();
    for (i = i1; i <= i2 - 4; i += 4) { 
        __m128d x_vec1 = _mm_loadu_pd(&x[i]); 
        __m128d y_vec1 = _mm_loadu_pd(&y[i]); 
        __m128d x_vec2 = _mm_loadu_pd(&x[i + 2]); 
        __m128d y_vec2 = _mm_loadu_pd(&y[i + 2]); 
        sum_vec = _mm_add_pd(sum_vec, _mm_mul_pd(x_vec1, y_vec1)); 
        sum_vec = _mm_add_pd(sum_vec, _mm_mul_pd(x_vec2, y_vec2)); 
    }

    for (; i < i2; i++) {
        s += x[i] * y[i];
    }

    double temp[2];
    _mm_storeu_pd(temp, sum_vec); 
    s += temp[0] + temp[1]; 

    s = reduce_sum_det(p, k, s);

    return s;
}

void apply_precenditiour_msr_matrix(int n, double *A, int* , double *v, double *r, int p, int k) {
    int i, i1, i2;
    thread_rows(n, p, k, i1, i2);
    
    for (i = i1; i <= i2 - 4; i += 4) {
        __m128d r_vec1 = _mm_loadu_pd(&r[i]);
        __m128d A_vec1 = _mm_loadu_pd(&A[i]);
        __m128d r_vec2 = _mm_loadu_pd(&r[i + 2]);
        __m128d A_vec2 = _mm_loadu_pd(&A[i + 2]);
        __m128d v_vec1 = _mm_div_pd(r_vec1, A_vec1); 
        __m128d v_vec2 = _mm_div_pd(r_vec2, A_vec2); 
        _mm_storeu_pd(&v[i], v_vec1); 
        _mm_storeu_pd(&v[i + 2], v_vec2); 
    }

    for (; i < i2; i++) {
        v[i] = r[i] / A[i];
    }

    reduce_sum(p);
}

void mult_sub_vector(int n, double *x, double *y, double t, int p, int k) {
    int l, l1, l2;
    thread_rows(n, p, k, l1, l2);

    __m128d t_vec = _mm_set1_pd(t);
    for (l = l1; l <= l2 - 4; l += 4) {
        __m128d y_vec1 = _mm_loadu_pd(&y[l]); 
        __m128d x_vec1 = _mm_loadu_pd(&x[l]); 
        __m128d y_vec2 = _mm_loadu_pd(&y[l + 2]); 
        __m128d x_vec2 = _mm_loadu_pd(&x[l + 2]); 
        __m128d t_y_vec1 = _mm_mul_pd(t_vec, y_vec1); 
        __m128d t_y_vec2 = _mm_mul_pd(t_vec, y_vec2); 
        __m128d result_vec1 = _mm_sub_pd(x_vec1, t_y_vec1); 
        __m128d result_vec2 = _mm_sub_pd(x_vec2, t_y_vec2); 
        _mm_storeu_pd(&x[l], result_vec1); 
        _mm_storeu_pd(&x[l + 2], result_vec2); 
    }

    for (; l < l2; l++) {
        x[l] -= t * y[l];
    }

    reduce_sum(p);
}
void mult_msr_matrix_vector(double *A, int *I, int n, double *x, double *y, int p, int k) {
    int i, j, l, J;
    int i1, i2;
    double s;
    i1 = n * k; i1 /= p;
    i2 = n * (k + 1); i2 /= p;

    for (i = i1; i < i2; i++) {
        s = A[i] * x[i]; // A_ii * x_i

        l = I[i + 1] - I[i];
        J = I[i];

        __m128d sum_vec = _mm_setzero_pd();
        for (j = 0; j <= l - 4; j += 4) {
            __m128d A_vec = _mm_loadu_pd(&A[J + j]); 
            __m128d x_vec = _mm_set_pd(x[I[J + j + 1]], x[I[J + j]]); 
            sum_vec = _mm_add_pd(sum_vec, _mm_mul_pd(A_vec, x_vec));
            __m128d A_vec1 = _mm_loadu_pd(&A[J + j + 2]); 
            __m128d x_vec1 = _mm_set_pd(x[I[J + j + 3]], x[I[J + j + 2]]); 
            sum_vec = _mm_add_pd(sum_vec, _mm_mul_pd(A_vec1, x_vec1));
        }

        double temp[2];
        _mm_storeu_pd(temp, sum_vec);
        s += temp[0] + temp[1];

        for (; j < l; j++) {
            s += A[J + j] * x[I[J + j]];
        }

        y[i] = s;
    }

    reduce_sum(p);
}

double scalar_product_slow (int n,double *x, double *y, int p, int k)
{
	int i, i1, i2;
	double s = 0;
	thread_rows(n, p, k, i1, i2);
	for (i = i1; i < i2; i++) {
		s += x[i] * y[i];
	}
    s = reduce_sum_det(p, k, s);
	return s;
}

void apply_precenditiour_msr_matrix_slow(int n, double *A, int* , double *v, double *r, int p, int k) 
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

void mult_sub_vector_slow(int n, double * x, double * y, double t, int p, int k)
{
    int l1,l2;
    thread_rows(n, p, k, l1, l2);
    for (int l = l1; l < l2; l ++)
        x[l] -= t * y[l];

  reduce_sum (p);
}

void mult_msr_matrix_vector_slow(double *A, int *I, int n, double *x, double *y, int p, int k) 
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

int minim_error_msr_matrix(int n, double *A, int *I, double *b, double *x, /*нач знач, а потом результат*/ double *r, double *u, double *v,	double eps,	int max_it, int p, int k) 
{
	double prec = 0, b_norm2 = 0, tau = 0;
	int it = 0;
    double c1 = 0, c2 = 0;
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
		// c_1 = (v, r)
		c1 = scalar_product(n, v, r, p, k);
		// c2 = (u, v)
		c2 = scalar_product(n, u, v, p, k);
		if (c1 < prec || c2 < prec) {
			break;
		}
		tau = c1/c2;
		//x -= tau * v;
		mult_sub_vector(n, x, v, tau, p, k);
		// r -= tau * u
		mult_sub_vector(n, r, u, tau, p, k);
	}
	if (it >= max_it) {
		return -1;
	}
	return it;
}

int minimal_error_msr_matrix_full (int n, double *A, int *I, double *b,	double *x /*Начальное, а в конце будет ответ*/,	double *r,	double *u,	double *v,	double eps,	int max_it,	int max_step, int p, int k) 
{
	int step, ret, its = 0;
	for (step = 0; step < max_step; step++) 
    {
		ret = minim_error_msr_matrix(n, A, I, b, x, r, u, v, eps, max_it, p, k); 
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


//xi = a+ i*hx, hx = (b−a)/nx, y j = c+ j*hy, hy = (d− c)/ny
double Pf(double* res, double x, double y, double a, double c, double hx, double hy, int nx, int ny) {
    int i, j;
    int l0, l1, l2;
    double x0, y0, z0;
    double x1, y1, z1;
    double x2, y2, z2;
    double dx0, dy0, dx1, dy1, dz1, dx2, dy2, dz2;
    double num, den;

    i = (x - a) / hx;
    j = (y - c) / hy;

    ij2l(nx, ny, i + 0, j + 0, l0);
    ij2l(nx, ny, i + 1, j + 1, l2);

    x0 = a + (i + 0) * hx;
    y0 = c + (j + 0) * hy;
    z0 = res[l0];

    if (hy * (x - x0) / hx + y0 - y >= 0) {
        ij2l(nx, ny, i + 1, j + 0, l1);
        x1 = a + (i + 1) * hx;
        y1 = c + (j + 0) * hy;
    } else {
        ij2l(nx, ny, i + 0, j + 1, l1);
        x1 = a + (i + 0) * hx;
        y1 = c + (j + 1) * hy;
    }
    z1 = res[l1];

    ij2l(nx, ny, i + 1, j + 1, l2);
    x2 = a + (i + 1) * hx;
    y2 = c + (j + 1) * hy;
    z2 = res[l2];

    dx0 = x - x0;
    dy0 = y - y0;
    dx1 = x1 - x0;
    dy1 = y1 - y0;
    dz1 = z1 - z0;
    dx2 = x2 - x0;
    dy2 = y2 - y0;
    dz2 = z2 - z0;

    num = dy1 * dx2 * z0 - dx1 * dy2 * z0 + dy0 * dz1 * dx2 - dx0 * dz1 * dy2 
          - dy0 * dx1 * dz2 + dx0 * dy1 * dz2;
    den = dy1 * dx2 - dx1 * dy2;

    return num / den;
}
double R1Calculation(double* x, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int rowIndex = 0, startRow = 0, endRow = 0;
    double residual = -1, currentValue = 0;
    thread_rows(ny - 1, p, k, startRow, endRow);
    for (rowIndex = startRow; rowIndex < endRow; ++rowIndex) 
    {
        for (int columnIndex = 0; columnIndex < nx; ++columnIndex) 
        {
            currentValue = fabs(f(a + (columnIndex + 2.0 / 3) * hx, c + (rowIndex + 1.0 / 3) * hy) - Pf(x, a + (columnIndex + 2.0 / 3) * hx, c + (rowIndex + 1.0 / 3) * hy, a, c, hx, hy, nx, ny));
            if (currentValue > residual) 
            {
                residual = currentValue;
            }
            currentValue = fabs(f(a + (columnIndex + 1.0 / 3) * hx, c + (rowIndex + 2.0 / 3) * hy) - Pf(x, a + (columnIndex + 1.0 / 3) * hx, c + (rowIndex + 2.0 / 3) * hy, a, c, hx, hy, nx, ny));
            if (currentValue > residual) 
            {
                residual = currentValue;
            }
        }
    }
    reduce_max(p, &residual, 1);
    return residual; 
}

double R2Calculation(double* x, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int rowIndex = 0, startRow = 0, endRow = 0;
    double residual = 0;
    double sq = hx * hx/2;
    
    thread_rows(ny - 1, p, k, startRow, endRow);
    
    for (rowIndex = startRow; rowIndex < endRow; ++rowIndex) 
    {
        for (int columnIndex = 0; columnIndex < nx; ++columnIndex) 
        {
            residual += fabs(f(a + (columnIndex + 2.0 / 3) * hx, c + (rowIndex + 1.0 / 3) * hy) - Pf(x, a + (columnIndex + 2.0 / 3) * hx, c + (rowIndex + 1.0 / 3) * hy, a, c, hx, hy, nx, ny));
            residual += fabs(f(a + (columnIndex + 1.0 / 3) * hx, c + (rowIndex + 2.0 / 3) * hy) - Pf(x, a + (columnIndex + 1.0 / 3) * hx, c + (rowIndex + 2.0 / 3) * hy, a, c, hx, hy, nx, ny));
        }
    }
    residual *= sq;
    reduce_sum(p,&residual, 1);
    return residual; 
}

double R3Calculation(double* x, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int rowIndex = 0, startRow = 0, endRow = 0, linearIndex = 0;
    double residual = -1, currentValue = 0;

    thread_rows(ny, p, k, startRow, endRow);

    for (rowIndex = startRow; rowIndex < endRow; ++rowIndex) 
    {
        for (int columnIndex = 0; columnIndex <= nx; ++columnIndex) 
        {
            ij2l(nx, ny, columnIndex, rowIndex, linearIndex);
            currentValue = fabs(f(a + columnIndex * hx, c + rowIndex * hy) - x[linearIndex]);
            if (currentValue > residual) 
            {
                residual = currentValue;
            }
        }
    }
    reduce_max(p, &residual, 1);
    return residual; 
}

double R4Calculation(double* x, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
    int rowIndex = 0, startRow = 0, endRow = 0, linearIndex = 0;
    double residual = 0;
    double sq = hx * hy;

    thread_rows(ny, p, k, startRow, endRow);

    for (rowIndex = startRow; rowIndex < endRow; ++rowIndex) 
    {
        for (int columnIndex = 0; columnIndex <= nx; ++columnIndex) 
        {
            ij2l(nx, ny, columnIndex, rowIndex, linearIndex);
            residual += fabs(f(a + columnIndex * hx, c + rowIndex * hy) - x[linearIndex]);
        }
    }
    residual *= sq;
    reduce_sum(p, &residual, 1);
    return residual; 
}

void ResidualCalculation(double *r1, double *r2, double *r3, double *r4, double* x, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double))
{
    *r1 = R1Calculation(x, a, c, hx, hy, nx, ny, p, k, f);
    *r2 = R2Calculation(x, a, c, hx, hy, nx, ny, p, k, f);
    *r3 = R3Calculation(x, a, c, hx, hy, nx, ny, p, k, f);
    *r4 = R4Calculation(x, a, c, hx, hy, nx, ny, p, k, f);
}
