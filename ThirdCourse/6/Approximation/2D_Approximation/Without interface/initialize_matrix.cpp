#include <iostream>
#include "initialize_matrix.hpp"

void ij2l (int nx, int /*ny*/, int i, int j, int &l) 
{
	l = i + j *(nx + 1);
}

void l2ij (int nx, int /*ny*/, int &i,	int &j,	int l) 
{
	j = l / (nx + 1);
	i = l - j*(nx + 1);
}

int IA_ij(int nx, int ny, double hx, double hy, int i, int j, int is, int js, int s, int */*I = nullptr*/, double*A = nullptr)
{
    int l, ls; 
    double sq = hx*hy;
    ij2l(nx, ny, i, j, l);
    ij2l(nx, ny, is, js, ls);

    if (i > 0 && i < nx && j>0 && j<ny)
    {
        if (l == ls)
        {
            A[s] = 6*sq/12;
        }
        else
        {
            A[s] = 2*sq/24;
        }
    }
    if(j == 0 && i>0 && i<nx)
    {
        if(l == ls)
        {
            A[s] = 3*sq/12;
        }
        else if(s == 0 || s == 1)
        {
            A[s] = sq/24;
        }
        else if(s == 2 || s == 3)
        {
            A[s] = 2*sq/24;
        }
        else
        {
            return -1;
        }
    }
    if(j == ny && i>0 && i<nx)
    {
        if(l == ls)
        {
            A[s] = 3*sq/12;
        }
        else if(s == 0 || s == 3)
        {
            A[s] = sq/24;
        }
        else if(s ==1 || s == 2)
        {
            A[s] = 2*sq/24;
        }
        else
        {
            return -1;
        }
    }
    if(i == 0&& j > 0 && j< ny)
    {
        if(l == ls)
        {
            A[s] = 3*sq/12;
        }
        else if(s == 0 || s == 3)
        {
            A[s] = 2*sq/24;
        }
        else if( s == 1 || s == 2)
        {
            A[s] = sq/24;
        }
        else 
        {
            return -1;
        }
    }
    if(i == nx && j > 0 && j < ny)
    {
        if(l == ls)
        {
            A[s] = 3*sq/12;
        }
        else if(s == 0 || s == 3)
        {
            A[s] = sq/24;
        }
        else if(s ==1 || s == 2)
        {
            A[s] = 2*sq/24;
        }
        else
        {
            return -1;
        }
    }
    if(i == 0 && j == 0)
    {
        if(l == ls)
        {
            A[s] = 2*sq/12;
        }
        else if(s == 0 || s == 1)
        {
            A[s] = sq/24;
        }
        else if(s == 2)
        {
            A[s] = 2*sq/24;
        }
        else
        {
            return -1;
        }
    }
    if(i == nx && j == ny)
    {
        if(l == ls)
        {
            A[s] = 2*sq/12;
        }
        else if(s == 0 || s == 2)
        {
            A[s] = sq/24;
        }
        else if(s == 1)
        {
            A[s] = 2*sq/24;
        }
        else
        {
            return -1;
        }
    }
    if((i == 0 && j == ny) || (i == nx && j == 0))
    {
        if(l == ls)
        {
            A[s] = sq/12;
        }
        else if(s == 0 || s == 1)
        {
            A[s] = sq/24;
        }
        else
        {
            return -1;
        }
    }

    return -1;
}

#define F(IS, JS, S) \
    (IA_ij(nx, ny, hx, hy, i, j, (IS), (JS), (S), I, A))

void thread_rows (int n, int p, int k, int &i1, int &i2)
{
	i1 = n*k;
	i1 /= p;
	i2 = n*(k+1);
	i2 /= p;
}

// Число внедиагональных элементов
// Общая длина msr матрицы = диаг + 1 + число внедиагональных
// (nx + 1) * (ny - 1) + 1 + get_len_msr(nx, ny)
int get_len_msr (int nx, int ny)
{
	return 6 * (nx - 1) * (ny - 1) + 4 * (2 * (nx - 1) + 2 * (ny - 1)) + 3 * 2 + 2 * 2;
}

//возвращает количество внедиагональных точки (i, j)
int get_off_diag(int nx, int ny, double hx, double hy, int i, int j, int *I, double *A)
{
    int s = 0;
    if(i < nx)
    {
        F(i + 1, j, s); s++;
    }
    if(j > 0)
    {
        F(i, j-1,s); s++;
    }
    if(i > 0 && j > 0)
    {
        F(i-1, j-1, s); s++;
    }
    if(i > 0)
    {
        F(i-1, j, s); s++;
    }
    if(j < ny)
    {
        F(i, j+1, s); s++;
    }
    if(i < nx && j < ny)
    {
        F(i+1, j+1, s); s++;
    }
    return s;//Amount of not diagonal elements
}

int get_len_msr_off_diag(int nx, int ny, double *A, int *I)
{
    double hx = 0, hy = 0;
    int i, j, res = 0;

    for(i = 0; i<=nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            res += get_off_diag(nx, ny, hx, hy, i, j, I, A);
        }
    }
    return res;
}

void fill_I(int nx, int ny, int *I)
{
    int i,j,l,N = (nx+1)*(ny+1);
    int r = N + 1;
    double hx = 0,hy = 0;

    for (l = 0; l < N; l++)
    {
        l2ij(nx,ny,i,j,l);
        int s = get_off_diag(nx, ny, hx, hy, i, j, 0, 0);
        I[l] = r;
        r += s;
    }
    I[l] = r;
}

int get_diag(int nx, int ny, double hx, double hy, int i, int j, int */*I*/, double* A)
{
    //tut mb nuzhna l
    return IA_ij(nx, ny, hx, hy, i, j, i, j, 0, nullptr, A);
}

int fill_IA(int nx, int ny, double hx, double hy, int *I, double* A, int p, int k)
{
    int i, j, l,l1,l2,N = (nx+1)*(ny+1), r,s,t;
    double err = 0; int len = 0;
    thread_rows(N,p,k,l1,l2);

    for (l = l1; l < l2; l++)
    {
        r = I[l];
        s = I[l+1]-I[l];
        l2ij(nx, ny, i, j, l);

        //Tut mb nado prosto A
        if(get_diag(nx, ny, hx, hy, i, j, I, A+l)!=0)
        {
            err = 1;
            break;
        }
        t = get_off_diag(nx, ny, hx, hy, i, j, I + r, A + r);
        if(t != s)
        {
            err = -1;
            break;
        }
        len += s;
    }
    reduce_sum(p, &err, 1);
    if(err< 0)
    {
        return -1;
    }
    reduce_sum_int(p, &len, 1);
    if(I[N] != N+1+len)
    {
        return -2;
    }

    return 0;
}

#define Fb(I, J) (f(x0 + (I) * hx, y0 + (J) * hy))

double F_IJ(int nx, int ny, double hx, double hy, double x0, double y0, int i, int j, double (*f)(double, double)) {
    double w = hx * hy / 192;
    if (i > 0 && i < nx && j > 0 && j < ny) {
        return w * (
                36 * Fb(i, j)
                + 20 * (
                    Fb(i + 0.5, j) + Fb(i, j - 0.5) + Fb(i - 0.5, j - 0.5)
                    + Fb(i - 0.5, j) + Fb(i, j + 0.5) + Fb(i + 0.5, j + 0.5)
                    )
                + 4 * (
                    Fb(i + 0.5, j - 0.5) + Fb(i - 0.5, j - 1) + Fb(i - 1, j - 0.5)
                    + Fb(i - 0.5, j + 0.5) + Fb(i + 0.5, j + 1) + Fb(i + 1, j + 0.5)
                    )
                + 2 * (
                    Fb(i + 1, j) + Fb(i, j - 1) + Fb(i - 1, j - 1) 
                    + Fb(i - 1, j) + Fb(i, j + 1) + Fb(i + 1, j + 1)
                    )
                );
    }
    if (i > 0 && i < nx && j == 0) {
        return w * (
                18 * Fb(i, j)
                + 10 * (Fb(i + 0.5, j) + Fb(i - 0.5, j))
                + 20 * (Fb(i, j + 0.5) + Fb(i + 0.5, j + 0.5))
                + 4 * (Fb(i - 0.5, j + 0.5) + Fb(i + 0.5, j + 1) + Fb(i + 1, j + 0.5))
                + 1 * (Fb(i - 1, j) + Fb(i + 1, j))
                + 2 * (Fb(i, j + 1) + Fb(i + 1, j + 1))
                );
    }
    if (i > 0 && i < nx && j == ny) {
        return w * (
                18 * Fb(i, j)
                + 10 * (Fb(i - 0.5, j) + Fb(i + 0.5, j))
                + 20 * (Fb(i, j - 0.5) + Fb(i - 0.5, j - 0.5))
                + 4 * (Fb(i + 0.5, j - 0.5) + Fb(i - 0.5, j - 1) + Fb(i - 1, j - 0.5))
                + 1 * (Fb(i - 1, j) + Fb(i + 1, j))
                + 2 * (Fb(i, j - 1) + Fb(i - 1, j - 1))
                );
    }
    if (i == 0 && j > 0 && j < ny) {
        return w * (
                18 * Fb(i, j)
                + 10 * (Fb(i, j - 0.5) + Fb(i, j + 0.5))
                + 20 * (Fb(i + 0.5, j) + Fb(i + 0.5, j + 0.5))
                + 4 * (Fb(i + 0.5, j - 0.5) + Fb(i + 0.5, j + 1) + Fb(i + 1, j + 0.5))
                + 1 * (Fb(i, j - 1) + Fb(i, j + 1))
                + 2 * (Fb(i + 1, j) + Fb(i + 1, j + 1))
                );
    }
    if (i == nx && j > 0 && j < ny) {
        return w * (
                18 * Fb(i, j)
                + 10 * (Fb(i, j - 0.5) + Fb(i, j + 0.5))
                + 20 * (Fb(i - 0.5, j) + Fb(i - 0.5, j + 0.5))
                + 4 * (Fb(i - 0.5, j) + Fb(i - 1, j - 0.5) + Fb(i - 0.5, j + 0.5))
                + 1 * (Fb(i, j - 1) + Fb(i, j + 1))
                + 2 * (Fb(i - 1, j) + Fb(i - 1, j - 1))
                );
    }
    if (i == 0 && j == 0) {
        return w * (
                12 * Fb(i, j)
                + 10 * (Fb(i + 0.5, j) + Fb(i, j + 0.5))
                + 20 * (Fb(i + 0.5, j + 0.5))
                + 4 * (Fb(i + 1, j + 0.5) + Fb(i + 0.5, j + 1))
                + 1 * (Fb(i + 1, j) + Fb(i, j + 1))
                + 2 * (Fb(i + 1, j + 1))
                );
    }
    if (i == nx && j == ny) {
        return w * (
                12 * Fb(i, j)
                + 10 * (Fb(i - 0.5, j) + Fb(i, j - 0.5))
                + 20 * (Fb(i - 0.5, j - 0.5))
                + 4 * (Fb(i - 0.5, j - 1) + Fb(i - 1, j - 0.5))
                + 1 * (Fb(i, j - 1) + Fb(i - 1, j))
                + 2 * (Fb(i - 1, j - 1))
                );
    }
    if (i == 0 && j == ny) {
        return w * (
                6 * Fb(i, j)
                + 10 * (Fb(i + 0.5, j) + Fb(i, j - 0.5))
                + 4 * (Fb(i + 0.5, j - 0.5))
                + 1 * (Fb(i + 1, j) + Fb(i, j - 1))
                );
    }
    if (i == nx && j == 0) {
        return w * (
                6 * Fb(i, j)
                + 10 * (Fb(i - 0.5, j) + Fb(i, j + 0.5))
                + 4 * (Fb(i - 0.5, j + 0.5))
                + 1 * (Fb(i - 1, j) + Fb(i, j + 1))
                );
    }

    return 1e308;
}

void fill_B(double a, double c, int nx, int ny, double hx, double hy, double* b, int p, int k, double (*f)(double, double))
{
    int l1 = 0, l2 = 0, l = 0;
    int N = (nx + 1) * (ny + 1);
    int i = 0, j = 0;
    l1 = N * k;
    l1 /= p;
    l2 = N * (k + 1);
    l2 /= p;
    for(l = l1; l < l2; ++l) {
        l2ij(nx, ny, i, j, l);
        b[l] = F_IJ(nx, ny, hx, hy, a, c, i, j, f); 
    }
    reduce_sum(p);
}
