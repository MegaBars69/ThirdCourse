#include <iostream>
#include <cmath>
#include "spline_approximation.hpp"

#include <cstring>
#define EPSILON pow(10,-15)
using namespace std;
double derivativeF(double (*f)(double), double x)
{
    double h = 1e-8;

    return (f(x+h) - f(x-h)) / (2*h);
}
void GaussSolve(double* A, double* b, int n) 
{
    int i,j, total = 3*n - 1;
    double aii, a;
    for (i = 0, j = 0; i < n; i++)
    {
        if(i == 1)
        {
            a = A[j-1];
            A[j] -= a * A[j-3];
            b[i] -= b[i-1]*a; 
        }
        else if(i > 0)
        {
            a = A[j-1];
            A[j] -= a * A[j-2];
            b[i] -= b[i-1]*a; 
        }
        
        aii = A[j];
        if(fabs(aii)< EPSILON)
        {
            cout<<"Method can't be applied"<<endl;
            return;
        }
        if(i != n-1)
        {
            A[j + 1] /= aii;
        }

        b[i] /= aii;
        if(i != 0)
        {
            j+=3;
        }
        else
        {
            j+=4;
        }
    
    }

    //GaussBackward
    for(i = n - 2, j = total - 3; i>= 0; i--, j-=3)
    {
        if(i == 0)
        {
            j--;
        }
        b[i] -= b[i+1]*A[j];
    }
}

void CreateMatrix(double* A, double* b, double*dF, double* X, double left_derivative, double right_derivative, int n)
{
    int i, j;
    double dx = 0, dx1, dx2, xi = X[0], xip1 = X[1],xip2 = X[2];
    for (i = 0, j = 0; i < n; i++, j += 3)
    {
        xip2 = X[i + 2];
        if(i == 0)
        {
            dx = X[1] - X[0];
            A[j] = 2*(X[2]-X[0]);
            A[j+1] = dx;
            A[j + 2] = 0;
            b[i] = 3*(dF[i]*(X[2] - X[1]) + dF[i+1]*dx) - left_derivative*(X[2] - X[1]);
        }
        else if(i == n-1)
        {
            dx2 = xip2 - xip1; 
            dx1 = xip1 - xi;
            A[j] = dx2;
            A[j + 1] = 2*(xip2 - xi);
            b[i] = 3*(dF[i] * dx2 + dF[i + 1]*dx1) - dx1*right_derivative;
        }
        else
        {
            dx2 = xip2 - xip1; 
            dx1 = xip1 - xi;
            A[j] = dx2;
            A[j + 1] = 2*(xip2 - xi);
            A[j + 2] = dx1;
            b[i] = 3*(dF[i] * dx2 + dF[i + 1]*dx1);
        }
        xi = xip1;
        xip1 = xip2;
    }
    
}

void CalculateParametrs(double* A, double*d, double*dF, double* X, double left_derivative, double right_derivative, int n)
{
    d[0] = left_derivative, d[n-1] = right_derivative;
    CreateMatrix(A, d+1, dF, X, left_derivative, right_derivative, n-2);
    GaussSolve(A, d+1, n-2);
}

void CalculateDiferences(double *dF, double *X, double *F, int n)
{
    double xi, xi_1, fi,fi_1;
    xi_1 = X[0];
    fi_1 = F[0];
    for (int i = 1; i < n; i++)
    {
        xi = X[i];
        fi = F[i];
        dF[i-1] = (fi - fi_1) / (xi - xi_1);
        xi_1 = xi;
        fi_1 = fi;
    }
    
}

void CalculateCoeficients(double* c, double* X, double* d, double* dF, double* F, int n)
{
    double dxi = 0;
    int i,j;
    for(i = 0, j = 0; i < n; i++, j += 4)
    {
        if(i != n-1)
        {
            c[j] = F[i];
            c[j + 1] = d[i];
            dxi = X[i+1] - X[i];
            c[j + 2] = (3*dF[i] - 2*d[i] - d[i+1])/dxi;
            dxi *= dxi;
            c[j + 3] = (d[i] + d[i + 1] - 2*dF[i])/dxi;
        }
        else
        {
            c[j] = F[i];
            c[j + 1] = 0;
            c[j + 2] = 0;
            c[j + 3] = 0;
        }
    }
}

int FindI(double x, double a, double b, int n) 
{
    if(fabs(x - a) < EPSILON || x < a)
    {
        return 0;
    }
    else if(fabs(x - b) < EPSILON || x > b)
    {
        return n-1;
    }
    

    double step = (b - a) / (n-1);

    int left = 0;
    int right = n;

    while (left < right) {
        int mid = left + (right - left) / 2;
        double midPoint = a + mid * step;

        if (x <= midPoint) {
            right = mid; 
        } else {
            left = mid + 1; 
        }
    }

    return left - 1; 
}


double SplineValue(double x, double a, double b, int n, double* c, double* X, int i)
{
    double Pf = 0, dx = 0;
    int j;
    if(i == -1)
    {
        i = FindI(x, a, b, n);
        if(i > n-1)
        {
            i = n-1;
        }
    }
    
    j = i*4;
    Pf += c[j];
    dx  = (x-X[i]);
    Pf += c[j + 1] * dx;
    Pf += c[j + 2] * dx*dx;
    Pf += c[j + 3] * dx*dx*dx;
    return Pf;
}
