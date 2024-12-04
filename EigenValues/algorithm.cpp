#include <iostream>
#include <cmath>
#include "algorithm.hpp"

double Trace(double* A, int n)
{
    double result = 0;
    for (int i = 0; i < n; i++)
    {
        result += A[i*n + i];
    }
    return result;
}

double Norm(double* A, double* results, int n)
{
    int m = n;
    int l = n%m;
    int k = (n-l)/m;
    int block_size_col, block_size_row;
    double sum_in_col = 0, final_norm = 0, el;

    for (int bj = 0; bj < k+1; bj++)
    {
        block_size_col = (bj < k ? m : l);
        for (int bi = 0; bi < k+1 ; bi++)
        {
            block_size_row = (bi < k ? m : l);

            double* pa = (A + m*n*bi + m*block_size_row*bj);
            for (int j = 0; j < block_size_col; j++)
            {
                sum_in_col = 0;
                for (int i = 0; i < block_size_row; i++)
                {
                    sum_in_col += fabs(pa[i*block_size_col + j]);
                }
                results[j] += sum_in_col;
            }
        }
        for (int j = 0; j < block_size_col; j++)
        {
            el = results[j];
            results[j] = 0;
            final_norm = (el < final_norm ? final_norm : el);
        }   

    }   
    return final_norm;
}

void ApplyLeftVector(double* X, double* A, int n, int k, bool inside)
{
    double s;
    double* px;
    for (int j = (inside ? k+1: 0); j < n; j++)
    {
        s = 0;
        px = X;
        for (int i = k + 1 ; i < n; i++,px++)
        {
            s += (*px) * (A[i*n + j]);
        }

        px = X;
        s *= 2;
        for (int i = k + 1; i < n; i++,px++)
        {
            A[i*n + j]-= s*(*px);
        }        
    }   
}

void AppyRightVector(double* X, double* A, int n, int k, bool inside)
{
    double s;
    double* px;
    for (int j = (inside ? k+1: 0); j < n; j++)
    {
        s = 0;
        px = X;
        for (int i = k + 1 ; i < n; i++,px++)
        {
            s += (*px) * (A[j*n + i]);
        }

        px = X;
        s *= 2;
        for (int i = k + 1; i < n; i++,px++)
        {
            A[j*n + i]-= s*(*px);
        }        
    }   
}


void TriDiagonalize(double* A, double* U, int n, double mera)
{
    double sk, ajk,akk;
    double new_diag_el, norm_xk, first_in_x;
    double* pu;

    for (int k = 0; k < n - 1; k++)
    {
        //Finding vector xk
        pu = U;
        sk = 0;

        for (int i = k+2; i < n; i++)
        {
            ajk = A[i*n + k];
            pu++;
            *pu = ajk;
            sk+= ajk*ajk;
            A[i*n + k] = 0;
            A[k*n + i] = 0;
        }

        akk = A[(k + 1)*n + k];
        new_diag_el = sqrt(sk + akk*akk);
        first_in_x = akk - new_diag_el;
        
        pu = U;
        *pu = first_in_x;
        
        norm_xk = sqrt(sk + first_in_x*first_in_x);

        if(!(fabs(norm_xk) < mera))
        {
            for (int i = 0; i < n - k - 1; i++, pu++)
            {
                *pu = (*pu)/norm_xk;
            }
        }

        //Applying Reflection
        A[(k + 1)*n + k] = new_diag_el;
        A[k*n + k + 1] = new_diag_el;

        ApplyLeftVector(U, A, n, k, true); 
        AppyRightVector(U, A, n, k, true); 
    }

}