#include "chebyshov_approximation.hpp"
#include <iostream>
#include <cmath>
#define PI 3.14159265358980

void ChebyshovAproximation(int n, double* f, double* a, double*g, double* g2, double* z)
{
    double* pg0 = g, *pg1 = g2;
    double gij = 0, res = 0;
    double zj;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == 0)
            {
                gij = f[j];
                pg0[j] = gij;
            }            
            else if(i == 1)
            {
                zj = cos(((2*j+1)*PI)/(2*n));
                gij = zj*f[j];
                pg1[j] = gij;
                z[j] = 2*zj;
            }
            else
            {
                if(i%2 == 0)
                {
                    gij = z[j]*pg1[j] - pg0[j];
                    pg0[j] = gij;
                }
                else
                {
                    gij = z[j]*pg0[j] - pg1[j];
                    pg1[j] = gij;
                }
            }
            res += gij;

        }
       
        if (i == 0)
        {
            a[i] = res/n;
        }
        else
        {
            a[i] = 2*res/n;
        }        
        
        res = 0;
    }
}

double ChebyshovValue(double X, double A, double B, int n, double* a)
{
    X/=2;
    double z = (2*X-(B+A))/(B-A);
    double T_i_2 = 1, T_i_1 = z, T_i=0;
    double Pf;
    
    Pf = a[0];
    if (n > 0)
    {
        Pf += a[1]*T_i_1;
    }
    z*=2;
    for (int i = 2; i < n; i++)
    {
        T_i = z*T_i_1 - T_i_2;  
        Pf += a[i]*T_i;
        T_i_2 = T_i_1;
        T_i_1 = T_i;  
    }

    return Pf;
}

