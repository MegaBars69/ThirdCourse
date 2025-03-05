#include "approximation.hpp"
#include <iostream>
#include <cmath>
#define PI 3.14159265358980

void ChebyshovAproximation(int n, double* f, double* a, double*g, double* z)
{
    double* pg0 = g, *pg1 = g + n;
    double gij = 0, res = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; i++)
        {
            if (i == 0)
            {
                gij = f[j];
                pg0[j] = gij;
            }            
            else if(i == 1)
            {
                z[j] = 2*cos((2*j-2)*(PI/n));
                gij = z[j]*f[j]/2;
                pg1[j] = gij;
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
    }
}

double ChebyshovValue(double X, double A, double B, int n, double* a)
{
    double z = 2*(2*X-(B+A))/(B-A);
    double T_i_2 = 1, T_i_1 = z/2, T_i=0;
    double Pf;
    Pf = a[0]*T_i_2;
    if (n > 0)
    {
        Pf += a[1]*T_i_1;
    }
    
    for (int i = 2; i < n; i++)
    {
        T_i = z*T_i_1 - T_i_2;  
        Pf += a[i]*T_i;
        T_i_2 = T_i_1;
        T_i_1 = T_i;  
    }

    return Pf;
}

