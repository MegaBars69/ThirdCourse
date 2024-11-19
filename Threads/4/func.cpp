#include <iostream>
#include <fstream>
#include <pthread.h>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include "func.hpp"
#define EPSILON pow(10, -16)

void UpdateElements(Args* a)
{
    double* pa = a->array;
    double a0,a1,a2, a3,cur;
    int n = a->m;
    int k = a->k;
    int p = a->p;
    int result = 0;
    double new_el = 0;

    if(a->n < 4)
    {
        return;
    }
    if (n == 1)
    {
        if ((k > 0) && (((k < p-1) && (a->l > 0)) || (k < p-2)))
        {
            new_el = (a->pan + pa[0] + a->na0 + a->na1)/4;
            pa[0] = new_el;
            result++;
        }
    }

    else if (n == 2)
    {
        if(k > 0 && k < p-1)
        {
            cur = pa[0];
            new_el = (a->pan + cur + pa[1] + a->na0)/4;
            pa[0] = new_el;
            result++;
            new_el = (cur + pa[1]+ a->na0 + a->na1)/4;
            pa[1] = new_el; 
            result++;
        }
        else if(k == 0 )
        {
            pa[1] = (pa[0] + pa[1] + a->na0 + a->na1)/4;
            result++;
        }
    }
    else if(n >= 3) 
    {
        if (k > 0)
        {
            a0 = a->pan;
            a1 = pa[0];
            a2 = pa[1];
            a3 = pa[2];
            new_el = (a0 + a1 +a2 + a3)/4;
            pa[0] = new_el;
            result++;

            a0 = a1;
            a1 = a2;
            a2 = a3;
        }
        else
        {
            a0 = pa[0];
            a1 = pa[1];
            a2 = pa[2];
        }

        for (int i = 3; i < n; i++)
        {
            a3 = pa[i];
            new_el = (a0+a1+a2+a3)/4;
            pa[i-2] = new_el;
            result++;

            a0 = a1;
            a1 = a2;
            a2 = a3;
        }
        if (k < p-1)
        {
            a3 = a->na0;
            new_el = (a0 + a1 +a2 + a3)/4;
            pa[n-2] = new_el;
            result++;

            a0 = a1;
            a1 = a2;
            a2 = a3;
            a3 = a->na1;
            
            new_el = (a0 + a1 +a2 + a3)/4;

            pa[n-1] = new_el;
            result++;
        }        
    }
    a->amount_of_changed = result; 
}


void* thread_func(void *arg)
{
    Args *a = (Args *)arg;

    UpdateElements(a);

    reduce_sum(a->p, &a->amount_of_changed, 1);

    return nullptr;    
}

