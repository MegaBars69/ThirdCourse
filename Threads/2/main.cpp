#include <iostream>
#include <fstream>
#include <pthread.h>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include "func.hpp"

int main(int argc, char* argv[])
{
    Args* a;
    if(argc == 1)
    {
        printf("Usage %s files ", argv[0]);
        return 1;
    }

    int p = argc - 1;
    a = new Args[p];

    int k;
    for (k = 0; k < p; k++)
    {
        a[k].k = k;
        a[k].p = p;
        a[k].filename = argv[k+1];
    }
    for (k = 1; k < p; k++)
    {
        if (pthread_create(&a[k].tid, nullptr, thread_func, a+k))
        {
            std::cerr << "Error creating thread " << k <<std::endl;
            return 1;
        }
    }
    a[0].tid = pthread_self();

    thread_func(a+0);

    for (k = 1; k < p; k++)
    {
        pthread_join(a[k].tid, nullptr);
    }
    
    if (procces_args(a) == 0)
    {
        printf("Result = %d \n", int(a[0].res));
    }
    
    delete[] a;
    return 0;
}