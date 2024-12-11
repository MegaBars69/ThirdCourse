#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <pthread.h>
#include <cmath>
#include "func.hpp"

int main(int argc, char* argv[])
{
    Args* a;
    int k = 0;
    int chunk = 20000;

    if(argc < 3)
    {
        printf("Usage ./a.out p n\n");
        return 1;
    }
    int n = atoi(argv[2]);
    int p = atoi(argv[1]);  

    if (p > n)
    {
        p = n; 
    }
    if (p <= 0)
    {
        printf(" p > 0 \n");
        return 4;
    }

    a = new Args[p];
        
    //Обработка массива. Распределение памяти.
    for (k = 0; k < p; k++)
    {
        a[k].k = k;
        a[k].p = p;
        a[k].n = n;
        a[k].chunk = chunk;
        a[k].L = 5 + k*chunk;
    }
    
    //Запуск потоков.
    for (k = 1; k < p; k++)
    {
        if (pthread_create(&a[k].tid, nullptr, thread_func, a+k))
        {
            std::cerr << "Error creating thread " << k <<std::endl;
            delete[] a;

            return 1;
        }
    }
    a[0].tid = pthread_self();

    thread_func(a+0);
 
    //Прибиваем потоки молотком.
    for (k = 1; k < p; k++)
    {
        pthread_join(a[k].tid, nullptr);
    }
    
    for (k = p-1; k >= 0; k--)
    {
        a[k].PrintAll();
    }
    
   
    
    
    delete[] a;
    return 0;
}