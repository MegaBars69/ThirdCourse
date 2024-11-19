#include <iostream>
#include <fstream>
#include <pthread.h>
#include <cmath>
#include "func.hpp"

int main(int argc, char* argv[])
{
    Args* a;
    double *array = nullptr, *pa = nullptr;
    FILE *fp = fopen(argv[3], "r");
    double el;
    int i = 0, k;

    if(argc == 1)
    {
        printf("Usage ./a.out p n filename\n");
        return 1;
    }
    int n = atoi(argv[2]);
    int p = atoi(argv[1]); 
    printf("RESULT %2d:",p);

    if (p > n)
    {
        p = n; 
    }
    if (p <= 0)
    {
        printf(" p > 0 \n");
        return 4;
    }
    
    int m = n/p;
    int l = n%p;

    if (n > 0)
    {
        array = new double[n];
    }
    else
    {
        printf(" n >=0");
        return 4;
    }        


    if(fp == nullptr)
    {
        printf("Can't open %s \n",argv[2]);
        delete[] array;
        return 1;
    }

    //Чтение массива
    while (true) 
    {
        int res = fscanf(fp, "%lf", &el);
        if(res == EOF)
        {  
            if(i != n)
            {
                printf("To few elements in file(<n)\n");
                delete[] array;
                return 3;
            }
            break;
        }
        else if (res== 0)
        {
            printf("Problem with reading an element \n");
            
            delete[] array;
            return 2;
        }

        if (i > n)
        {
            printf("To many elements in file(>n)\n");
            delete[] array;
            return 3;
        }
        array[i] = el;
        i++;
    
    }
    a = new Args[p];
    pa = array;
    //Обработка массива. Распределение памяти.
    for (k = 0; k < p; k++,pa += m)
    {
        a[k].k = k;
        a[k].p = p;
        a[k].n = n;
        a[k].l = l;
        a[k].m = (k < p-1 ? m : m+l);
        if (k > 0)
        {
            a[k].pan = *(pa-1);
        }
        if (k < p-1)
        {
            a[k].na0 = *(pa + m);
            if (m > 1 || (((((k < p-1) && (l > 0)) || (k < p-2)) && (m==1))))
            {
                a[k].na1 = *(pa + m + 1);
            }     
        }
     
        a[k].array = pa;
    }

    //Запуск потоков.
    for (k = 1; k < p; k++)
    {
        if (pthread_create(&a[k].tid, nullptr, thread_func, a+k))
        {
            std::cerr << "Error creating thread " << k <<std::endl;
            
            delete[] array;
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
    
    /*for (k = 0; k < p; k++)
    {
        a[k].PrintAll();
    }*/
    
    
    for (int i = 0; i < n; i++)
    {
        printf(" %8.2e", array[i]);
    }
    printf("\n");
        
    
    delete[] array;
    return 0;
}